// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the blood flow problem.
 */

#ifndef DUMUX_BlOOD_FLOW_SPATIALPARAMS_HH
#define DUMUX_BlOOD_FLOW_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the blood flow problem.
 */
template<class GridGeometry, class Scalar>
class BloodFlowSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, BloodFlowSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = BloodFlowSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    BloodFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        radius_ = getParam<Scalar>("SpatialParams.Radius");
    }

    /*!
     * \brief Returns the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \param ipGlobal The integration point
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& ipGlobal) const
    {
        return (1 + ipGlobal[2] + 0.5*ipGlobal[2]*ipGlobal[2])/(M_PI*radius(0)*radius(0));
    }

    //! we evaluate the permeability directly at the scvf since we have an analytical expression for it
    static constexpr bool evaluatePermeabilityAtScvfIP()
    { return true; }

    /*!
     * \brief Returns the radius of the circular pipe for the current sub-control volume in [m].
     *
     * \param eIdxGlobal the index of the element
     */
    Scalar radius(unsigned int eIdxGlobal) const
    {
        return radius_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$.
     *
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 1.0;
    }

    /*!
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto r = radius(eIdx);
        return M_PI*r*r;
    }

private:
    Scalar radius_;
};

} // end namespace Dumux

#endif
