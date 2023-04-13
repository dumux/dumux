// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the matrix and fracture problem.
 */

#ifndef DUMUX_FRACTURE_TEST_SPATIAL_PARAMS_HH
#define DUMUX_FRACTURE_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the matrix and fracture problem.
 */
template<class GridGeometry, class Scalar>
class MatrixFractureSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, MatrixFractureSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = MatrixFractureSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    MatrixFractureSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                                const std::string& paramGroup = "")
    : ParentType(gridGeometry)
    {
        permeability_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability");
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity", 1.0);
        extrusion_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture", 1.0);
    }

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param element The current finite element
     * \param scv The sub control volume
     * \param elemSol The current element solution vector
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosity_;
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
    { return extrusion_; }

private:
    Scalar permeability_;
    Scalar porosity_;
    Scalar extrusion_;
};

} // end namespace Dumux

#endif
