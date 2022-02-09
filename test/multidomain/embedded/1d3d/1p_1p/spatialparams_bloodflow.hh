// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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
