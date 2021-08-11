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
 * \ingroup TwoPTests
 * \brief spatial params for the flow problem of the tracer conservation test
 */

#ifndef DUMUX_TWOP_FLOW_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TWOP_FLOW_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief spatial params for the flow problem of the tracer conservation test
 */
template<class GridGeometry, class Scalar>
class TwoPFlowTestSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, TwoPFlowTestSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = TwoPFlowTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PcKrSwCurve = FluidMatrix::LinearMaterialDefault<Scalar>;
    using PermeabilityType = Scalar;

    TwoPFlowTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSw_("SpatialParams")
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1e-10; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        return makeFluidMatrixInteraction(pcKrSw_);
    }

    auto fluidMatrixInteraction(const Element& element) const
    {
        return makeFluidMatrixInteraction(pcKrSw_);
    }

    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    const PcKrSwCurve pcKrSw_;
};

} // end namespace Dumux

#endif
