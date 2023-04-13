// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The spatial params 2p rotational symmetry test
 */
#ifndef DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_SPATIAL_PARAMS_HH
#define DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class TwoPRotationalSymmetrySpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, TwoPRotationalSymmetrySpatialParams<GridGeometry, Scalar>>
{
    using ThisType = TwoPRotationalSymmetrySpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;
public:
    using PermeabilityType = Scalar;

    TwoPRotationalSymmetrySpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    ,pcKrSwCurve_("SpatialParams")
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1e-11; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }




    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSwCurve_); }

    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
