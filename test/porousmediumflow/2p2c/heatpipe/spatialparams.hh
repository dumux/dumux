// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux-Lecture contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HEATPIPE_SPATIAL_PARAMS_HH
#define DUMUX_HEATPIPE_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>

#include "krpcheatpipe.hh"

namespace Dumux
{

template<class FVGridGeometry, class Scalar>
class HeatPipeSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<FVGridGeometry, Scalar, HeatPipeSpatialParams<FVGridGeometry, Scalar>>

{
    using ThisType = HeatPipeSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<FVGridGeometry, Scalar, ThisType>;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::KrPcHeatPipeDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    HeatPipeSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        permeability_ = getParam<Scalar>("Problem.Permeability");
        porosity_ = 0.4;
        Scalar p0 = std::pow((porosity_/permeability_), 0.5);

        typename PcKrSwCurve::BasicParams params(0.15, 0.0, p0);
        pcKrSwCurve_ = std::make_unique<PcKrSwCurve>(params);
    }

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(*pcKrSwCurve_); }


    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    PermeabilityType permeability_;
    Scalar porosity_;
    std::unique_ptr<const PcKrSwCurve> pcKrSwCurve_;
};

}

#endif
