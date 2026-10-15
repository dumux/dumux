// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Spatial parameters of the 2p1cni heat pipe test.
 */
#ifndef DUMUX_TEST_2P1C_HEATPIPE_SPATIALPARAMS_HH
#define DUMUX_TEST_2P1C_HEATPIPE_SPATIALPARAMS_HH

#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/heatpipelaw.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCTests
 * \brief Homogeneous medium with Leverett's capillary pressure and cubic relative permeabilities.
 *
 * The cubic relative permeabilities of the HeatPipeLaw hold for arguments up to 0.95,
 * which covers the saturation range of the steady state.
 */
template<class GridGeometry, class Scalar>
class TwoPOneCHeatPipeSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, TwoPOneCHeatPipeSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = TwoPOneCHeatPipeSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::HeatPipeLaw<Scalar>;

public:
    using PermeabilityType = Scalar;

    TwoPOneCHeatPipeSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    , pcKrSwCurve_(typename PcKrSwCurve::Params(getParam<Scalar>("SpatialParams.SurfaceTension"),
                                                std::sqrt(porosity_/permeability_)))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSwCurve_); }

    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::liquidPhaseIdx; }

private:
    Scalar permeability_;
    Scalar porosity_;
    PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
