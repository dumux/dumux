// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Spatial parameters for the Buckley-Leverett two-phase test.
 */
#ifndef DUMUX_TEST_TWOP_BUCKLEYLEVERETT_SPATIALPARAMS_HH
#define DUMUX_TEST_TWOP_BUCKLEYLEVERETT_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class BuckleyLeverettSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, BuckleyLeverettSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = BuckleyLeverettSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    BuckleyLeverettSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    , temperature_(getParam<Scalar>("SpatialParams.Temperature"))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSwCurve_); }

    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

    Scalar constantPorosity() const
    { return porosity_; }

    /*!
     * \brief Defines the temperature \f$[K]\f$ at the given position
     * \param globalPos The global position
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return temperature_; }

private:
    const PcKrSwCurve pcKrSwCurve_;
    Scalar permeability_;
    Scalar porosity_;
    const Scalar temperature_;
};

} // end namespace Dumux

#endif
