// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief The properties for the two-phase pore network model.
 */
#ifndef DUMUX_PNM_2P_THETA_PROPERTIES_HH
#define DUMUX_PNM_2P_THETA_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porenetwork/2p/spatialparams.hh>
#include <dumux/porenetwork/2p/invasionstate.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include <dumux/common/fvproblem.hh>

#include <dumux/common/properties.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/porenetwork/common/utilities.hh>

#include "problem_porenetwork.hh"
#include "spatialparams_porenetwork.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoP>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DrainageProblem> { using type = DrainageProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DrainageProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

#if USETHETAREGULARIZATION
//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::DrainageProblem>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
    using Traits = PoreNetwork::PNMTwoPDefaultGridFVCTraits<Problem,
                                                            FluxVariablesCache,
                                                            Dumux::PoreNetwork::TwoPInvasionState<Problem, Dumux::PoreNetwork::StateSwitchMethod::theta>>;
public:
    using type = PoreNetwork::PNMTwoPGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache, Traits>;
};
#endif

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DrainageProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalRules = PoreNetwork::FluidMatrix::MultiShapeTwoPLocalRules<Scalar>;
public:
    using type = PoreNetwork::TwoPDrainageSpatialParams<GridGeometry, Scalar, LocalRules>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DrainageProblem> { using type = Dune::FoamGrid<1, 1>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DrainageProblem>
{
public:
    using type = void;
};

} //end namespace Dumux::Properties

#endif
