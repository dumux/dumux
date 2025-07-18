// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Pore Scale Simulations
 * \ingroup Basic Pore Structure
 * \brief The properties file for the evaluation of a pore geometry
 */

#include <dune/grid/spgrid.hh>
#include <dune/common/hybridutilities.hh>

#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/fcstaggered.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

#ifndef DUMUX_PERIODIC_PORE_STRUCTURE_PROPERTIES_HH
#define DUMUX_PERIODIC_PORE_STRUCTURE_PROPERTIES_HH

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct FlowPoreScaleModel {};
struct FlowMassModel { using InheritsFrom = std::tuple<FlowPoreScaleModel, NavierStokesMassOneP, CCTpfaModel>; };
struct FlowMomentumModel { using InheritsFrom = std::tuple<FlowPoreScaleModel, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
} // end namespace TTag

///////////////////////
/// FLOW Properties ///
///////////////////////

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FlowPoreScaleModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FlowPoreScaleModel>
{
    static constexpr auto dim = 2;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using HostGrid = Dune::SPGrid<Scalar, dim>;
    using type = Dune::SubGrid<dim, HostGrid>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FlowMomentumModel>
{ using type = FlowPoreStructureProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FlowMassModel>
{ using type = FlowPoreStructureProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FlowPoreScaleModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FlowPoreScaleModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FlowPoreScaleModel> { static constexpr bool value = true; };

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FlowPoreScaleModel>
{
    using Traits = MultiDomainTraits<TTag::FlowMomentumModel, TTag::FlowMassModel>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Properties

#endif
