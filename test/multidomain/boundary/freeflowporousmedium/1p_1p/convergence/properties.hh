// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The properties for a simple Darcy test (cell-centered finite volume method)
 */
#ifndef DUMUX_DARCYSTOKES_PROPERTIES_HH
#define DUMUX_DARCYSTOKES_PROPERTIES_HH

#ifndef NAVIER_STOKES_MODEL
#define NAVIER_STOKES_MODEL NavierStokesMomentum
#endif

#ifndef MOMENTUM_DISCRETIZATION_MODEL
#define MOMENTUM_DISCRETIZATION_MODEL FaceCenteredStaggeredModel
#endif

#ifndef MASS_DISCRETIZATION_MODEL
#define MASS_DISCRETIZATION_MODEL CCTpfaModel
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>

#include <dumux/multidomain/boundary/freeflowporousmedium/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "spatialparams.hh"
#include "problem_darcy.hh"
#include "problem_freeflow.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DarcyOneP { using InheritsFrom = std::tuple<OneP, MASS_DISCRETIZATION_MODEL>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOneP>
{ using type = DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOneP>
{
#ifndef GRIDTYPE
    using type = Dune::YaspGrid<2>;
#else
    using type = GRIDTYPE;
#endif
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOneP>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ConvergenceTestSpatialParams<GridGeometry, Scalar>;
};

// Create new type tags
namespace TTag {
struct FreeFlowOneP {};
struct FreeFlowOnePMass { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMassOneP, MASS_DISCRETIZATION_MODEL>; };
struct FreeFlowOnePMomentum { using InheritsFrom = std::tuple<FreeFlowOneP, NAVIER_STOKES_MODEL, MOMENTUM_DISCRETIZATION_MODEL>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowOneP>
{
#ifndef GRIDTYPE
    using Coords = Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>;
    using type = Dune::YaspGrid<2, Coords>;
#else
    using type = GRIDTYPE;
#endif
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMomentum>
{ using type = FreeFlowSubProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMass>
{ using type = FreeFlowSubProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOneP>
{
    using Traits = MultiDomainTraits<TTag::FreeFlowOnePMomentum, TTag::FreeFlowOnePMass, TTag::DarcyOneP>;
    using type = FreeFlowPorousMediumCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FreeFlowOneP>
{
    using Traits = MultiDomainTraits<TTag::FreeFlowOnePMomentum, TTag::FreeFlowOnePMass, TTag::DarcyOneP>;
    using type = FreeFlowPorousMediumCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
