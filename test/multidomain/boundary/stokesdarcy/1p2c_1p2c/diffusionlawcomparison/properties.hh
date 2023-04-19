// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief A simple Darcy test problem (cell-centered finite volume method) for
 *        the comparison of different diffusion laws.
 */

#ifndef DUMUX_STOKESDARCY_PROPERTIES_HH
#define DUMUX_STOKESDARCY_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/flux/maxwellstefanslaw.hh>

#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingdata.hh>
#include <dumux/multidomain/boundary/stokesdarcy/couplingmanager.hh>
#include <dumux/multidomain/staggeredtraits.hh>

#include "problem_darcy.hh"
#include "problem_stokes.hh"
#include "./../spatialparams.hh"

#ifndef DIFFUSIONTYPE
#define DIFFUSIONTYPE FicksLaw<TypeTag>
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DarcyOnePTwoC { using InheritsFrom = std::tuple<OnePNC, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePTwoC> { using type = Dumux::DarcySubProblem<TypeTag>; };

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOnePTwoC>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::liquidPhaseIdx; // simulate the water phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Use moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::DarcyOnePTwoC> { static constexpr bool value = true; };

// Do not replace one equation with a total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DarcyOnePTwoC> { static constexpr int value = 3; };

//! Use a model with constant tortuosity for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::DarcyOnePTwoC>
{ using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOnePTwoC> { using type = Dune::YaspGrid<2>; };

// Set the diffusion type
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::DarcyOnePTwoC> { using type = DIFFUSIONTYPE; };

// Set the spatial parameters type
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOnePTwoC>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePSpatialParams<GridGeometry, Scalar>;
};

// Create new type tags
namespace TTag {
struct StokesOnePTwoC { using InheritsFrom = std::tuple<NavierStokesNC, StaggeredFreeFlowModel>; };
} // end namespace TTag

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOnePTwoC>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::liquidPhaseIdx; // simulate the water phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOnePTwoC>
{ using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOnePTwoC> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };

// Use moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::StokesOnePTwoC> { static constexpr bool value = true; };

// Set the grid type
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::StokesOnePTwoC> { using type = DIFFUSIONTYPE; };

// Do not replace one equation with a total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::StokesOnePTwoC> { static constexpr int value = 3; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::StokesOnePTwoC>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyOnePTwoC>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOnePTwoC>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::StokesOnePTwoC, Properties::TTag::StokesOnePTwoC, TypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
