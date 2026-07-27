// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the turbulent channel flow test using the zero-equation RANS model.
 */
#ifndef DUMUX_TEST_RANS_ZEROEQ_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_RANS_ZEROEQ_CHANNEL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RANSZeroEqChannelTest {};
struct RANSZeroEqChannelTestMomentum { using InheritsFrom = std::tuple<RANSZeroEqChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct RANSZeroEqChannelTestMass { using InheritsFrom = std::tuple<RANSZeroEqChannelTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property: the momentum problem is layered on top of Dumux::ZeroEqProblem
// (the zero-equation RANS closure mixin) instead of directly on Dumux::NavierStokesMomentumProblem
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSZeroEqChannelTestMomentum>
{ using type = RANSZeroEqChannelTestProblem<TypeTag, Dumux::ZeroEqProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::RANSZeroEqChannelTestMass>
{ using type = RANSZeroEqChannelTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

// the fluid system: air, matching releases/3.10:test/freeflow/rans/ (PipeLauferProblem)
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSZeroEqChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type: a graded (wall-refined) tensor-product-coordinate YaspGrid, matching
// the mesh grading approach used by releases/3.10:test/freeflow/rans/ (see params.input:
// Positions0/1, Cells0/1, Grading1) - a plain equidistant grid does not resolve the near-wall
// region finely enough for a physically meaningful turbulent boundary layer.
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSZeroEqChannelTest>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSZeroEqChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSZeroEqChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSZeroEqChannelTest>
{ static constexpr bool value = true; };

// Set the coupling manager property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RANSZeroEqChannelTest>
{
    using Traits = MultiDomainTraits<TTag::RANSZeroEqChannelTestMomentum, TTag::RANSZeroEqChannelTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
