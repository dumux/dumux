// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the turbulent channel flow test using the two-equation,
 *        low-Reynolds-number k-epsilon RANS model.
 */
#ifndef DUMUX_TEST_RANS_LOWREKEPSILON_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_RANS_LOWREKEPSILON_CHANNEL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/rans/lowrekepsilon/model.hh>
#include <dumux/freeflow/rans/lowrekepsilon/momentumproblem.hh>
#include <dumux/freeflow/rans/lowrekepsilon/massproblem.hh>

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RANSLowReKEpsilonChannelTest {};
struct RANSLowReKEpsilonChannelTestMomentum { using InheritsFrom = std::tuple<RANSLowReKEpsilonChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct RANSLowReKEpsilonChannelTestMass { using InheritsFrom = std::tuple<RANSLowReKEpsilonChannelTest, NavierStokesMassOneLowReKEpsilon, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property: the momentum problem is layered on top of
// Dumux::LowReKEpsilonMomentumProblem (the RANS closure mixin), the mass problem on top of
// Dumux::LowReKEpsilonMassProblem (which needs the momentum problem's concrete type as an
// explicit second template parameter, see dumux/freeflow/rans/lowrekepsilon/massproblem.hh).
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSLowReKEpsilonChannelTestMomentum>
{ using type = RANSLowReKEpsilonChannelTestProblem<TypeTag, Dumux::LowReKEpsilonMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::RANSLowReKEpsilonChannelTestMass>
{
private:
    using MomentumProblem = GetPropType<TTag::RANSLowReKEpsilonChannelTestMomentum, Properties::Problem>;
public:
    using type = RANSLowReKEpsilonChannelTestProblem<TypeTag, Dumux::LowReKEpsilonMassProblem<TypeTag, MomentumProblem>>;
};

// the fluid system: air, matching the k-omega/one-equation tests
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSLowReKEpsilonChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type: a graded (wall-refined) tensor-product-coordinate YaspGrid, matching the
// other RANS tests (see params.input: Positions0/1, Cells0/1, Grading1).
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSLowReKEpsilonChannelTest>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSLowReKEpsilonChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSLowReKEpsilonChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSLowReKEpsilonChannelTest>
{ static constexpr bool value = true; };

// Set the coupling manager property: the existing, unmodified momentum<->mass coupling
// manager - no dedicated turbulence coupling manager is needed, see
// whatisimplemented.md/proposedimplementation.md.
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RANSLowReKEpsilonChannelTest>
{
    using Traits = MultiDomainTraits<TTag::RANSLowReKEpsilonChannelTestMomentum, TTag::RANSLowReKEpsilonChannelTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
