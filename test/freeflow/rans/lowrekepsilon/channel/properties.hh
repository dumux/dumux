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
#if NONISOTHERMAL
struct RANSLowReKEpsilonChannelTestMass { using InheritsFrom = std::tuple<RANSLowReKEpsilonChannelTest, NavierStokesMassOneLowReKEpsilonNI, CCTpfaModel>; };
#else
struct RANSLowReKEpsilonChannelTestMass { using InheritsFrom = std::tuple<RANSLowReKEpsilonChannelTest, NavierStokesMassOneLowReKEpsilon, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property: the momentum problem is based on Dumux::LowReKEpsilonMomentumProblem and the
// mass problem on Dumux::LowReKEpsilonMassProblem, which takes the type of the momentum problem as
// second template parameter (see dumux/freeflow/rans/lowrekepsilon/massproblem.hh)
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

// the fluid system: air
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSLowReKEpsilonChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type: a YaspGrid with tensor-product coordinates, which allows grading
// the mesh towards the walls (see params.input: Positions, Cells, Grading)
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

// Set the coupling manager property: k and epsilon-tilde are solved as additional equations of the
// mass sub-model, so the standard momentum-mass coupling manager is sufficient
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RANSLowReKEpsilonChannelTest>
{
    using Traits = MultiDomainTraits<TTag::RANSLowReKEpsilonChannelTestMomentum, TTag::RANSLowReKEpsilonChannelTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
