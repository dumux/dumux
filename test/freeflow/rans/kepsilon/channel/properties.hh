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
 *        wall-function k-epsilon RANS model.
 */
#ifndef DUMUX_TEST_RANS_KEPSILON_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_RANS_KEPSILON_CHANNEL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/rans/kepsilon/model.hh>
#include <dumux/freeflow/rans/kepsilon/momentumproblem.hh>
#include <dumux/freeflow/rans/kepsilon/massproblem.hh>

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RANSKEpsilonChannelTest {};
struct RANSKEpsilonChannelTestMomentum { using InheritsFrom = std::tuple<RANSKEpsilonChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct RANSKEpsilonChannelTestMass { using InheritsFrom = std::tuple<RANSKEpsilonChannelTest, NavierStokesMassOneKEpsilon, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property: the momentum problem is layered on top of Dumux::KEpsilonMomentumProblem
// (the k-epsilon RANS closure mixin), the mass problem on top of Dumux::KEpsilonMassProblem
// (which needs the momentum problem's concrete type as an explicit second template parameter,
// see dumux/freeflow/rans/kepsilon/massproblem.hh).
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSKEpsilonChannelTestMomentum>
{ using type = RANSKEpsilonChannelTestProblem<TypeTag, Dumux::KEpsilonMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::RANSKEpsilonChannelTestMass>
{
private:
    using MomentumProblem = GetPropType<TTag::RANSKEpsilonChannelTestMomentum, Properties::Problem>;
public:
    using type = RANSKEpsilonChannelTestProblem<TypeTag, Dumux::KEpsilonMassProblem<TypeTag, MomentumProblem>>;
};

// the fluid system: air, matching the other RANS tests and (before them)
// releases/3.10:test/freeflow/rans/ (PipeLauferProblem)
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSKEpsilonChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type: a graded (wall-refined) tensor-product-coordinate YaspGrid, matching the
// other RANS tests (see params.input: Positions0/1, Cells0/1, Grading1).
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSKEpsilonChannelTest>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSKEpsilonChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSKEpsilonChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSKEpsilonChannelTest>
{ static constexpr bool value = true; };

// Set the coupling manager property: the existing, unmodified momentum<->mass coupling
// manager - no dedicated turbulence coupling manager is needed, see
// whatisimplemented.md/proposedimplementation.md.
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RANSKEpsilonChannelTest>
{
    using Traits = MultiDomainTraits<TTag::RANSKEpsilonChannelTestMomentum, TTag::RANSKEpsilonChannelTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
