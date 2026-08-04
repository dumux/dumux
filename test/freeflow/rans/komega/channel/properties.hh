// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the turbulent channel flow test using the two-equation k-omega
 *        RANS model.
 */
#ifndef DUMUX_TEST_RANS_KOMEGA_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_RANS_KOMEGA_CHANNEL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/rans/komega/model.hh>
#include <dumux/freeflow/rans/komega/momentumproblem.hh>
#include <dumux/freeflow/rans/komega/massproblem.hh>

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RANSKOmegaChannelTest {};
struct RANSKOmegaChannelTestMomentum { using InheritsFrom = std::tuple<RANSKOmegaChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if NONISOTHERMAL
struct RANSKOmegaChannelTestMass { using InheritsFrom = std::tuple<RANSKOmegaChannelTest, NavierStokesMassOneKOmegaNI, CCTpfaModel>; };
#else
struct RANSKOmegaChannelTestMass { using InheritsFrom = std::tuple<RANSKOmegaChannelTest, NavierStokesMassOneKOmega, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property: the momentum problem is layered on top of Dumux::KOmegaMomentumProblem
// (the k-omega RANS closure mixin), the mass problem on top of Dumux::KOmegaMassProblem
// (which needs the momentum problem's concrete type as an explicit second template parameter,
// see dumux/freeflow/rans/komega/massproblem.hh).
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSKOmegaChannelTestMomentum>
{ using type = RANSKOmegaChannelTestProblem<TypeTag, Dumux::KOmegaMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::RANSKOmegaChannelTestMass>
{
private:
    using MomentumProblem = GetPropType<TTag::RANSKOmegaChannelTestMomentum, Properties::Problem>;
public:
    using type = RANSKOmegaChannelTestProblem<TypeTag, Dumux::KOmegaMassProblem<TypeTag, MomentumProblem>>;
};

// the fluid system: air, matching the one-equation test and (before it)
// releases/3.10:test/freeflow/rans/ (PipeLauferProblem)
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSKOmegaChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type: a graded (wall-refined) tensor-product-coordinate YaspGrid, matching the
// one-equation/zero-equation tests (see params.input: Positions0/1, Cells0/1, Grading1).
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSKOmegaChannelTest>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSKOmegaChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSKOmegaChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSKOmegaChannelTest>
{ static constexpr bool value = true; };

// Set the coupling manager property: the existing, unmodified momentum<->mass coupling
// manager - no dedicated turbulence coupling manager is needed, since k and omega are
// solved as extra equations on the mass sub-model rather than on a separate sub-domain.
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RANSKOmegaChannelTest>
{
    using Traits = MultiDomainTraits<TTag::RANSKOmegaChannelTestMomentum, TTag::RANSKOmegaChannelTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
