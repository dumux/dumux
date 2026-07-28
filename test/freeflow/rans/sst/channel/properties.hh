// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the turbulent channel flow test using the two-equation SST
 *        RANS model.
 */
#ifndef DUMUX_TEST_RANS_SST_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_RANS_SST_CHANNEL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/rans/sst/model.hh>
#include <dumux/freeflow/rans/sst/momentumproblem.hh>
#include <dumux/freeflow/rans/sst/massproblem.hh>

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RANSSSTChannelTest {};
struct RANSSSTChannelTestMomentum { using InheritsFrom = std::tuple<RANSSSTChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if NONISOTHERMAL
struct RANSSSTChannelTestMass { using InheritsFrom = std::tuple<RANSSSTChannelTest, NavierStokesMassOneSSTNI, CCTpfaModel>; };
#else
struct RANSSSTChannelTestMass { using InheritsFrom = std::tuple<RANSSSTChannelTest, NavierStokesMassOneSST, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property: the momentum problem is layered on top of Dumux::SSTMomentumProblem
// (the SST RANS closure mixin), the mass problem on top of Dumux::SSTMassProblem (which needs
// the momentum problem's concrete type as an explicit second template parameter, see
// dumux/freeflow/rans/sst/massproblem.hh).
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSSSTChannelTestMomentum>
{ using type = RANSSSTChannelTestProblem<TypeTag, Dumux::SSTMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::RANSSSTChannelTestMass>
{
private:
    using MomentumProblem = GetPropType<TTag::RANSSSTChannelTestMomentum, Properties::Problem>;
public:
    using type = RANSSSTChannelTestProblem<TypeTag, Dumux::SSTMassProblem<TypeTag, MomentumProblem>>;
};

// the fluid system: air, matching the k-omega test
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSSSTChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type: a graded (wall-refined) tensor-product-coordinate YaspGrid, matching the
// other RANS tests (see params.input: Positions0/1, Cells0/1, Grading1).
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSSSTChannelTest>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSSSTChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSSSTChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSSSTChannelTest>
{ static constexpr bool value = true; };

// Set the coupling manager property: the existing, unmodified momentum<->mass coupling
// manager - no dedicated turbulence coupling manager is needed, see
// whatisimplemented.md/proposedimplementation.md.
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RANSSSTChannelTest>
{
    using Traits = MultiDomainTraits<TTag::RANSSSTChannelTestMomentum, TTag::RANSSSTChannelTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
