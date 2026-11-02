// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the turbulent channel flow test using the one-equation
 *        (Spalart-Allmaras) RANS model.
 */
#ifndef DUMUX_TEST_RANS_ONEEQ_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_RANS_ONEEQ_CHANNEL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/rans/oneeq/model.hh>
#include <dumux/freeflow/rans/oneeq/momentumproblem.hh>
#include <dumux/freeflow/rans/oneeq/massproblem.hh>

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RANSOneEqChannelTest {};
struct RANSOneEqChannelTestMomentum { using InheritsFrom = std::tuple<RANSOneEqChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if NONISOTHERMAL
struct RANSOneEqChannelTestMass { using InheritsFrom = std::tuple<RANSOneEqChannelTest, NavierStokesMassOnePOneEqNI, CCTpfaModel>; };
#else
struct RANSOneEqChannelTestMass { using InheritsFrom = std::tuple<RANSOneEqChannelTest, NavierStokesMassOnePOneEq, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property: the momentum problem is based on Dumux::OneEqMomentumProblem and the
// mass problem on Dumux::RANSMassOneEqProblem, which takes the type of the momentum problem as
// second template parameter (see dumux/freeflow/rans/oneeq/massproblem.hh)
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSOneEqChannelTestMomentum>
{ using type = RANSOneEqChannelTestProblem<TypeTag, Dumux::OneEqMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::RANSOneEqChannelTestMass>
{
private:
    using MomentumProblem = GetPropType<TTag::RANSOneEqChannelTestMomentum, Properties::Problem>;
public:
    using type = RANSOneEqChannelTestProblem<TypeTag, Dumux::RANSMassOneEqProblem<TypeTag, MomentumProblem>>;
};

// the fluid system: air
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSOneEqChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type: a YaspGrid with tensor-product coordinates, which allows grading
// the mesh towards the walls (see params.input: Positions, Cells, Grading)
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSOneEqChannelTest>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSOneEqChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSOneEqChannelTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSOneEqChannelTest>
{ static constexpr bool value = true; };

// Set the coupling manager property: the working viscosity is solved as additional equations of the
// mass sub-model, so the standard momentum-mass coupling manager is sufficient
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RANSOneEqChannelTest>
{
    using Traits = MultiDomainTraits<TTag::RANSOneEqChannelTestMomentum, TTag::RANSOneEqChannelTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
