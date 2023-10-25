// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesNCTests
 * \brief The properties of the test for the compositional staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_DENSITY_FLOW_NC_TEST_PROPERTIES_HH
#define DUMUX_DENSITY_FLOW_NC_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/mass/1pnc/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DensityDrivenFlowTest {};
struct DensityDrivenFlowMomentum { using InheritsFrom = std::tuple<DensityDrivenFlowTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if !NONISOTHERMAL
struct DensityDrivenFlowMass { using InheritsFrom = std::tuple<DensityDrivenFlowTest, NavierStokesMassOnePNC, CCTpfaModel>; };
#else
struct DensityDrivenFlowMass { using InheritsFrom = std::tuple<DensityDrivenFlowTest, NavierStokesMassOnePNCNI, CCTpfaModel>; };
#endif
} // end namespace TTag

// Select the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DensityDrivenFlowTest>
{
    using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    static constexpr int phaseIdx = H2OAir::liquidPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DensityDrivenFlowTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DensityDrivenFlowMomentum>
{ using type = DensityDrivenFlowProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DensityDrivenFlowMass>
{ using type = DensityDrivenFlowProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr bool value = true; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr bool value = true; };

// Set the coupling manager
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DensityDrivenFlowTest>
{
private:
    using Traits = MultiDomainTraits<TTag::DensityDrivenFlowMomentum, TTag::DensityDrivenFlowMass>;
public:
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
