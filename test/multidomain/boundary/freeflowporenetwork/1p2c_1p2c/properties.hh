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
#ifndef DUMUX_FREEFLOW_PORE_NETWORK_ONEPTWOC_PROPERTIES_HH
#define DUMUX_FREEFLOW_PORE_NETWORK_ONEPTWOC_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/porenetwork/1pnc/model.hh>
#include <dumux/freeflow/navierstokes/mass/1pnc/model.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>

#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/multidomain/boundary/freeflowporenetwork/couplingmanager.hh>

#include "problem_pnm.hh"
#include "problem_freeflow.hh"

namespace Dumux::Properties {

//////////////////////////////////// PNM ////////////////////////////////////
// Create new type tags
namespace TTag {
#if ISOTHERMAL
struct PNMOnePNCModel { using InheritsFrom = std::tuple<PNMOnePNC>; };
#else
struct PNMOnePNCModel { using InheritsFrom = std::tuple<PNMOnePNCNI>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePNCModel> { using type = Dumux::PNMOnePNCProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePNCModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2OAir = FluidSystems::H2OAir<Scalar>;
    static constexpr int phaseIdx = H2OAir::gasPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePNCModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = Dumux::PoreNetwork::TransmissibilityBruus<Scalar>;
public:
    using type =  Dumux::PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePNCModel> { using type = Dune::FoamGrid<1, 2>; };

// default usesMoles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::PNMOnePNCModel> { static constexpr bool value = true; };

//! Set as default that no component mass balance is replaced by the total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::PNMOnePNCModel>
{
    static constexpr auto value = 3;
};

//////////////////////////////////// FF ////////////////////////////////////
// Create new type tags
namespace TTag {
struct FreeFlowOnePNC {};
struct FreeFlowOnePNCMomentum { using InheritsFrom = std::tuple<FreeFlowOnePNC, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if ISOTHERMAL
struct FreeFlowOnePNCMass { using InheritsFrom = std::tuple<FreeFlowOnePNC, NavierStokesMassOnePNC, CCTpfaModel>; };
#else
struct FreeFlowOnePNCMass { using InheritsFrom = std::tuple<FreeFlowOnePNC, NavierStokesMassOnePNCNI, CCTpfaModel>; };
#endif
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowOnePNC>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2OAir = FluidSystems::H2OAir<Scalar>;
    static constexpr int phaseIdx = H2OAir::gasPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowOnePNC>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, 2>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePNCMomentum>
{ using type = Dumux::FreeFlowOnePNCTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePNCMass>
{ using type = Dumux::FreeFlowOnePNCTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FreeFlowOnePNC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowOnePNC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowOnePNC> { static constexpr bool value = true; };
//Use mole fraction formulation
template<class TypeTag>
struct UseMoles<TypeTag, TTag::FreeFlowOnePNC> { static constexpr bool value = true; };

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::FreeFlowOnePNC>
{ static constexpr int value = 3; };

////////////////////////////////CouplingManager///////////////////////////////
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PNMOnePNCModel>
{
    using Traits = MultiDomainTraits<Properties::TTag::FreeFlowOnePNCMomentum, Properties::TTag::FreeFlowOnePNCMass, TTag::PNMOnePNCModel>;
    using type = Dumux::FreeFlowPoreNetworkCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FreeFlowOnePNC>
{
    using Traits = MultiDomainTraits<Properties::TTag::FreeFlowOnePNCMomentum, Properties::TTag::FreeFlowOnePNCMass, TTag::PNMOnePNCModel>;
    using type = Dumux::FreeFlowPoreNetworkCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
