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
#ifndef DUMUX_PNMSTOKES_PROPERTIES_HH
#define DUMUX_PNMSTOKES_PROPERTIES_HH

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
#include "../spatialparams_freeflow.hh" //TODO: Do we need this or are the default spatial params enough?

namespace Dumux::Properties {

/////////////////////////// PNM ///////////////////////////
// Create new type tags
namespace TTag {
struct PNMOnePNCNIModel { using InheritsFrom = std::tuple<PNMOnePNCNI>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePNCNIModel> { using type = Dumux::PNMOnePNCNIProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePNCNIModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    static constexpr int phaseIdx = H2OAir::liquidPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePNCNIModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = Dumux::PoreNetwork::TransmissibilityBruus<GetPropType<TypeTag, Properties::Scalar>>;
public:
    using type =  Dumux::PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePNCNIModel> { using type = Dune::FoamGrid<1, 2>; };

// default usesMoles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::PNMOnePNCNIModel> { static constexpr bool value = true; };

//! Set as default that no component mass balance is replaced by the total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::PNMOnePNCNIModel>  //TODO: Check how to set this
{
    static constexpr auto value = 3;
};

/////////////////////////// FreeFlow ///////////////////////////
// Create new type tags
namespace TTag {
struct FreeFlowOnePNCNI {};
struct FreeFlowOnePNCNIMass { using InheritsFrom = std::tuple<FreeFlowOnePNCNI, NavierStokesMassOnePNCNI, CCTpfaModel>; }; // non-isothermal
struct FreeFlowOnePNCNIMomentum { using InheritsFrom = std::tuple<FreeFlowOnePNCNI, NavierStokesMomentum, FaceCenteredStaggeredModel>; }; // NONIsothermal not needed here as only on mass grid
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowOnePNCNI>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    static constexpr int phaseIdx = H2OAir::liquidPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowOnePNCNI>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, 2>>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FreeFlowOnePNCNI>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ChannelSpatialParams<GridGeometry, Scalar>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePNCNIMomentum>
{ using type = Dumux::FreeFlowOnePNCNITestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePNCNIMass>
{ using type = Dumux::FreeFlowOnePNCNITestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FreeFlowOnePNCNI> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowOnePNCNI> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowOnePNCNI> { static constexpr bool value = true; };

//Use mole fraction formulation
template<class TypeTag>
struct UseMoles<TypeTag, TTag::FreeFlowOnePNCNI> { static constexpr bool value = true; };

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::FreeFlowOnePNCNI> //TODO: Check how to set this
{ static constexpr int value = 3; };

////////////////////////////////CouplingManager///////////////////////////////
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PNMOnePNCNIModel>
{
    using Traits = MultiDomainTraits<Properties::TTag::FreeFlowOnePNCNIMomentum, Properties::TTag::FreeFlowOnePNCNIMass, TTag::PNMOnePNCNIModel>;
    using type = Dumux::FreeFlowPoreNetworkCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FreeFlowOnePNCNI>
{
    using Traits = MultiDomainTraits<Properties::TTag::FreeFlowOnePNCNIMomentum, Properties::TTag::FreeFlowOnePNCNIMass, TTag::PNMOnePNCNIModel>;
    using type = Dumux::FreeFlowPoreNetworkCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
