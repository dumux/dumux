// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MPNCTests
 * \brief The properties of the problem where hot, pure liquid water is injected from the left hand
 * side into a initially isotherm domain.
 */
#ifndef DUMUX_COMBUSTION_PROPERTIES_ONE_COMPONENT_HH
#define DUMUX_COMBUSTION_PROPERTIES_ONE_COMPONENT_HH

#include <dune/grid/onedgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/mpnc/model.hh>
#include <dumux/porousmediumflow/mpnc/pressureformulation.hh>

#include <dumux/material/solidstates/compositionalsolidstate.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/simplefluidlumping.hh>

#include "spatialparams.hh"
#include "combustionfluidsystem.hh"
#include "combustionlocalresidual.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct CombustionOneComponent { using InheritsFrom = std::tuple<MPNCNonequil>; };
struct CombustionOneComponentBox { using InheritsFrom = std::tuple<CombustionOneComponent, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::CombustionOneComponent> { using type = Dune::OneDGrid; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::CombustionOneComponent>
{ using type = CombustionProblemOneComponent<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CombustionOneComponent>
{
    using GridGeometry = GetPropType<TypeTag, GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = CombustionSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::CombustionOneComponent>
{ using type = FluidSystems::CombustionFluidsystem<GetPropType<TypeTag, Properties::Scalar>>; };

//! Set the default pressure formulation: either pw first or pn first
template<class TypeTag>
struct PressureFormulation<TypeTag, TTag::CombustionOneComponent>
{
public:
    static const MpNcPressureFormulation value = MpNcPressureFormulation::mostWettingFirst;
};

//! Custom model traits to deactivate diffusion for this test
template<int numP, int numC, MpNcPressureFormulation formulation, bool useM>
struct CombustionModelTraits : public MPNCModelTraits<numP, numC, formulation, useM>
{
    static constexpr bool enableMolecularDiffusion() { return false; }
};

// We use different model traits for the equilibrium part because we want to deactivate diffusion
template<class TypeTag>
struct EquilibriumModelTraits<TypeTag, TTag::CombustionOneComponent>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = CombustionModelTraits< FluidSystem::numPhases,
                                        FluidSystem::numComponents,
                                        getPropValue<TypeTag, Properties::PressureFormulation>(),
                                        getPropValue<TypeTag, Properties::UseMoles>() >;
};

template<class TypeTag>
struct FluidState<TypeTag, TTag::CombustionOneComponent>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};
//#################
//changes from the default settings which also assume chemical non-equilibrium
//set the number of energyequations we want to use
template<class TypeTag>
struct NumEnergyEqFluid<TypeTag, TTag::CombustionOneComponent> { static constexpr int value = 1; };
template<class TypeTag>
struct NumEnergyEqSolid<TypeTag, TTag::CombustionOneComponent> { static constexpr int value = 1; };

// by default chemical non equilibrium is enabled in the nonequil model, switch that off here
template<class TypeTag>
struct EnableChemicalNonEquilibrium<TypeTag, TTag::CombustionOneComponent> { static constexpr bool value = false; };
//#################

template<class TypeTag>
struct SolidSystem<TypeTag, TTag::CombustionOneComponent>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentOne = Dumux::Components::Constant<1, Scalar>;
    using ComponentTwo = Dumux::Components::Constant<2, Scalar>;
    static constexpr int numInertComponents = 2;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo, numInertComponents>;
};

template<class TypeTag>
struct SolidState<TypeTag, TTag::CombustionOneComponent>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = CompositionalSolidState<Scalar, SolidSystem>;
};

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::CombustionOneComponent>
{ using type = CombustionEnergyLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties

#endif
