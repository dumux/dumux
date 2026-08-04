// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief The one-equation (Spalart-Allmaras) RANS turbulence closure, fused as a second
 *        transported equation (the working viscosity ν̃) onto the single-phase Navier-Stokes
 *        mass balance.
 *
 * The transport equation is fused onto the *mass* sub-model rather than living on a
 * dedicated third coupled sub-domain.
 */
#ifndef DUMUX_RANS_ONEEQ_MASS_MODEL_HH
#define DUMUX_RANS_ONEEQ_MASS_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/energy/model.hh>
#include <dumux/freeflow/rans/common/thermalconductivitymodel.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "localresidual.hh"
#include "iofields.hh"
#include "advectiveflux.hh"

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Traits for the one-equation (Spalart-Allmaras) turbulence closure fused onto the
 *        single-phase Navier-Stokes mass model.
 */
struct OneEqMassModelTraits : public NavierStokesMassOnePModelTraits
{
    //! One mass balance equation, plus one turbulence transport equation (ν̃).
    static constexpr int numEq() { return 2; }

    using Indices = OneEqMassIndices;
};

namespace Properties {

namespace TTag {
//! The type tag for the single-phase Navier-Stokes mass model fused with the one-equation
//! (Spalart-Allmaras) turbulence closure.
struct NavierStokesMassOnePOneEq { using InheritsFrom = std::tuple<NavierStokesMassOneP>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOnePOneEq>
{ using type = OneEqMassModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassOnePOneEq>
{ using type = OneEqMassLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOnePOneEq>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesMassOnePVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = OneEqMassVolumeVariables<Traits>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOnePOneEq> { using type = OneEqMassIOFields; };

//////////////////////////////////////////////////////////////////////////////////////
// Nonisothermal variant.
//////////////////////////////////////////////////////////////////////////////////////

namespace TTag {
struct NavierStokesMassOnePOneEqNI { using InheritsFrom = std::tuple<NavierStokesMassOnePOneEq>; };
} // end namespace TTag

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{ using type = NavierStokesEnergyIOFields<OneEqMassIOFields>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{ using type = NavierStokesEnergyModelTraits<OneEqMassModelTraits>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using BaseTraits = NavierStokesMassOnePVolumeVariablesTraits<PV, FSY, FST, MT>;
    using ETCM = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using HCT = GetPropType<TypeTag, Properties::HeatConductionType>;
    struct NITraits : public BaseTraits
    {
        using EffectiveThermalConductivityModel = ETCM;
        using HeatConductionType = HCT;
    };
public:
    using type = OneEqMassVolumeVariables<NITraits>;
};

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{ using type = RANSThermalConductivityModel; };

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{ using type = FouriersLaw<TypeTag>; };

template<class TypeTag>
struct FluxVariables<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    struct DiffusiveFluxTypes
    { using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>; };

    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
public:
    using type = NavierStokesMassOnePFluxVariables<
        Problem, ModelTraits, DiffusiveFluxTypes, ElementVolumeVariables, ElementFluxVariablesCache
    >;
};

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{
    struct type : public GetPropType<TypeTag, Properties::HeatConductionType>::Cache
    {};
};

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    static constexpr bool diffusionIsSolDependent = false;
    static constexpr bool heatConductionIsSolDependent
        = getPropValue<TypeTag, Properties::SolutionDependentHeatConduction>();

    using type = FreeFlowScalarFluxVariablesCacheFiller<
        Problem, ModelTraits, diffusionIsSolDependent, heatConductionIsSolDependent
    >;
};

template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::NavierStokesMassOnePOneEqNI>
{ static constexpr bool value = true; };

} // end namespace Properties
} // end namespace Dumux

#endif
