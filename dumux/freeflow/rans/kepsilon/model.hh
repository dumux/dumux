// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief The two-equation, high-Reynolds-number wall-function k-epsilon RANS turbulence
 *        closure, fused as two extra transported equations (turbulent kinetic energy k,
 *        dissipation rate epsilon) onto the single-phase Navier-Stokes mass balance.
 *
 * The two transport equations are fused onto the *mass* sub-model, following the same
 * strategy already used for the one-equation and k-omega models.
 */
#ifndef DUMUX_RANS_KEPSILON_MASS_MODEL_HH
#define DUMUX_RANS_KEPSILON_MASS_MODEL_HH

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
 * \brief Traits for the k-epsilon turbulence closure fused onto the single-phase Navier-Stokes
 *        mass model.
 */
struct KEpsilonMassModelTraits : public NavierStokesMassOnePModelTraits
{
    //! One mass balance equation, plus two turbulence transport equations (k, epsilon).
    static constexpr int numEq() { return 3; }

    using Indices = KEpsilonMassIndices;
};

namespace Properties {

namespace TTag {
//! The type tag for the single-phase Navier-Stokes mass model fused with the two-equation
//! k-epsilon turbulence closure.
struct NavierStokesMassOneKEpsilon { using InheritsFrom = std::tuple<NavierStokesMassOneP>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOneKEpsilon>
{ using type = KEpsilonMassModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassOneKEpsilon>
{ using type = KEpsilonMassLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOneKEpsilon>
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
    using type = KEpsilonMassVolumeVariables<Traits>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOneKEpsilon> { using type = KEpsilonMassIOFields; };

//////////////////////////////////////////////////////////////////////////////////////
// Nonisothermal variant. As with every other RANS model, the wall temperature is enforced
// via the same simple weak-Neumann treatment already used for k (not a Jayatilleke
// wall-function energy flux) - a deliberate simplification, not an oversight.
//////////////////////////////////////////////////////////////////////////////////////

namespace TTag {
//! The type tag for the single-phase, non-isothermal Navier-Stokes mass model fused with the
//! two-equation, wall-function k-epsilon turbulence closure.
struct NavierStokesMassOneKEpsilonNI { using InheritsFrom = std::tuple<NavierStokesMassOneKEpsilon>; };
} // end namespace TTag

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
{ using type = NavierStokesEnergyIOFields<KEpsilonMassIOFields>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
{ using type = NavierStokesEnergyModelTraits<KEpsilonMassModelTraits>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
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
    using type = KEpsilonMassVolumeVariables<NITraits>;
};

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
{ using type = RANSThermalConductivityModel; };

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
{ using type = FouriersLaw<TypeTag>; };

template<class TypeTag>
struct FluxVariables<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
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
struct FluxVariablesCache<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
{
    struct type : public GetPropType<TypeTag, Properties::HeatConductionType>::Cache
    {};
};

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
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
struct SolutionDependentHeatConduction<TypeTag, TTag::NavierStokesMassOneKEpsilonNI>
{ static constexpr bool value = true; };

} // end namespace Properties
} // end namespace Dumux

#endif
