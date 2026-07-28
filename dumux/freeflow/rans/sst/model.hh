// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief The two-equation SST (Shear Stress Transport, Menter 1994) RANS turbulence closure,
 *        fused as two extra transported equations (turbulent kinetic energy k, specific
 *        dissipation rate omega) onto the single-phase Navier-Stokes mass balance.
 *
 * SST blends the k-omega closure (near the wall) with a transformed k-epsilon closure (in the
 * free stream) via the F1 blending function, and optionally limits the eddy viscosity using
 * Bradshaw's assumption (the runtime-selectable RANS.SSTModelVersion = "SST" vs. "BSL" variant,
 * see Dumux::SSTModel in dumux/freeflow/turbulencemodel.hh). See turbulenceequations.md for the
 * physics, and whatisimplemented.md/proposedimplementation.md for why these equations live on
 * the *mass* sub-model, following the same strategy already validated for k-omega.
 */
#ifndef DUMUX_RANS_SST_MASS_MODEL_HH
#define DUMUX_RANS_SST_MASS_MODEL_HH

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
 * \brief Traits for the two-equation SST turbulence closure fused onto the single-phase
 *        Navier-Stokes mass model.
 */
struct SSTMassModelTraits : public NavierStokesMassOnePModelTraits
{
    //! One mass balance equation, plus two turbulence transport equations (k, omega).
    static constexpr int numEq() { return 3; }

    using Indices = SSTMassIndices;
};

namespace Properties {

namespace TTag {
//! The type tag for the single-phase Navier-Stokes mass model fused with the two-equation
//! SST turbulence closure.
struct NavierStokesMassOneSST { using InheritsFrom = std::tuple<NavierStokesMassOneP>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOneSST>
{ using type = SSTMassModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassOneSST>
{ using type = SSTMassLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOneSST>
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
    using type = SSTMassVolumeVariables<Traits>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOneSST> { using type = SSTMassIOFields; };

//////////////////////////////////////////////////////////////////////////////////////
// Nonisothermal variant - see whatisimplemented.md's Phase 8 section.
//////////////////////////////////////////////////////////////////////////////////////

namespace TTag {
struct NavierStokesMassOneSSTNI { using InheritsFrom = std::tuple<NavierStokesMassOneSST>; };
} // end namespace TTag

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOneSSTNI>
{ using type = NavierStokesEnergyIOFields<SSTMassIOFields>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOneSSTNI>
{ using type = NavierStokesEnergyModelTraits<SSTMassModelTraits>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOneSSTNI>
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
    using type = SSTMassVolumeVariables<NITraits>;
};

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::NavierStokesMassOneSSTNI>
{ using type = RANSThermalConductivityModel; };

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::NavierStokesMassOneSSTNI>
{ using type = FouriersLaw<TypeTag>; };

template<class TypeTag>
struct FluxVariables<TypeTag, TTag::NavierStokesMassOneSSTNI>
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
struct FluxVariablesCache<TypeTag, TTag::NavierStokesMassOneSSTNI>
{
    struct type : public GetPropType<TypeTag, Properties::HeatConductionType>::Cache
    {};
};

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::NavierStokesMassOneSSTNI>
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
struct SolutionDependentHeatConduction<TypeTag, TTag::NavierStokesMassOneSSTNI>
{ static constexpr bool value = true; };

} // end namespace Properties
} // end namespace Dumux

#endif
