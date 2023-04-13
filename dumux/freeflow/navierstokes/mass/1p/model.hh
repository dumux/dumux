// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 *
 * \brief A single-phase, isothermal Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}}) = \nabla \cdot (\mu (\nabla \textbf{v} + \nabla \textbf{v}^{\text{T}}))
   - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * By setting the runtime parameter <code>Problem.EnableInertiaTerms</code> to <code>false</code> the Stokes
 * equation can be solved. In this case the term
 * \f[
 *    \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}})
 * \f]
 * is neglected.
 *
 * The <B> mass balance equation </B>
 * \f[
       \frac{\partial \varrho}{\partial t} + \nabla \cdot (\varrho \textbf{v}) - q = 0
 * \f]
 *
 * closes the system.
 *
 */

#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_1P_MODEL_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_1P_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/freeflow/spatialparams.hh>
#include <dumux/freeflow/navierstokes/iofields.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/navierstokes/energy/model.hh>
#include <dumux/freeflow/navierstokes/scalarfluxvariablescachefiller.hh>

#include <dumux/flux/fourierslaw.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "indices.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Traits for the single-phase flow Navier-Stokes mass model
 */
struct NavierStokesMassOnePModelTraits
{
    //! There are as many momentum balance equations as dimensions
    //! and one mass balance equation.
    static constexpr int numEq() { return 1; }

    //! The number of phases is 1
    static constexpr int numFluidPhases() { return 1; }

    //! The number of components is 1
    static constexpr int numFluidComponents() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! The one-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return false; }

    //! The model is isothermal
    static constexpr bool enableEnergyBalance() { return false; }

    //! The model does not include a turbulence model
    static constexpr bool usesTurbulenceModel() { return false; }

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::none; }

    //! the indices
    using Indices = NavierStokesMassOnePIndices;
};

/*!
 * \ingroup NavierStokesModel
 * \brief Traits class for the volume variables of the Navier-Stokes model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam MT The model traits
 */
template<class PV,
         class FSY,
         class FST,
         class MT>
struct NavierStokesMassOnePVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using ModelTraits = MT;
};

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal Navier-Stokes model
struct NavierStokesMassOneP{ using InheritsFrom = std::tuple<ModelProperties>; };
struct NavierStokesMassOnePNI{ using InheritsFrom = std::tuple<NavierStokesMassOneP>; };
} // end namespace TTag


template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOneP>
{ using type = NavierStokesMassOnePModelTraits; };

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::NavierStokesMassOneP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = Dumux::ImmiscibleFluidState<Scalar, FluidSystem>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassOneP>
{ using type = NavierStokesMassOnePLocalResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOneP>
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
    using type = NavierStokesMassOnePVolumeVariables<Traits>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::NavierStokesMassOneP>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    struct DiffusiveFluxTypes {}; // no diffusion
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
public:
    using type = NavierStokesMassOnePFluxVariables<
        Problem, ModelTraits, DiffusiveFluxTypes, ElementVolumeVariables, ElementFluxVariablesCache
    >;
};

// ! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOneP> { using type = NavierStokesIOFields; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::NavierStokesMassOneP>
{
public:
    using type = struct EmptyCouplingManager {};
};

// Set the default spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::NavierStokesMassOneP>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FreeFlowDefaultSpatialParams<GridGeometry, Scalar>;
};
///////////////////////////////////////////////////////////////////////////
// Properties for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! Add temperature to the output
template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOnePNI>
{ using type = NavierStokesEnergyIOFields<NavierStokesIOFields>; };

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOnePNI>
{ using type = NavierStokesEnergyModelTraits<NavierStokesMassOnePModelTraits>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOnePNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using BaseTraits = NavierStokesMassOnePVolumeVariablesTraits<PV, FSY, FST, MT>;
    using ETCM = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using HCT = GetPropType<TypeTag, Properties::HeatConductionType>;
    struct NITraits : public BaseTraits
    {
        using EffectiveThermalConductivityModel = ETCM;
        using HeatConductionType = HCT;
    };
public:
    using type = NavierStokesMassOnePVolumeVariables<NITraits>;
};

//! Use the average for effective conductivities
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::NavierStokesMassOnePNI>
{
    struct type
    {
        template<class VolVars>
        static auto effectiveThermalConductivity(const VolVars& volVars)
        {
            return volVars.fluidThermalConductivity();
        }
    };
};

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::NavierStokesMassOnePNI>
{ using type = FouriersLaw<TypeTag>; };

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::NavierStokesMassOnePNI>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    struct DiffusiveFluxTypes
    {
        using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;
    };

    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;

public:
    using type = NavierStokesMassOnePFluxVariables<
        Problem, ModelTraits, DiffusiveFluxTypes, ElementVolumeVariables, ElementFluxVariablesCache
    >;
};

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::NavierStokesMassOnePNI>
{
    struct type : public GetPropType<TypeTag, Properties::HeatConductionType>::Cache
    {};
};

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::NavierStokesMassOnePNI>
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    static constexpr bool diffusionIsSolDependent = false; // no diffusion;
    static constexpr bool heatConductionIsSolDependent
        = getPropValue<TypeTag, Properties::SolutionDependentHeatConduction>();

    using type = FreeFlowScalarFluxVariablesCacheFiller<
        Problem, ModelTraits, diffusionIsSolDependent, heatConductionIsSolDependent
    >;
};

template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::NavierStokesMassOnePNI>
{ static constexpr bool value = true; };

} // end namespace Properties
} // end namespace Dumux

#endif
