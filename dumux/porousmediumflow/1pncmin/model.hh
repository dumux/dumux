// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
* \file
* \ingroup OnePNCMinModel
* \brief A single-phase, multi-component model considering mineralization processes.
*
* This model implements one-phase n-component flow of a compressible fluid composed of
* the n components \f$\kappa \f$ in combination with mineral precipitation and dissolution
* of the solid phases. The standard multi-phase Darcy
* approach is used as the equation for the conservation of momentum.
* For details on Darcy's law see dumux/flux/darcyslaw.hh.
*
* By inserting Darcy's law into the equations for the conservation of the
* components, one gets one transport equation for each component,
* \f[
 \frac{\partial ( \phi \varrho_f X^\kappa )}
{\partial t} -  \nabla \cdot \left\{ \varrho_f X^\kappa
\frac{k_{r}}{\mu} \mathbf{K}
(\nabla  p - \varrho_{f}  \mathbf{g}) \right\}
- \nabla \cdot \left\{{\bf D_{pm}^\kappa} \varrho_{f} \nabla  X^\kappa \right\}
-  q_\kappa = 0 \qquad \kappa \in \{w, a,\cdots \},
* \f]
* where:
* * \f$ \phi \f$ is the  porosity,
* * \f$ \varrho_f \f$ is the mass density of the fluid,
* * \f$ X^\kappa \f$ is the mass fraction of component \f$ \kappa \f$ in the fluid,
* * \f$ k_{r} \f$ is the relative permeability,
* * \f$ \mu \f$ represents the dynamic viscosity,
* * \f$ \mathbf{K} \f$ is the intrinsic permeability tensor,
* * \f$ p \f$ is the pressure,
* * \f$ \mathbf{g} \f$ is the gravitational acceleration vector,
* * \f$ {\bf D_{pm}^\kappa} \f$ is the effective diffusivity in the porous medium,
* * \f$ q_\kappa \f$ is a source or sink term.
*
* The solid or mineral phases are assumed to consist of a single component.
* Their mass balance consists of only a storage and a source term,
* \f[
 \frac{\partial ( \varrho_\lambda \phi_\lambda )} {\partial t} = q_\lambda,
* \f]
* where:
* * \f$ \varrho_\lambda\f$ represents the mass density of the solid phase,
* * \f$ \phi_\lambda \f$ is the porosity of the solid,
* * \f$ q_\lambda \f$ is a source or sink term representing.
*
* The primary variables are the pressure \f$p\f$ and the mole fractions of the
* dissolved components \f$x^k\f$. The primary variable of the solid phases is the volume
* fraction
\f$\phi_\lambda = \frac{V_\lambda}{V_{total}}\f$,
*
* where:
* * \f$ V_\lambda \f$ is the volume of phase \f$ \lambda \f$,
* * \f$ V_{\text{total}} \f$ is the total volume of the system.
*
* The source an sink terms link the mass balances of the n-transported component to the
* solid phases. The porosity \f$\phi\f$ is updated according to the reduction of the initial
* (or solid-phase-free porous medium) porosity \f$\phi_0\f$ by the accumulated volume
* fractions of the solid phases,
* \f$ \phi = \phi_0 - \sum (\phi_\lambda),\f$
*
* where:
* * \f$ \phi \f$ represents the remaining porosity in the system,
* * \f$ \phi_0 \f$ is the initial porosity,
* * \f$ \phi_\lambda \f$ denotes the volume fraction of phase \f$ \lambda \f$.
* Additionally, the permeability is updated depending on the current porosity.
*/

#ifndef DUMUX_1PNCMIN_MODEL_HH
#define DUMUX_1PNCMIN_MODEL_HH

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/1pnc/indices.hh>
#include <dumux/porousmediumflow/1pnc/volumevariables.hh>

#include <dumux/material/solidstates/compositionalsolidstate.hh>

#include <dumux/porousmediumflow/mineralization/model.hh>
#include <dumux/porousmediumflow/mineralization/localresidual.hh>
#include <dumux/porousmediumflow/mineralization/volumevariables.hh>
#include <dumux/porousmediumflow/mineralization/iofields.hh>

#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
// Create new type tags
namespace TTag {
struct OnePNCMin { using InheritsFrom = std::tuple<OnePNC>; };
struct OnePNCMinNI { using InheritsFrom = std::tuple<OnePNCMin>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags for the isothermal 1pncmin model
//////////////////////////////////////////////////////////////////

//! use the mineralization volume variables together with the 1pnc vol vars
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::OnePNCMin>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    using BaseTraits = OnePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;

    using CDM = GetPropType<TypeTag, Properties::CompositionalDispersionModel>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    template<class BaseTraits, class CDM, class DT, class EDM>
    struct NCTraits : public BaseTraits
    {
        using CompositionalDispersionModel = CDM;
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
    };

    using NonMinVolVars = OnePNCVolumeVariables<NCTraits<BaseTraits, CDM, DT, EDM>>;
public:
    using type = MineralizationVolumeVariables<NCTraits<BaseTraits, CDM, DT, EDM>, NonMinVolVars>;
};

// Use the mineralization local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePNCMin> { using type = MineralizationLocalResidual<TypeTag>; };

//! Use non-mineralization model traits with 1pnc traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OnePNCMin>
{
private:
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using NonMinTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
public:
    using type = MineralizationModelTraits<NonMinTraits, SolidSystem::numComponents, SolidSystem::numInertComponents>;
};

//! The two-phase model uses the immiscible fluid state
template<class TypeTag>
struct SolidState<TypeTag, TTag::OnePNCMin>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = CompositionalSolidState<Scalar, SolidSystem>;
};

//! Use the mineralization vtk output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::OnePNCMin> { using type = MineralizationIOFields<OnePNCIOFields>; };

//////////////////////////////////////////////////////////////////
// Properties for the non-isothermal 1pncmin model
//////////////////////////////////////////////////////////////////

//! Non-isothermal vtk output
template<class TypeTag>
struct IOFields<TypeTag, TTag::OnePNCMinNI>
{
    using MineralizationIOF = MineralizationIOFields<OnePNCIOFields>;
    using type = EnergyIOFields<MineralizationIOF>;
};

//! The non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OnePNCMinNI>
{
private:
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using OnePNCTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
    using IsothermalTraits = MineralizationModelTraits<OnePNCTraits, SolidSystem::numComponents, SolidSystem::numInertComponents>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! use the mineralization volume variables together with the 1pnc vol vars
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::OnePNCMinNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    using BaseTraits = OnePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;

    using CDM = GetPropType<TypeTag, Properties::CompositionalDispersionModel>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;
    template<class BaseTraits, class CDM, class DT, class EDM, class ETCM>
    struct NCNITraits : public BaseTraits
    {
        using CompositionalDispersionModel = CDM;
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
        using EffectiveThermalConductivityModel = ETCM;
    };
    using NonMinVolVars = OnePNCVolumeVariables<NCNITraits<BaseTraits, CDM, DT, EDM, ETCM>>;
public:
    using type = MineralizationVolumeVariables<NCNITraits<BaseTraits, CDM, DT, EDM, ETCM>, NonMinVolVars>;
};
//! Use the average for effective conductivities
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::OnePNCMinNI>
{ using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Properties
} // end namespace Dumux

#endif
