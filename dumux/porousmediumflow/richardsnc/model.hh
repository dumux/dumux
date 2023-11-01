// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsNCModel
 * \brief Base class for all models which use the Richards,
 *        n-component fully implicit model.
 *
 * This extension of Richards' equation, allows for
 * the wetting phase to consist of multiple components:
 *\f{eqnarray*}
 && \frac{\partial (\sum_w \varrho_w X_w^\kappa \phi S_w )}
 {\partial t}
 - \sum_w  \nabla \cdot \left\{ \varrho_w X_w^\kappa
 \frac{k_{rw}}{\mu_w} \mathbf{K}
 (\nabla  p_w - \varrho_{w}  \mathbf{g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_w \nabla \cdot \left\{{\bf D_{w, pm}^\kappa} \varrho_{w} \nabla  X^\kappa_{w} \right\}
 - \sum_w q_w^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 w \in \{w, g\},
 \f}
 * where:
 * * \f$ \phi \f$ is the porosity of the porous medium,
 * * \f$ S_w \f$ represents the saturation of the wetting phase,
 * * \f$ \varrho_w \f$ is the mass density of the wetting phase,
 * * \f$ k_{rw} \f$ is the relative permeability of the wetting phase,
 * * \f$ \mu_w \f$ is the dynamic viscosity of the wetting phase,
 * * \f$ \mathbf{K} \f$ is the intrinsic permeability tensor,
 * * \f$ p_w \f$ is the pressure of the wetting phase,
 * * \f$ \mathbf{g} \f$ is the gravitational acceleration vector,
 * * \f$ \bf D_{w,pm}^{k} \f$ is the effective diffusivity of component \f$ \kappa \f$ in the wetting phase,
 * * \f$ X_w^k \f$ is the mass fraction of component \f$ \kappa \f$ in the wetting phase,
 * * \f$ q_w \f$ is a source or sink term in the wetting phase.
 */

#ifndef DUMUX_RICHARDSNC_MODEL_HH
#define DUMUX_RICHARDSNC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/richards/model.hh>

#include "volumevariables.hh"
#include "indices.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup RichardsNCModel
 * \brief Specifies a number properties of the Richards n-components model.
 *
 * \tparam nComp the number of components to be considered.
 * \tparam useMol whether to use mass or mole balances
 */
template<int nComp, bool useMol, int repCompEqIdx = nComp>
struct RichardsNCModelTraits
{
    using Indices = RichardsNCIndices;

    static constexpr int numEq() { return nComp; }
    static constexpr int numFluidPhases() { return 1; }
    static constexpr int numFluidComponents() { return nComp; }
    static constexpr int replaceCompEqIdx() { return repCompEqIdx; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }
    static constexpr bool enableCompositionalDispersion() { return false; }
    static constexpr bool enableThermalDispersion() { return false; }

    static constexpr bool useMoles() { return useMol; }
};

/*!
 * \ingroup RichardsNCModel
 * \brief Traits class for the Richards n-components model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 * \tparam DT The diffusion type
 * \tparam EDM The effective diffusivity model
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT, class DT, class EDM>
struct RichardsNCVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
    using DiffusionType = DT;
    using EffectiveDiffusivityModel = EDM;
};

namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
// Create new type tags
namespace TTag {
struct RichardsNC { using InheritsFrom = std::tuple<PorousMediumFlow>; };
struct RichardsNCNI { using InheritsFrom = std::tuple<RichardsNC>; };
} // end namespace TTag
//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the model traits class
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::RichardsNC>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = RichardsNCModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>(), getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
};
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::RichardsNC> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; };

//! Define that per default mole fractions are used in the balance equations
template<class TypeTag>
struct UseMoles<TypeTag, TTag::RichardsNC> { static constexpr bool value = true; };

//! Use the dedicated local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::RichardsNC> { using type = CompositionalLocalResidual<TypeTag>; };

//! We set the replaceCompIdx to 0, i.e. the first equation is substituted with
//! the total mass balance, i.e. the phase balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::RichardsNC> { static constexpr int value = 0; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::RichardsNC>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using NCTraits = RichardsNCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;

public:
    using type = RichardsNCVolumeVariables<NCTraits>;
};

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the liquid phase fluid system with simple H2O.
 */
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RichardsNC>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::SimpleH2O<Scalar>, Components::Constant<1, Scalar>>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::RichardsNC>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::RichardsNC> { using type = RichardsNCIOFields; };

//! The model after Millington (1961) is used for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::RichardsNC> { using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//! average is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::RichardsNCNI> { using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>; };

//////////////////////////////////////////////////////////////////
// Property values for non-isothermal Richards n-components model
//////////////////////////////////////////////////////////////////

//! set non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::RichardsNCNI>
{
private:
    using IsothermalTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::RichardsNCNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using BaseTraits = RichardsNCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;

    using ETCM = GetPropType< TypeTag, Properties::ThermalConductivityModel>;
    template<class BaseTraits, class ETCM>
    struct NCNITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };
public:
    using type = RichardsNCVolumeVariables<NCNITraits<BaseTraits, ETCM>>;
};


} // end namespace Properties
} // end namespace Dumux

#endif
