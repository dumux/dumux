// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsModel
 * \brief This model implements a variant of the Richards'
 *        equation for quasi-twophase flow.
 *
 * In the unsaturated zone, Richards' equation
 * is frequently used to approximate the water distribution
 * above the groundwater level (in the unsaturated zone):
 \f[
 \frac{\partial (\phi S_w \varrho_w)}{\partial t}
 -
 \nabla \cdot \left\lbrace
 \varrho_w \frac{k_{rw}}{\mu_w} \; \mathbf{K} \;
 \left( \nabla
 p_w - \varrho_w \textbf{g}
 \right)
 \right\rbrace
 =
 q_w,
 \f]
 *
 * where:
 * * \f$ \phi \f$ is the porosity of the porous medium,
 * * \f$ S_w \f$ represents the water saturation,
 * * \f$ \varrho_w \f$ is the water density,
 * * \f$ k_{rw} \f$ is the relative permeability of the water phase,
 * * \f$ \mu_w \f$ is the dynamic viscosity of the water phase,
 * * \f$ \mathbf{K} \f$ is the intrinsic permeability tensor,
 * * \f$ p_w \f$ is the liquid water pressure,
 * * \f$ \mathbf{g} \f$ is the gravitational acceleration vector,
 * * \f$ q_w \f$ is a source or sink term.
 *
 * It can be derived from the two-phase flow equations.
 * In contrast to the full two-phase model, the Richards model assumes
 * gas as the nonwetting fluid and that it exhibits a much lower
 * viscosity than the (liquid) wetting phase. (For example at
 * atmospheric pressure and at room temperature, the viscosity of air
 * is only about \f$1\%\f$ of the viscosity of liquid water.) As a
 * consequence, the mobility (\f$\frac{k_{r}}{\mu}\f$) is
 * typically much larger for the gas phase than for the wetting
 * phase. For this reason, the Richards model assumes that
 * gas phase mobility is infinitely large. This implies that
 * the pressure of the gas phase is equivalent to the static pressure
 * distribution and that therefore, mass conservation only needs to be
 * considered for the wetting phase.
 *
 * The model thus chooses the absolute pressure of the wetting phase
 * \f$p_w\f$ as its only primary variable. The wetting phase
 * saturation is calculated using the inverse of the capillary
 * pressure, i.e.
 \f[
 S_w = p_c^{-1}(p_g - p_w)
 \f]
 * holds, where \f$p_g\f$ is a given reference gas pressure. Nota bene,
 * that the last step is assumes that the capillary
 * pressure-saturation curve can be uniquely inverted, so it is not
 * possible to set the capillary pressure to zero when using the
 * Richards model!
 */

#ifndef DUMUX_RICHARDS_MODEL_HH
#define DUMUX_RICHARDS_MODEL_HH

#include <type_traits>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"
#include "localresidual.hh"
#include "velocityoutput.hh"
#include "balanceequationopts.hh"

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Specifies a number properties of the Richards model.
 *
 * \tparam enableDiff specifies if diffusion of water in air is to be considered.
 */
struct RichardsModelTraits
{
    using Indices = RichardsIndices;

    static constexpr int numEq() { return 1; }
    static constexpr int numFluidPhases() { return 2; }
    static constexpr int numFluidComponents() { return 1; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }

    //! The Richards model has some assumptions on the fluid systems
    //! that can be verified with this trait
    template<class FluidSystem>
    static constexpr bool fluidSystemIsCompatible()
    {
        return !FluidSystem::isGas(FluidSystem::phase0Idx)
            && FluidSystem::isGas(FluidSystem::phase1Idx);
    }

    //! The Richards model has some assumptions on the fluid systems
    //! that can be verified with this trait
    template<class FluidSystem>
    static constexpr auto checkFluidSystem(const FluidSystem& fs)
    {
        struct FluidSystemCheck {
            static_assert(fluidSystemIsCompatible<FluidSystem>(),
                "Richards model currently assumes the first phase to be liquid and the second phase to be gaseous.");
        };
        return FluidSystemCheck{};
    }
};

/*!
 * \ingroup RichardsModel
 * \brief Traits class for the Richards model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct RichardsVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
};

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal Richards model.
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
// Create new type tags
namespace TTag {
struct Richards { using InheritsFrom = std::tuple<PorousMediumFlow>; };
struct RichardsNI { using InheritsFrom = std::tuple<Richards>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////

//! The local residual operator
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Richards> { using type = RichardsLocalResidual<TypeTag>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::Richards> { using type = RichardsIOFields; };

template<class TypeTag>
struct VelocityOutput<TypeTag, TTag::Richards>
{
    using type = RichardsVelocityOutput<
        GetPropType<TypeTag, Properties::GridVariables>,
        GetPropType<TypeTag, Properties::FluxVariables>
    >;
};

//! The model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Richards> { using type = RichardsModelTraits; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Richards>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = RichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = RichardsVolumeVariables<Traits>;
};

//! Use the model after Millington (1961) for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::Richards>
{ using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the H2O-Air fluid system with Simple H2O (constant density and viscosity).
 */
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Richards>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::H2OAir<Scalar,
                                      Components::SimpleH2O<Scalar>,
                                      FluidSystems::H2OAirDefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::Richards>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

//! Set a richards specific class for the balance equation options
template<class TypeTag>
struct BalanceEqOpts<TypeTag, TTag::Richards>
{ using type = RichardsBalanceEquationOptions<GetPropType<TypeTag, Properties::FluidSystem>>; };

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::RichardsNI>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

/////////////////////////////////////////////////////
// Property values for non-isothermal Richars model
/////////////////////////////////////////////////////

//! set non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::RichardsNI>
{
private:
    using IsothermalTraits = RichardsModelTraits;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the vtk output fields specific to th non-isothermal model
template<class TypeTag>
struct IOFields<TypeTag, TTag::RichardsNI> { using type = EnergyIOFields<RichardsIOFields>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::RichardsNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using BaseTraits = RichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;

    using ETCM = GetPropType< TypeTag, Properties::ThermalConductivityModel>;
    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };

public:
    using type = RichardsVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
