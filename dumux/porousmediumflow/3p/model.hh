// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePModel
 * \brief Adaption of the fully implicit scheme to the three-phase flow model.
 *
 * This model implements three-phase flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$.
 * The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum.
 * For details on Darcy's law see dumux/flux/darcyslaw.hh.
 *
 * By inserting Darcy's law into the equations for the conservation
 * of the phase mass, one gets
 \f[
 \frac{\partial (\phi \varrho_\alpha S_\alpha )}{\partial t}
 -
 \nabla \cdot \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\nabla  p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 \right\} - q_\alpha = 0,
 \f]
 * where:
 * * \f$ \phi \f$ is the porosity of the porous medium,
 * * \f$ S_\alpha \f$ represents the saturation of phase \f$ \alpha \f$,
 * * \f$ \varrho_\alpha \f$ is the mass density of phase \f$ \alpha \f$,
 * * \f$ k_{r\alpha} \f$ is the relative permeability of phase \f$ \alpha \f$,
 * * \f$ \mu_\alpha \f$ is the dynamic viscosity of phase \f$ \alpha \f$,
 * * \f$ \mathbf{K} \f$ is the intrinsic permeability tensor,
 * * \f$ p_\alpha \f$ is the pressure of phase \f$ \alpha \f$,
 * * \f$ \mathbf{g} \f$ is the gravitational acceleration vector,
 * * \f$ q_\alpha \f$ is a source or sink term.
 *
 * The model uses commonly applied auxiliary conditions like
 * \f$S_w + S_n + S_g = 1\f$ for the saturations.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are gas phase pressure \f$p_g\f$,
 * water saturation \f$S_w\f$ and NAPL saturation \f$S_n\f$.
 */
#ifndef DUMUX_3P_MODEL_HH
#define DUMUX_3P_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup ThreePModel
 * \brief Specifies a number properties of three-phase models.
 */
struct ThreePModelTraits
{
    using Indices = ThreePIndices;

    static constexpr int numEq() { return 3; }
    static constexpr int numFluidPhases() { return 3; }
    static constexpr int numFluidComponents() { return 3; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
};

/*!
 * \ingroup ThreePModel
 * \brief Traits class for the two-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct ThreePVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
};

namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
namespace TTag {
//! The type tags for the isothermal three-phase model
struct ThreeP { using InheritsFrom = std::tuple<PorousMediumFlow>; };
//! The type tags for the non-isothermal three-phase model
struct ThreePNI { using InheritsFrom = std::tuple<ThreeP>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Properties for the isothermal 3p model
//////////////////////////////////////////////////////////////////

//! Set the model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ThreeP>
{
 private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p model!");
    static_assert(FluidSystem::numComponents == 3, "Only fluid systems with 3 components are supported by the 3p model!");
 public:
    using type = ThreePModelTraits;
};

//! The local residual function of the conservation equations
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ThreeP> { using type = ImmiscibleLocalResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ThreeP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    using Traits = ThreePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = ThreePVolumeVariables<Traits>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state.
 *
 *  The fluid state should be chosen appropriately for the model ((non-)isothermal, equilibrium, ...).
 *  This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::ThreeP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::ThreeP> { using type = ThreePIOFields; };

/////////////////////////////////////////////////
// Properties for the non-isothermal 3p model
/////////////////////////////////////////////////

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::ThreePNI>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

//! Set non-isothermal output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::ThreePNI> { using type = EnergyIOFields<ThreePIOFields>; };

//! Set non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ThreePNI>
{
private:
   using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
   static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p model!");
   static_assert(FluidSystem::numComponents == 3, "Only fluid systems with 3 components are supported by the 3p model!");
public:
   using type = PorousMediumFlowNIModelTraits<ThreePModelTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ThreePNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using BaseTraits = ThreePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;

    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };
public:
    using type = ThreePVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
