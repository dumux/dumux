// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePModel
 * \brief A single-phase, isothermal flow model using the fully implicit scheme.
 *
 * Single-phase, isothermal flow model, which uses a standard Darcy approach as the
 * equation for the conservation of momentum. For details on Darcy's law see dumux/flux/darcyslaw.hh.
 *
 * Furthermore, it solves the mass continuity equation
 * \f[
 \frac{\partial (\phi \varrho) }{\partial t} + \nabla \cdot \left\lbrace
 - \varrho \frac{\textbf K}{\mu} \left( \nabla  p -\varrho {\textbf g} \right) \right\rbrace = q,
 * \f]
* where:
 * * \f$ \phi \f$ is the porosity of the porous medium,
 * * \f$ \varrho \f$ is the mass density,
 * * \f$ \textbf{K} \f$ is the intrinsic permeability tensor,
 * * \f$ \mu \f$ represents the dynamic viscosity,
 * * \f$  p \f$ is the pressure,
 * * \f$ \textbf{g} \f$ is the gravitational acceleration vector,
 * * and \f$ q \f$ is a source or sink term.
 *
 *
 * The model supports compressible as well as incompressible fluids.
 */

#ifndef DUMUX_1P_MODEL_HH
#define DUMUX_1P_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup OnePModel
 * \brief Specifies a number properties of single-phase models.
 */
struct OnePModelTraits
{
    //! Per default, we use the indices without offset
    using Indices = OnePIndices<>;

    static constexpr int numEq() { return 1; }
    static constexpr int numFluidPhases() { return 1; }
    static constexpr int numFluidComponents() { return 1; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
    static constexpr bool enableThermalDispersion() { return false; }
};

/*!
 * \ingroup OnePModel
 * \brief Traits class for the volume variables of the single-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct OnePVolumeVariablesTraits
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
// Create new type tags
namespace TTag {
//! The type tags for the isothermal single phase model
struct OneP { using InheritsFrom = std::tuple<PorousMediumFlow>; };

//! The type tags for the non-isothermal single phase model
struct OnePNI { using InheritsFrom = std::tuple<OneP>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// Properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
template<class TypeTag>
struct IOFields<TypeTag, TTag::OneP> { using type = OnePIOFields; };                          //!< default I/O fields specific to this model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OneP> { using type = ImmiscibleLocalResidual<TypeTag>; }; //!< the local residual function
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::OneP> { using type = OnePModelTraits; };                //!< states some specifics of the one-phase model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OneP> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; }; //!< default the actually used traits to the base traits

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::OneP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    using Traits = OnePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = OnePVolumeVariables<Traits>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state.
 *
 * This should be chosen appropriately for the model ((non-)isothermal,
 * equilibrium, ...). This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::OneP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

///////////////////////////////////////////////////////////////////////////
// Properties for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! Add temperature to the output
template<class TypeTag>
struct IOFields<TypeTag, TTag::OnePNI> { using type = EnergyIOFields<OnePIOFields>; };

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OnePNI> { using type = PorousMediumFlowNIModelTraits<OnePModelTraits>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::OnePNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using BaseTraits = OnePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;

    using ETCM = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };

public:
    using type = OnePVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

//! Use the average for effective conductivities
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::OnePNI>
{ using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Properties
} // end namespace Dumux

#endif
