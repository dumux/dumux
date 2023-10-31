// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ExtendedRichardsModel
 * \brief This model implements a variant of the extended Richards'
 *        equation for quasi-twophase flow (see e.g. Vanderborght et al. 2017).
 *
 * The extended Richards' equation
 \f[
 \frac{\partial (\phi S_w \varrho_w) }{\partial t}
 +
 \frac{\partial (\phi (1-S_w)\varrho_n X_n^w ) }{\partial t}
 -
 \nabla \cdot \left\lbrace
 \varrho_w \frac{k_{rw}}{\mu_w} \; \mathbf{K} \;
 \left( \text{\nabla}
 p_w - \varrho_w \textbf{g}
 \right)
 +
 {\bf D_{n, pm}^w} \varrho_n \nabla  X^w_n
 \right\rbrace
 =
 q_w,
 \f]
 * where:
 * * \f$ \phi \f$ is the porosity of the porous medium,
 * * \f$ S_w \f$ represents the saturation of the wetting-phase,
 * * \f$ \varrho_w \f$ is the mass density of the wetting phase,
 * * \f$ \varrho_n \f$ is the mass density of the non-wetting phase,
 * * \f$ k_{rw} \f$ is the relative permeability of the wetting phase,
 * * \f$ \mu_w \f$ is the dynamic viscosity of the wetting phase,
 * * \f$ \mathbf{K} \f$ is the intrinsic permeability tensor,
 * * \f$ p_w \f$ is the pressure of the wetting phase,
 * * \f$ \mathbf{g} \f$ is the gravitational acceleration vector,
 * * \f$ \bf D_{n,pm}^{w} \f$ is the effective diffusivity of water in the non-wetting phase,
 * * \f$ X_n^w \f$ is the mass fraction of water in the non-wetting phase,
 * * \f$ q_w \f$ is a source or sink term in the wetting phase,
 *
 * additionally models water vapor diffusion in the gas phase.
 * The model is derived based on the two-phase flow equations
 * based on the assumption that the gas phase does not move but
 * and remains at constant pressure.
 */

#ifndef DUMUX_RICHARDSEXTENDED_MODEL_HH
#define DUMUX_RICHARDSEXTENDED_MODEL_HH

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
#include <dumux/porousmediumflow/richards/balanceequationopts.hh>
#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/porousmediumflow/richards/velocityoutput.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"
#include "localresidual.hh"

namespace Dumux {

/*!
 * \ingroup ExtendedRichardsModel
 * \brief Specifies a number properties of the extended Richards model.
 */
struct ExtendedRichardsModelTraits : public RichardsModelTraits
{
    using Indices = ExtendedRichardsIndices;

    static constexpr bool enableMolecularDiffusion() { return true; }
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
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT, class DT, class EDM>
struct ExtendedRichardsVolumeVariablesTraits
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
struct ExtendedRichards { using InheritsFrom = std::tuple<Richards>; };
struct ExtendedRichardsNI { using InheritsFrom = std::tuple<ExtendedRichards>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////

//! The local residual operator
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ExtendedRichards> { using type = ExtendedRichardsLocalResidual<TypeTag>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::ExtendedRichards>
{
    using type = ExtendedRichardsIOFields;
};

//! The model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ExtendedRichards> { using type = ExtendedRichardsModelTraits; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ExtendedRichards>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using Traits = ExtendedRichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;
public:
    using type = ExtendedRichardsVolumeVariables<Traits>;
};

//! Use the model after Millington (1961) for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::ExtendedRichards>
{ using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//! The primary variables vector for the richards model
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ExtendedRichards>
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                     GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

/////////////////////////////////////////////////////
// Property values for non-isothermal Richars model
/////////////////////////////////////////////////////

//! set non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ExtendedRichardsNI>
{
private:
    using IsothermalTraits = ExtendedRichardsModelTraits;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the vtk output fields specific to th non-isothermal model
template<class TypeTag>
struct IOFields<TypeTag, TTag::ExtendedRichardsNI>
{
    using type = EnergyIOFields<ExtendedRichardsIOFields>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ExtendedRichardsNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using BaseTraits = ExtendedRichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;

    using ETCM = GetPropType< TypeTag, Properties::ThermalConductivityModel>;
    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };

public:
    using type = ExtendedRichardsVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::ExtendedRichardsNI>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
