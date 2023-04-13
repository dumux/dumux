// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMOnePNCModel
 * \brief Adaption of the fully implicit scheme to the one-phase n-component pore network model.
 *
 * A mass balance equation is formulated for each pore body \f$i\f$ and each component \f$\kappa\f$:
 *
 * \f[
 *	V_i \frac{\partial (\varrho_{i} X^\kappa)}{\partial t} + \sum_j (\varrho X^\kappa Q)_{ij} = (V q^\kappa)_i ~.
 * \f]
 *
 * \f$V_i\f$ is the pore body volume, and the advective mass flow \f$(\varrho Q X^\kappa)_{ij}\f$ through throat \f$ij\f$ can be based on the fluid phase density
 * \f$\varrho\f$ either of the upstream pore body \f$i\f$ or \f$j\f$ (upwinding) or on the respective averaged value. \f$q_i\f$ is a mass sink or source
 * term defined on pore body \f$i\f$.
 *
 * Per default, the volume flow rate \f$Q_{ij}\f$ follows a linear Hagen-Poiseuille-type law (PoreNetworkModel::CreepingFlow) which is only valid for \f$Re < 1\f$:
 *
 * \f[
 *	Q_{ij} = g_{ij} (p_{i} - p_{j} + \Psi)  ~.
 * \f]
 *
 * \f$g_{ij}\f$ is a suitable throat conductance value (see e.g. PoreNetwork::TransmissibilityPatzekSilin) while \f$p_i\f$ and \f$p_j\f$ are averaged pore body pressures.
 *
 * The (optional) influence of gravity is given by
 *
 * \f[
 *	\Psi = \varrho \mathbf{g} (\mathbf{x_i} - \mathbf{x_j}) ~,
 * \f]
 *
 * where \f$\mathbf{x_i} - \mathbf{x_j}\f$  is the distance vector between the centers of pore bodies \f$i\f$ and \f$j\f$ and \f$\mathbf{g}\f$ is the gravitational acceleration.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables are the pressure \f$p\f$ and the mass or mole fraction of dissolved components \f$X^\kappa\f$ or \f$x^\kappa\f$.
 */

#ifndef DUMUX_PNM1PNC_MODEL_HH
#define DUMUX_PNM1PNC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/flux/porenetwork/fickslaw.hh>

#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/porenetwork/properties.hh>
#include <dumux/porenetwork/1p/model.hh>
#include <dumux/porenetwork/1p/spatialparams.hh>

#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/1pnc/iofields.hh>

#include "iofields.hh"
#include "volumevariables.hh"

namespace Dumux::Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit single-phase mulit-component problems
// Create new type tags
namespace TTag {
struct PNMOnePNC { using InheritsFrom = std::tuple<PoreNetworkModel, OnePNC>; };

//! The type tags for the corresponding non-isothermal problems
struct PNMOnePNCNI { using InheritsFrom = std::tuple<PNMOnePNC>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! The spatial parameters to be employed.
//! Use PNMOnePSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMOnePNC>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::OnePDefaultSpatialParams<GridGeometry, Scalar>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type = PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

//! Set as default that no component mass balance is replaced by the total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::PNMOnePNC>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    static constexpr auto value = FluidSystem::numComponents;
};

//! We use fick's law as the default for the diffusive fluxes
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::PNMOnePNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
public:
    using type = PoreNetwork::PNMFicksLaw<Scalar, ModelTraits::numFluidPhases(), ModelTraits::numFluidComponents()>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMOnePNC>
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

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    template<class BaseTraits, class DT, class EDM>
    struct NCTraits : public BaseTraits
    {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
    };

public:
    using type = PoreNetwork::OnePNCVolumeVariables<NCTraits<BaseTraits, DT, EDM>>;
};

//!< Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMOnePNC>
{ using type = PoreNetwork::OnePNCIOFields<GetPropType<TypeTag, Properties::FluidSystem>>; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::PNMOnePNC> { static constexpr bool value = true; };

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//! model traits of the non-isothermal model.
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PNMOnePNCNI>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using CompositionalDispersionModel = GetPropType<TypeTag, Properties::CompositionalDispersionModel>;
    using ThermalDispersionModel = GetPropType<TypeTag, Properties::CompositionalDispersionModel>;
    using IsothermalTraits = OnePNCModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>(),
                                                                           getPropValue<TypeTag, Properties::EnableCompositionalDispersion>(),
                                                                           getPropValue<TypeTag, Properties::EnableThermalDispersion>(),
                                                                           getPropValue<TypeTag, Properties::ReplaceCompEqIdx>(),
                                                                           CompositionalDispersionModel>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits, ThermalDispersionModel>;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMOnePNCNI>
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

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;
    template<class BaseTraits, class DT, class EDM, class ETCM>
    struct NCNITraits : public BaseTraits
    {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
        using EffectiveThermalConductivityModel = ETCM;
    };

public:
    using type = PoreNetwork::OnePNCVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMOnePNCNI>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IsothermalFields = PoreNetwork::OnePNCIOFields<FluidSystem>;
public:
    using type = EnergyIOFields<IsothermalFields>;
};

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMOnePNCNI>
{
    using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>;
}; //!< Use the average for effective conductivities

// template<class TypeTag>
// struct HeatConductionType<TypeTag, TTag::PNMOnePNCNI>
// {
//     TODO uncomment this as soon as there is a generalized approach for component enthalpies in all fluid systems
//     using type = PoreNetwork::PNMFouriersLaw<GetPropType<TypeTag, MolecularDiffusionType>>;
// }; //!< Use Fourier's law and also consider enthalpy transport by molecular diffusion

} // end namespace Dumux::Properties
#endif
