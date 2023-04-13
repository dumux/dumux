// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \brief A two-phase-flow, isothermal pore-network model using the fully implicit scheme.
 *
 * A mass balance equation is formulated for each pore body \f$i\f$ and each phase \f$\alpha\f$:
 *
 * \f[
 *	V_i \frac{\partial (\varrho_\alpha S_\alpha)_i}{\partial t} + \sum_j (\varrho_\alpha  Q_\alpha)_{ij} = (V q_\alpha)_i ~.
 * \f]
 *
 * \f$V_i\f$ is the pore body volume, and the advective mass flow \f$(\varrho_\alpha Q_\alpha)_{ij}\f$ through throat \f$ij\f$ can be based on the fluid phase density
 * \f$\varrho\f$ either of the upstream pore body \f$i\f$ or \f$j\f$ (upwinding) or on the respective averaged value. \f$q_\alpha\f$ is a mass sink or source
 * term defined on pore body \f$i\f$.
 *
 * Per default, the volume flow rate \f$Q_{\alpha,ij}\f$ follows a linear Hagen-Poiseuille-type law (PoreNetworkModel::CreepingFlow) which is only valid for \f$Re < 1\f$:
 *
 * \f[
 *	Q_{\alpha,ij} = g_{\alpha, ij} (p_{\alpha, i} - p_{\alpha, j} + \Psi_\alpha)  ~.
 * \f]
 *
 * \f$g_{\alpha,ij}\f$ is a suitable throat conductance value that takes into account the presence/saturation of the individual phases while \f$p_{\alpha,i}\f$ and \f$p_{\alpha,j}\f$ are averaged pore body phase pressures.
 *
 * The (optional) influence of gravity is given by
 *
 * \f[
 *	\Psi_\alpha = \varrho_\alpha \mathbf{g} (\mathbf{x_i} - \mathbf{x_j}) ~,
 * \f]
 *
 * where \f$\mathbf{x_i} - \mathbf{x_j}\f$  is the distance vector between the centers of pore bodies \f$i\f$ and \f$j\f$ and \f$\mathbf{g}\f$ is the gravitational acceleration.
 *
 * The primary variables are the wetting phase pressure and the the nonwetting phase saturation (\f$p_w\f$ and \f$S_n\f$) or the nonwetting phase pressure and the the wetting phase saturation (\f$p_n\f$ and \f$S_w\f$),
 * depending on the chose formulation (see TwoPModel).
 */

#ifndef DUMUX_PNM2P_MODEL_HH
#define DUMUX_PNM2P_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/flux/porenetwork/advection.hh>
#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/porenetwork/properties.hh>

#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/2p/model.hh>

#include "fluxvariablescache.hh"
#include "gridfluxvariablescache.hh"
#include "iofields.hh"
#include "volumevariables.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {
// \{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit two-phase problems
// Create new type tags
namespace TTag {
struct PNMTwoP { using InheritsFrom = std::tuple<PoreNetworkModel, TwoP>; };

//! The type tags for the corresponding non-isothermal problems
struct PNMTwoPNI { using InheritsFrom = std::tuple<PNMTwoP>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal two-phase model
///////////////////////////////////////////////////////////////////////////

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMTwoP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of non-wetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;

    using Traits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;
public:
    using type = PoreNetwork::TwoPVolumeVariables<Traits>;
};

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMTwoP>
{ using type = PoreNetwork::TwoPFluxVariablesCache<GetPropType<TypeTag, Properties::AdvectionType>>; };

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PNMTwoP>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
    using Traits = PoreNetwork::PNMTwoPDefaultGridFVCTraits<Problem, FluxVariablesCache>;
public:
    using type = PoreNetwork::PNMTwoPGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache, Traits>;
};

//! The spatial parameters to be employed.
//! Use PoreNetwork::TwoPDefaultSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMTwoP>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalRules = PoreNetwork::FluidMatrix::MultiShapeTwoPLocalRules<Scalar>;
public:
    using type = PoreNetwork::TwoPDefaultSpatialParams<GridGeometry, Scalar, LocalRules>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMTwoP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using S = PoreNetwork::TransmissibilityPatzekSilin<Scalar>;
    using W = PoreNetwork::WettingLayerTransmissibility::RansohoffRadke<Scalar>;
    using N = PoreNetwork::NonWettingPhaseTransmissibility::BakkeOren<Scalar>;

public:
    using type = PoreNetwork::CreepingFlow<Scalar, S, W, N>;
};

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::PNMTwoP> { using type = Dumux::EnergyLocalResidual<TypeTag> ; };

template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMTwoP> { using type = PoreNetwork::TwoPIOFields; };

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//! model traits of the non-isothermal model.
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PNMTwoPNI>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IsothermalTraits = TwoPModelTraits<getPropValue<TypeTag, Properties::Formulation>()>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMTwoPNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of non-wetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;

    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };

public:
    using type = PoreNetwork::TwoPVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

//! Set the vtk output fields specific to the non-isothermal two-phase model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMTwoPNI> { using type = EnergyIOFields<PoreNetwork::TwoPIOFields>; };

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMTwoPNI>
{
    using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>;
}; //!< Use the average for effective conductivities

} // end namespace Dumux::Properties

#endif
