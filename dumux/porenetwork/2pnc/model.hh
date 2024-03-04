// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPNCModel
 * \brief A two-phase multi-compositional pore-network model using fully implicit scheme.
 *
 * A mass balance equation is formulated for each pore body \f$i\f$ and each component \f$\kappa\f$:
 *
 * \f[
 *	V_i \frac{\partial (\sum_{\alpha} (x^\kappa_\alpha \varrho_\alpha S_\alpha)_i) }{\partial t} + \sum_\alpha \sum_j (x^\kappa_\alpha \varrho_\alpha  Q_\alpha)_{ij}
 *  + \sum_\alpha \sum_j (A_\alpha j^\kappa_\alpha)_{ij}= (V q^\kappa)_i ~.
 * \f]
 *
 * \f$V_i\f$ is the pore body volume, \f$x^\kappa_\alpha\f$ resepresents the mole fraction of component \f$\kappa\f$  in phase \f$\alpha\f$ , the advective mass flow \f$(\varrho_\alpha Q_\alpha)_{ij}\f$
 * through throat \f$ij\f$ can be based on the fluid phase density \f$\varrho\f$ either of the upstream pore body \f$i\f$ or \f$j\f$ (upwinding) or on the respective averaged value.
 * \f$j^\kappa_\alpha\f$ represents the diffuive flux vector and \f$A_\alpha\f$ the area of throat cross section occupied by phase \f$\alpha\f$. \f$q_\alpha\f$ is a mass sink or source term defined on pore body \f$i\f$.
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
 *	\Psi_\alpha = \varrho_\alpha \mathbf{g} \cdot (\mathbf{x}_i - \mathbf{x}_j),
 * \f]
 *
 * where \f$\mathbf{x}_i - \mathbf{x}_j\f$  is the distance vector between the centers of pore bodies \f$i\f$ and \f$j\f$ and \f$\mathbf{g}\f$ is the gravity vector.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables are the pressure \f$p_\alpha\f$ for phase \f$\alpha\f$ and the mole or mass fraction of dissolved components in phase \f$\alpha\f$ \f$x^\kappa_\alpha\f$ or \f$X^\kappa_\alpha\f$.
 */

#ifndef DUMUX_PNM_2P_NC_MODEL_HH
#define DUMUX_PNM_2P_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porenetwork/properties.hh>
#include <dumux/flux/porenetwork/fickslaw.hh>

#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include <dumux/flux/porenetwork/advection.hh>
#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include <dumux/porenetwork/2p/fluxvariablescache.hh>
#include <dumux/porenetwork/2p/gridfluxvariablescache.hh>
#include <dumux/porenetwork/2p/spatialparams.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>


#include "iofields.hh"
#include "volumevariables.hh"

namespace Dumux::Properties{
// \{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit two-phase n-components problems
// Create new type tags
namespace TTag {
struct PNMTwoPNC { using InheritsFrom = std::tuple<PoreNetworkModel, TwoPNC>; };

//! The type tags for the corresponding non-isothermal problems
struct PNMTwoPNCNI { using InheritsFrom = std::tuple<PNMTwoPNC>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMTwoPNC>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;

    template<class BaseTraits, class DT, class EDM>
    struct NCTraits : public BaseTraits
    {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
    };
public:
    using type = PoreNetwork::TwoPNCVolumeVariables<NCTraits<BaseTraits, DT, EDM>>;
};

//! The primary variables vector for the 2pnc model
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::PNMTwoPNC>
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                     GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! The spatial parameters to be employed.
//! Use PNMTwoPSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMTwoPNC>
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
struct AdvectionType<TypeTag, TTag::PNMTwoPNC>
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
struct IOFields<TypeTag, TTag::PNMTwoPNC> { using type = PoreNetwork::TwoPNCIOFields; };

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PNMTwoPNC>
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

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMTwoPNC>
{ using type = PoreNetwork::TwoPFluxVariablesCache<GetPropType<TypeTag, Properties::AdvectionType>>; };

//! We use fick's law as the default for the diffusive fluxes
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::PNMTwoPNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
public:
    using type = PoreNetwork::PNMFicksLaw<Scalar, ModelTraits::numFluidPhases(), ModelTraits::numFluidComponents()>;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PNMTwoPNCNI>
{
private:
    //! we use the number of components specified by the fluid system here
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static_assert(FluidSystem::numPhases == 2, "Only fluid systems with 2 fluid phases are supported by the 2p-nc model!");
    using IsothermalTraits = TwoPNCModelTraits<FluidSystem::numComponents,
                                               getPropValue<TypeTag, Properties::UseMoles>(),
                                               getPropValue<TypeTag, Properties::SetMoleFractionsForFirstPhase>(),
                                               getPropValue<TypeTag, Properties::Formulation>(), getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
}; //!< The model traits of the non-isothermal model

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMTwoPNCNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

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
    using type = PoreNetwork::TwoPNCVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
};

//! Set non-isothermal output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMTwoPNCNI> { using type = EnergyIOFields<PoreNetwork::TwoPNCIOFields>; };

// Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMTwoPNCNI>
{
    using type = ThermalConductivitySomertonTwoP<GetPropType<TypeTag, Properties::Scalar>>;
};

} // end namespace Dumux::Properties

#endif
