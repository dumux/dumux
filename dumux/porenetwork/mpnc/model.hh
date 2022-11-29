// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Base class for all models which use the one-phase,
 *        fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase flow model.
 */

#ifndef DUMUX_PNM_MPNC_MODEL_HH
#define DUMUX_PNM_MPNC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porenetwork/properties.hh>
#include <dumux/flux/porenetwork/advection.hh>
#include <dumux/flux/porenetwork/fickslaw.hh>
#include <dumux/flux/porenetwork/fourierslaw.hh>

#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include <dumux/porousmediumflow/mpnc/model.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/material/spatialparams/porenetwork/porenetwork2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include <dumux/porousmediumflow/compositional/localresidual.hh>

#include "iofields.hh"
#include "volumevariables.hh"

#include "fluxvariablescache.hh"
#include "gridfluxvariablescache.hh"

namespace Dumux
{
// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit single-phase problems
// Create new type tags
namespace TTag {
struct PNMMPNC { using InheritsFrom = std::tuple<PoreNetworkModel, MPNC>; };

//! The type tags for the corresponding non-isothermal problems
struct PNMMPNCNI { using InheritsFrom = std::tuple<PNMMPNC>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMMPNC>
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

    using Traits = Dumux::MPNCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;
    static_assert(MT::numFluidPhases() == 2, "Only two-phase flow is supported");
public:
    using type = Dumux::PoreNetwork::MPNCVolumeVariables<Traits>;
};

//! We use fick's law as the default for the diffusive fluxes
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::PNMMPNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
public:
    using type = Dumux::PoreNetwork::PNMFicksLaw<Scalar, ModelTraits::numFluidPhases(), ModelTraits::numFluidComponents()>;
};

//////////////////////////////////////////////////////////////////
// Properties for the isothermal 2pnc model
//////////////////////////////////////////////////////////////////
//! The primary variables vector for the 2pnc model
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::PNMMPNC>
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                     GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    using type = Dumux::SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! The spatial parameters to be employed.
//! Use PNMTwoPSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMMPNC>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalRules = Dumux::PoreNetwork::FluidMatrix::MultiShapeTwoPLocalRules<Scalar>;
public:
    using type = Dumux::PoreNetwork::TwoPDefaultSpatialParams<GridGeometry, Scalar, LocalRules>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMMPNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using S = Dumux::PoreNetwork::TransmissibilityPatzekSilin<Scalar>;
    using W = Dumux::PoreNetwork::WettingLayerTransmissibility::RansohoffRadke<Scalar>;
    using N = Dumux::PoreNetwork::NonWettingPhaseTransmissibility::BakkeOren<Scalar>;

public:
    using type = Dumux::PoreNetwork::CreepingFlow<Scalar, S, W, N>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMMPNC> { using type = Dumux::PoreNetwork::MPNCIOFields; };

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::PNMMPNC> { using type = Dumux::EnergyLocalResidual<TypeTag> ; };

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PNMMPNC>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using Traits = Dumux::PoreNetwork::PNMMPNCDefaultGridFVCTraits<Problem, FluxVariablesCache>;
public:
    using type = Dumux::PoreNetwork::PNMMPNCGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache, Traits>;
};

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMMPNC>
{ using type = Dumux::PoreNetwork::MPNCFluxVariablesCache<GetPropType<TypeTag, Properties::AdvectionType>>; };

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PNMMPNCNI>
{
private:
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IsothermalTraits = MPNCModelTraits<FluidSystem::numPhases,
                                             FluidSystem::numComponents,
                                             getPropValue<TypeTag, Properties::PressureFormulation>(),
                                             getPropValue<TypeTag, Properties::UseMoles>(),
                                             getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
}; //!< The model traits of the non-isothermal model

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMMPNCNI>
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
    using BaseTraits = Dumux::MPNCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;

    using ETCM = GetPropType< TypeTag, Properties::ThermalConductivityModel>;
    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };
    using Traits = NITraits<BaseTraits, ETCM>;
public:
    using type = Dumux::PoreNetwork::MPNCVolumeVariables<Traits>;
};

//! Set non-isothermal output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMMPNCNI> { using type = EnergyIOFields<Dumux::PoreNetwork::MPNCIOFields>; };

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMMPNCNI>
{
    using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>;
}; //!< Use the average for effective conductivities

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::PNMMPNCNI>
{
private:
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
public:
    using type = Dumux::PoreNetwork::PNMFouriersLaw<typename FluxVariables::MolecularDiffusionType>;
};

} // end namespace
}


#endif // DUMUX_PNM_MPNC_MODEL_HH
