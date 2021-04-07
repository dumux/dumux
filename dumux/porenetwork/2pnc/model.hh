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

#ifndef DUMUX_PNM2PNC_MODEL_HH
#define DUMUX_PNM2PNC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porenetwork/properties.hh>

#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/material/spatialparams/porenetwork/porenetwork2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include <dumux/flux/porenetwork/advection.hh>
#include <dumux/flux/porenetwork/fickslaw.hh>
#include <dumux/porenetwork/2p/fluxvariablescache.hh>
#include <dumux/porenetwork/2p/gridfluxvariablescache.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>

#include "iofields.hh"
#include "volumevariables.hh"

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
    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = Dumux::TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = Dumux::TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

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

template<class TypeTag>
struct Formulation<TypeTag, TTag::PNMTwoPNC>
{ static constexpr auto value = TwoPFormulation::p0s1; };

//! We use fick's law as the default for the diffusive fluxes
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::PNMTwoPNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
public:
    using type = Dumux::PoreNetwork::PNMFicksLaw<Scalar, ModelTraits::numFluidPhases(), ModelTraits::numFluidComponents()>;
};


///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! The spatial parameters to be employed.
//! Use PNMTwoPSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMTwoPNC>
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
struct AdvectionType<TypeTag, TTag::PNMTwoPNC>
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
struct IOFields<TypeTag, TTag::PNMTwoPNC> { using type = Dumux::PoreNetwork::TwoPNCIOFields; };

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::PNMTwoPNC> { using type = Dumux::EnergyLocalResidual<TypeTag> ; };

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PNMTwoPNC>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using Traits = Dumux::PoreNetwork::PNMTwoPDefaultGridFVCTraits<Problem, FluxVariablesCache, Indices, Labels>;
public:
    using type = Dumux::PoreNetwork::PNMTwoPGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache, Traits>;
};

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMTwoPNC>
{ using type = Dumux::PoreNetwork::TwoPFluxVariablesCache<GetPropType<TypeTag, Properties::AdvectionType>>; };

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
    using IsothermalTraits = Dumux::TwoPNCModelTraits<FluidSystem::numComponents,
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
    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = Dumux::TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = Dumux::TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

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
    using type = Dumux::PoreNetwork::TwoPNCVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
};

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMTwoPNCNI>
{
    using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>;
}; //!< Use the average for effective conductivities

//! Set non-isothermal output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMTwoPNCNI> { using type = Dumux::EnergyIOFields<Dumux::PoreNetwork::TwoPNCIOFields>; };


} // end namespace
}


#endif
