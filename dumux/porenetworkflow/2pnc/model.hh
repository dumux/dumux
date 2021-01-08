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
#include <dumux/flux/porenetwork/advection.hh>
#include <dumux/porenetworkflow/properties.hh>

#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>

#include <dumux/porenetworkflow/2p/model.hh>
#include <dumux/material/spatialparams/porenetwork/porenetwork2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porenetworkflow/2p/volumevariables.hh>
#include <dumux/porenetworkflow/2p/iofields.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>

#include "iofields.hh"
#include "volumevariables.hh"
#include "gridfluxvariablescache.hh"
#include "indices.hh"

namespace Dumux
{

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal two phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

// Create new type tags
namespace TTag {
struct PNMTwoPNC { using InheritsFrom = std::tuple<PoreNetworkModel, TwoPNC>; };

//! The type tags for the corresponding non-isothermal problems
struct PNMTwoPNCNI { using InheritsFrom = std::tuple<PNMTwoPNC>; };
} // end namespace TTag


//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::PNMTwoPNC>
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                     GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

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
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of non-wetting phase saturations
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
    using type = PNMTwoPNCVolumeVariables<NCTraits<BaseTraits, DT, EDM>>;
};


//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMTwoPNC> { using type = PNMTwoPNCIOFields; };

//! Default formulation is pw-Sn, overwrite if necessary
template<class TypeTag>
struct Formulation<TypeTag, TTag::PNMTwoPNC>
{ static constexpr auto value = TwoPFormulation::p0s1; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::PNMTwoPNC> { static constexpr bool value = true; };                         //!< Use mole fractions in the balance equations by default

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMTwoPNC>
{ using type = PNMTwoPFluxVariablesCache<GetPropType<TypeTag, Properties::AdvectionType>>; };

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
    using Traits = PNMTwoPDefaultGridFVCTraits<Problem, FluxVariablesCache, Indices, Labels>;
public:
    using type = PNMTwoPGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache, Traits>;
};

//! The spatial parameters to be employed.
//! Use PNMTwoPSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMTwoPNC>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalRules = RegularizedTwoPLocalRulesCubeJoekarNiasar<Scalar>;
public:
    using type = PNMTwoPDefaultSpatialParams<GridGeometry, Scalar, LocalRules>;

};

template<class TypeTag>
struct SinglePhaseTransmissibilityLaw<TypeTag, TTag::PNMTwoPNC>
{
    using type = TransmissibilityPatzekSilin<GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct WettingLayerTransmissibilityLaw<TypeTag, TTag::PNMTwoPNC>
{
    using type = WettingLayerTransmissibility::RansohoffRadke<GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct NonWettingPhaseTransmissibilityLaw<TypeTag, TTag::PNMTwoPNC>
{
    using type = NonWettingPhaseTransmissibility::BakkeOren<GetPropType<TypeTag, Properties::Scalar>>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMTwoPNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using S = GetPropType<TypeTag, Properties::SinglePhaseTransmissibilityLaw>;
    using W = GetPropType<TypeTag, Properties::WettingLayerTransmissibilityLaw>;
    using N = GetPropType<TypeTag, Properties::NonWettingPhaseTransmissibilityLaw>;

public:
    using type = Dumux::PoreNetworkCreepingFlow<Scalar, S, W, N>;
};


//! The labels
template<class TypeTag>
struct Labels<TypeTag, TTag::PNMTwoPNC> { using type = Dumux::Labels; };

//! Set whether a pc of 0 should be applied to pores with sw = 1
template<class TypeTag>
struct ZeroPc<TypeTag, TTag::PNMTwoPNC> { static constexpr bool value = true; };

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::PNMTwoPNC> { using type = Dumux::EnergyLocalResidual<TypeTag> ; };

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//! model traits of the non-isothermal model.
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PNMTwoPNCNI>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static_assert(FluidSystem::numPhases == 2, "Only fluid systems with 2 fluid phases are supported by the 2p-nc model!");
    using IsothermalTraits = TwoPNCModelTraits<FluidSystem::numComponents,
                                   getPropValue<TypeTag, Properties::UseMoles>(),
                                   getPropValue<TypeTag, Properties::SetMoleFractionsForFirstPhase>(),
                                   getPropValue<TypeTag, Properties::Formulation>(), getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

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
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of non-wetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;

    template<class BaseTraits,class DT, class EDM,  class ETCM>
    struct NCNITraits : public BaseTraits {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
        using EffectiveThermalConductivityModel = ETCM;
    };

public:
    using type = PNMTwoPNCVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
    //using type = PNMTwoPNCVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

//! Set the vtk output fields specific to the non-isothermal twop model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMTwoPNCNI> { using type = EnergyIOFields<PNMTwoPNCIOFields>; };

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMTwoPNCNI>
{
    using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>;
}; //!< Use the average for effective conductivities


} // end namespace
}


#endif
