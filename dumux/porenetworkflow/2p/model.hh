// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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

#ifndef DUMUX_PNM2P_MODEL_HH
#define DUMUX_PNM2P_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/flux/porenetwork/advection.hh>
#include <dumux/porenetworkflow/properties.hh>

#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/material/spatialparams/porenetwork/porenetwork2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/regularizedlocalrules.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/regularizedlocalrulesforcube.hh>

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include "fluxvariablescache.hh"
#include "gridfluxvariablescache.hh"
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
struct PNMTwoP { using InheritsFrom = std::tuple<PoreNetworkModel, TwoP>; };

//! The type tags for the corresponding non-isothermal problems
struct PNMTwoPNI { using InheritsFrom = std::tuple<PNMTwoP>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
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

    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of non-wetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;

    using Traits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;
public:
    using type = PNMTwoPVolumeVariables<Traits>;
};

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMTwoP>
{ using type = PNMTwoPFluxVariablesCache<GetPropType<TypeTag, Properties::AdvectionType>>; };

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::PNMTwoP>
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
struct SpatialParams<TypeTag, TTag::PNMTwoP>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    // using LocalRules = RegularizedTwoPLocalRules<Scalar, getPropValue<TypeTag, Properties::ZeroPc>()>;
    using LocalRules = RegularizedTwoPLocalRulesCubeJoekarNiasar<Scalar>;
public:
    using type = PNMTwoPDefaultSpatialParams<GridGeometry, Scalar, LocalRules>;
};

template<class TypeTag>
struct SinglePhaseTransmissibilityLaw<TypeTag, TTag::PNMTwoP>
{
    using type = TransmissibilityPatzekSilin<GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct WettingLayerTransmissibilityLaw<TypeTag, TTag::PNMTwoP>
{
    using type = WettingLayerTransmissibility::RansohoffRadke<GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct NonWettingPhaseTransmissibilityLaw<TypeTag, TTag::PNMTwoP>
{
    using type = NonWettingPhaseTransmissibility::BakkeOren<GetPropType<TypeTag, Properties::Scalar>>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMTwoP>
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
struct Labels<TypeTag, TTag::PNMTwoP> { using type = Dumux::Labels; };

//! Set whether a pc of 0 should be applied to pores with sw = 1
template<class TypeTag>
struct ZeroPc<TypeTag, TTag::PNMTwoP> { static constexpr bool value = true; };

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::PNMTwoP> { using type = Dumux::EnergyLocalResidual<TypeTag> ; };

template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMTwoP> { using type = PNMTwoPIOFields; };

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
    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of non-wetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;

    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };

public:
    using type = PNMTwoPVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

//! Set the vtk output fields specific to the non-isothermal twop model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMTwoPNI> { using type = EnergyIOFields<PNMTwoPIOFields>; };

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMTwoPNI>
{
    using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>;
}; //!< Use the average for effective conductivities


} // end namespace
}


#endif
