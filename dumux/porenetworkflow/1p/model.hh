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
 * \brief A one-phase-flow, isothermal pore-network model using the fully implicit scheme.
 */

#ifndef DUMUX_PNM1P_MODEL_HH
#define DUMUX_PNM1P_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/flux/porenetwork/advection.hh>
#include <dumux/porenetworkflow/properties.hh>

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/spatialparams/porenetwork/porenetwork1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>

#include "iofields.hh"
#include "volumevariables.hh"
#include "fluxvariablescache.hh"

namespace Dumux {
/*!
 * \ingroup PoreNetworkOnePModel
 * \brief A one-phase, isothermal pore-network model using the fully implicit scheme.
 *        TODO docme
 */

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
struct PNMOneP { using InheritsFrom = std::tuple<PoreNetworkModel, OneP>; };

//! The type tags for the corresponding non-isothermal problems
struct PNMOnePNI { using InheritsFrom = std::tuple<PNMOneP>; };
} // end namespace TTag

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMOneP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    using Traits = Dumux::OnePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = Dumux::PoreNetwork::OnePVolumeVariables<Traits>;
};

//! The spatial parameters to be employed.
//! Use PNMOnePSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMOneP>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::PoreNetwork::OnePDefaultSpatialParams<GridGeometry, Scalar>;
};

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PoreNetworkModel>
{ using type = Dumux::PoreNetwork::OnePFluxVariablesCache<GetPropType<TypeTag, Properties::AdvectionType>>; };

template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMOneP> { using type = Dumux::PoreNetwork::OnePIOFields; };

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOneP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = Dumux::PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type = Dumux::PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMOnePNI>
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
    using type = Dumux::PoreNetwork::OnePVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMOnePNI> { using type = EnergyIOFields<Dumux::PoreNetwork::OnePIOFields>; }; //!< Add temperature to the output

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PNMOnePNI> { using type = PorousMediumFlowNIModelTraits<OnePModelTraits>; }; //!< The model traits of the non-isothermal model

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMOnePNI> { using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>; }; //!< Use the average for effective conductivities

} // end namespace Properies

} //namespace Dumux


#endif
