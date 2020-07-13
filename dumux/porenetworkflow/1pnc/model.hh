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
 * \brief  Adaption of the fully implicit scheme to the one-phase n-component pore network model.
 */

#ifndef DUMUX_PNM1PNC_MODEL_HH
#define DUMUX_PNM1PNC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/flux/porenetwork/fickslaw.hh>
#include <dumux/porenetworkflow/properties.hh>

#include <dumux/porenetworkflow/1p/model.hh>
#include <dumux/material/spatialparams/porenetwork/porenetwork1p.hh>

#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/1pnc/iofields.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

#include <dumux/material/spatialparams/porenetwork/porenetwork1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>

#include <dumux/material/fluidstates/immiscible.hh>

#include "iofields.hh"
#include "volumevariables.hh"

namespace Dumux
{

namespace Properties
{
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
    using type = Dumux::PoreNetwork::OnePDefaultSpatialParams<GridGeometry, Scalar>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = Dumux::PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type = Dumux::PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
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
    using type = Dumux::PoreNetwork::PNMFicksLaw<Scalar, ModelTraits::numFluidPhases(), ModelTraits::numFluidComponents()>;
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
    using BaseTraits = Dumux::OnePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    template<class BaseTraits, class DT, class EDM>
    struct NCTraits : public BaseTraits
    {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
    };

public:
    using type = Dumux::PoreNetwork::OnePNCVolumeVariables<NCTraits<BaseTraits, DT, EDM>>;
};

//!< Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMOnePNC>
{ using type = Dumux::PoreNetwork::OnePNCIOFields<GetPropType<TypeTag, Properties::FluidSystem>>; };

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
    using IsothermalTraits = OnePNCModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>(), getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
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
    using BaseTraits = Dumux::OnePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;

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
    using type = Dumux::PoreNetwork::OnePNCVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMOnePNCNI>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IsothermalFields = Dumux::PoreNetwork::OnePNCIOFields<FluidSystem>;
public:
    using type = EnergyIOFields<IsothermalFields>;
};

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::PNMOnePNCNI>
{
    using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>;
}; //!< Use the average for effective conductivities

}

}

#endif
