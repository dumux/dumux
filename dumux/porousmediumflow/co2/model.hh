// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CO2Model
 * \brief Adaption of the non-isothermal two-phase two-component flow model to problems with CO2
 *
 *  See @ref TwoPTwoCModel for reference to the equations used.
 *  The CO2 model is derived from the 2p2c model, however, in the CO2 model the phase switch criterion
 *  is different from the 2p2c model.
 *  The phase switch occurs when the equilibrium concentration
 *  of a component in a phase is exceeded, instead of the sum of the components in the virtual phase
 *  (the phase which is not present) being greater that unity as done in the 2p2c model.
 *  The CO2VolumeVariables do not use a constraint solver for calculating the mole fractions as is the
 *  case in the 2p2c model. Instead, mole fractions are calculated in the FluidSystem with a given
 *  temperature, pressure and salinity.
 *  The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false,
 *  but make sure that the according units are used in the problem setup. useMoles is set to false by default.
 */

#ifndef DUMUX_TWOP_TWOC_CO2_MODEL_HH
#define DUMUX_TWOP_TWOC_CO2_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>
#include "volumevariables.hh"

namespace Dumux {
namespace Properties {

// Create new type tags
namespace TTag {
struct TwoPTwoCCO2 { using InheritsFrom = std::tuple<TwoPTwoC>; };
struct TwoPTwoCCO2NI { using InheritsFrom = std::tuple<TwoPTwoCNI>; };
} // end namespace TTag

//! the co2 volume variables use the same traits as the 2p2c model
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPTwoCCO2>
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
    using type = TwoPTwoCCO2VolumeVariables<NCTraits<BaseTraits, DT, EDM>>;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPTwoCCO2NI>
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
    using type = TwoPTwoCCO2VolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
