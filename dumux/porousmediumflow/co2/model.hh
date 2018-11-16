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
 * \brief Adaption of the fully implicit scheme to the CO2Model model.
 */
#ifndef DUMUX_TWOP_TWOC_CO2_MODEL_HH
#define DUMUX_TWOP_TWOC_CO2_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>
#include "primaryvariableswitch.hh"
#include "volumevariables.hh"

/*!
 * \file
 * \ingroup CO2Model
 * \brief Adaption of the non-isothermal two-phase two-component flow model to problems with CO2
 *
 *   TODO: Put a doxgyen link refernce here
 *   See TwoPTwoCModel for reference to the equations used.
 *   The CO2 model is derived from the 2p2c model. In the CO2 model the phase switch criterion
 *   is different from the 2p2c model.
 *   The phase switch occurs when the equilibrium concentration
 *   of a component in a phase is exceeded, instead of the sum of the components in the virtual phase
 *   (the phase which is not present) being greater that unity as done in the 2p2c model.
 *   The CO2VolumeVariables do not use a constraint solver for calculating the mole fractions as is the
 *   case in the 2p2c model. Instead mole fractions are calculated in the FluidSystem with a given
 *   temperature, pressure and salinity.
 *   The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 *   problem file. Make sure that the according units are used in the problem setup. useMoles is set to false by default.
 *
 */
namespace Dumux {

/*!
 * \ingroup CO2Model
 * \brief Traits class for the two-phase two-component CO2 model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct TwoPTwoCCO2VolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
};

namespace Properties {

// Create new type tags
namespace TTag {
struct TwoPTwoCCO2 { using InheritsFrom = std::tuple<TwoPTwoC>; };
struct TwoPTwoCCO2NI { using InheritsFrom = std::tuple<TwoPTwoCNI>; };
} // end namespace TTag

//! the CO2 privarswitch and VolumeVariables properties
template<class TypeTag>
struct PrimaryVariableSwitch<TypeTag, TTag::TwoPTwoCCO2> { using type = TwoPTwoCCO2PrimaryVariableSwitch; };
template<class TypeTag>
struct PrimaryVariableSwitch<TypeTag, TTag::TwoPTwoCCO2NI> { using type = TwoPTwoCCO2PrimaryVariableSwitch; };

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
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    using Traits = TwoPTwoCCO2VolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = TwoPTwoCCO2VolumeVariables< Traits >;
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
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    using Traits = TwoPTwoCCO2VolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = TwoPTwoCCO2VolumeVariables< Traits >;
};

} // end namespace Properties
} // end namespace Dumux

#endif
