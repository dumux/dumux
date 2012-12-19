// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by                                                   *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Adaption of the BOX scheme to the non-isothermal two-phase two-component flow model without constraint solver.
 */
#ifndef DUMUX_CO2NI_MODEL_HH
#define DUMUX_CO2NI_MODEL_HH

#include <dumux/implicit/co2/co2model.hh>

namespace Dumux {
/*!
 * \ingroup CO2NIModel
 * \brief Adaption of the BOX scheme to the non-isothermal two-phase two-component flow model.
 *   See TwoPTwoCNI model for reference to the equations.
 *   The CO2NI model is derived from the CO2 model. In the CO2 model the phase switch criterion
 *   is different from the 2p2c model.
 *   The phase switch occurs when the equilibrium concentration
 *   of a component in a phase is exceeded instead of the sum of the components in the virtual phase
 *   (the phase which is not present) being greater that unity as done in the 2p2c model.
 *   The CO2VolumeVariables do not use a constraint solver for calculating the mole fractions as is the
 *   case in the 2p2c model. Instead mole fractions are calculated in the FluidSystem with a given
 *   temperature, pressure and salinity.
 *
 */
template<class TypeTag>
class CO2NIModel : public CO2Model<TypeTag>
{
};

}

#include <dumux/implicit/2p2cni/2p2cnipropertydefaults.hh>

#endif
