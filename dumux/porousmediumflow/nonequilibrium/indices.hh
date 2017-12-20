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
 * \brief The primary variable and equation indices for the MpNc model.
 */
#ifndef DUMUX_NONEQUILIBRIUM_INDICES_HH
#define DUMUX_NONEQUILIBRIUM_INDICES_HH

#include "model.hh"

namespace Dumux
{

/*!
 * \ingroup NonEquilibrium
 * \ingroup ImplicitIndices
 * \brief The primary variable and equation indices for the MpNc model.
 */
template <class TypeTag, int BasePVOffset = 0>
class NonEquilbriumIndices: public GET_PROP_TYPE(TypeTag, EquilibriumIndices)
{
public:
     using FluidSystem  = typename GET_PROP_TYPE(TypeTag, FluidSystem);
     enum { numPhases = FluidSystem::numPhases };
     enum { numEnergyEqFluid = GET_PROP_VALUE(TypeTag, NumEnergyEqFluid) };
     enum { numEnergyEqSolid = GET_PROP_VALUE(TypeTag, NumEnergyEqSolid) };
     static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);

    /*!
     * \brief Index for the temperature of the wetting phase in a vector of primary
     *        variables.
     */
    static const unsigned int temperature0Idx = BasePVOffset + numEq - numEnergyEqFluid - numEnergyEqSolid;

    /*!
     * \brief Index for the temperature of the solid phase in a vector of primary
     *        variables.
     */
    static const unsigned int temperatureSolidIdx = BasePVOffset + numEq - numEnergyEqSolid;
    /*!
     * \brief Compatibility with non kinetic models
     */
    static const unsigned int temperatureIdx = temperature0Idx;
    /*!
     * \brief Equation index of the energy equation.
     */
    static const unsigned int energyEq0Idx = BasePVOffset + numEq - numEnergyEqFluid - numEnergyEqSolid;
    /*!
     * \brief Compatibility with non kinetic models
     */
    static const unsigned int energyEqIdx = energyEq0Idx;

    /*!
     * \brief Equation index of the energy equation.
     */
    static const unsigned int energyEqSolidIdx = BasePVOffset + numEq - numEnergyEqSolid;
};

}

#endif
