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
 * \brief The indices for the thermal non-equilibrium part of the MpNc model.
 */
#ifndef DUMUX_MPNC_INDICES_ENERGY_KINETIC_HH
#define DUMUX_MPNC_INDICES_ENERGY_KINETIC_HH

#include <dumux/porousmediumflow/mpnc/implicit/indices.hh>

namespace Dumux
{

/*!
 * \brief The indices required for the energy equation. Specialization for the case of
 *        *3* energy balance equations.
 */
template <int PVOffset>
struct MPNCEnergyIndices<PVOffset, /*enableEnergy=*/true, /*numEnergyEquations=*/3>
{
public:
    /*!
     * \brief This module defines one new primary variable.
     */
    static const unsigned int numPrimaryVars = 3;

    /*!
     * \brief Index for the temperature of the wetting phase in a vector of primary
     *        variables.
     */
    static const unsigned int temperature0Idx = PVOffset + 0;

    /*!
     * \brief Compatibility with non kinetic models
     */
    static const unsigned int temperatureIdx = temperature0Idx;
    /*!
     * \brief Equation index of the energy equation.
     */
    static const unsigned int energyEq0Idx = PVOffset + 0;
    /*!
     * \brief Compatibility with non kinetic models
     */
    static const unsigned int energyEqIdx = energyEq0Idx;
};

/*!
 * \brief The indices required for the energy equation. Specialization for the case of
 *        *2* energy balance equations.
 */
template <int PVOffset>
struct MPNCEnergyIndices<PVOffset, /*enableEnergy=*/true, /*numEnergyEquations=*/2>
{
public:
    /*!
     * \brief This module defines one new primary variable.
     */
    static const unsigned int numPrimaryVars = 2;

    /*!
     * \brief Index for the temperature of the wetting phase in a vector of primary
     *        variables.
     */
    static const unsigned int temperature0Idx = PVOffset + 0;

    /*!
     * \brief Compatibility with non kinetic models
     */
    static const unsigned int temperatureIdx = temperature0Idx;
    /*!
     * \brief Equation index of the energy equation.
     */
    static const unsigned int energyEq0Idx = PVOffset + 0;
    /*!
     * \brief Equation index of the energy equation.
     */
    static const unsigned int energyEqSolidIdx = energyEq0Idx + numPrimaryVars - 1  ;

    /*!
     * \brief Index for storing e.g. temperature fluidState
     */
    static const unsigned int  temperatureFluidIdx = 0 ;
    static const unsigned int  temperatureSolidIdx = 1 ;


    /*!
     * \brief Compatibility with non kinetic models
     */
    static const unsigned int energyEqIdx = energyEq0Idx;


};

/*!
 * \brief The indices required for the energy equation.
 */
template <int PVOffset, bool isNonIsothermal>
struct MPNCEnergyIndices<PVOffset, isNonIsothermal, /*numEnergyEquations=*/1 >
{
public:
    /*!
     * \brief This module defines one new primary variable.
     */
    static const unsigned int numPrimaryVars = 1;

    /*!
     * \brief Index for the temperature in a vector of primary
     *        variables.
     */
    static const unsigned int temperatureIdx = PVOffset + 0;
    /*!
     * \brief Equation index of the energy equation.
     */
    static const unsigned int energyEqIdx = PVOffset + 0;
};

/*!
 * \brief The indices for the energy equation.
 *
 * This is a dummy class for the isothermal case.
 */
template <int PVOffset>
struct MPNCEnergyIndices<PVOffset, /*isNonIsothermal=*/false, /*numEnergyEquations*/0>
{
public:
    /*!
     * \brief This module does not define any primary variables in the
     *        isothermal case.
     */
    static const unsigned int numPrimaryVars = 0;

    /*!
     * \brief Equation index of the temperature primary variable. This
     *        is a dummy value which hopefully makes the simulation
     *        crash if used.
     */
    static const unsigned int temperatureIdx = -1;

    /*!
     * \brief Equation index of the energy equation. This is a dummy
     *        value which hopefully makes the simulation crash if used.
     */
    static const unsigned int energyEqIdx = -1;
};



}

#endif
