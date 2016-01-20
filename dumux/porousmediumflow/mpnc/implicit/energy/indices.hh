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
 * \brief The indices for the non-isothermal part of the compositional
 *        multi-phase model.
 */
#ifndef DUMUX_MPNC_INDICES_ENERGY_HH
#define DUMUX_MPNC_INDICES_ENERGY_HH

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitIndices
 * \brief The indices for the energy equation.
 *
 * This is a dummy class for the isothermal case.
 */
template <int PVOffset, bool enableEnergy/*=false*/, int numEnergyEquations/*=0*/>
struct MPNCEnergyIndices
{
    static_assert(((numEnergyEquations<1) and  not enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");

    static_assert( (numEnergyEquations < 1) ,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");
public:
    /*!
     * \brief This module does not define any primary variables in the
     *        isothermal case.
     */
    static const unsigned int numPrimaryVars = 0;
};

/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitIndices
 * \brief The indices required for the energy equation.
 */
template <int PVOffset>
struct MPNCEnergyIndices<PVOffset, /*isNonIsothermal=*/true, /*numEnergyEquations=*/ 1 >
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

}

#endif
