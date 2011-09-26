/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief The indices for the non-isothermal part of the compositional
 *        multi-phase model.
 */
#ifndef DUMUX_MPNC_INDICES_ENERGY_HH
#define DUMUX_MPNC_INDICES_ENERGY_HH

namespace Dumux
{
/*!
 * \brief The indices for the energy equation.
 *
 * This is a dummy class for the isothermal case.
 */
template <int PVOffset, bool enableEnergy/*=false*/, bool kineticEnergyTransfer/*=false*/>
struct MPNCEnergyIndices
{
    static_assert(!(kineticEnergyTransfer && !enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(!kineticEnergyTransfer,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");
public:
    /*!
     * \brief This module does not define any primary variables in the
     *        isothermal case.
     */
    static const int NumPrimaryVars = 0;
};

/*!
 * \brief The indices required for the energy equation.
 */
template <int PVOffset>
struct MPNCEnergyIndices<PVOffset, /*isNonIsothermal=*/true, /*kineticEnergyTransfer*/false>
{
public:
    /*!
     * \brief This module defines one new primary variable.
     */
    static const int NumPrimaryVars = 1;

    /*!
     * \brief Index for the temperature in a vector of primary
     *        variables.
     */
    static const int temperatureIdx = PVOffset + 0;
    /*!
     * \brief Equation index of the energy equation.
     */
    static const int energyEqIdx = PVOffset + 0;
};

}

#endif
