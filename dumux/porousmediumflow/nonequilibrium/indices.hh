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
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief The primary variable and equation indices for the MpNc model.
 */
#ifndef DUMUX_NONEQUILIBRIUM_INDICES_HH
#define DUMUX_NONEQUILIBRIUM_INDICES_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief The primary variable and equation indices for the MpNc model.
 */
template <class EquilibriumIndices, class FluidSystem, int numEnergyEquationFluid, int numEnergyEquationSolid, int numEquationBalance, int BasePVOffset = 0>
class NonEquilbriumIndices: public EquilibriumIndices
{
public:
     enum { numPhases = FluidSystem::numPhases };
     enum { numEnergyEqFluid = numEnergyEquationFluid };
     enum { numEnergyEqSolid = numEnergyEquationSolid };
     /*! \todo Replacing the sum below with GET_PROP_VALUE(TypeTag, NumEq)
      *        yields a compilation error with clang, due to complex
      *        interdependencies of MPNC and NonEquilibrium type tags and
      *        Indices classes. This should be fixed.
      */
     static const unsigned int numEq = numEquationBalance
                                     + numEnergyEqFluid
                                     + numEnergyEqSolid;

    /*!
     * \brief Index for the temperature of the wetting phase in a vector of primary
     *        variables.
     */
    static constexpr unsigned int temperature0Idx = BasePVOffset + numEq - numEnergyEqFluid - numEnergyEqSolid;

    /*!
     * \brief Index for the temperature of the solid phase in a vector of primary
     *        variables.
     */
    static constexpr unsigned int temperatureSolidIdx = BasePVOffset + numEq - numEnergyEqSolid;
    /*!
     * \brief Compatibility with non kinetic models
     */
    static constexpr unsigned int temperatureIdx = temperature0Idx;
    /*!
     * \brief Equation index of the energy equation.
     */
    static constexpr unsigned int energyEq0Idx = BasePVOffset + numEq - numEnergyEqFluid - numEnergyEqSolid;
    /*!
     * \brief Compatibility with non kinetic models
     */
    static constexpr unsigned int energyEqIdx = energyEq0Idx;

    /*!
     * \brief Equation index of the energy equation.
     */
    static constexpr unsigned int energyEqSolidIdx = BasePVOffset + numEq - numEnergyEqSolid;
};

} // end namespace Dumux

#endif
