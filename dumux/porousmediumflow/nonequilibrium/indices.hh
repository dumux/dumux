// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NonEquilibriumModel
 * \brief The primary variable and equation indices for the MpNc model.
 */

#ifndef DUMUX_NONEQUILIBRIUM_INDICES_HH
#define DUMUX_NONEQUILIBRIUM_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NonEquilibriumModel
 * \brief The primary variable and equation indices for the MpNc model.
 */
template <class EquilibriumIndices, int numEnergyEqFluid, int numEnergyEqSolid, int numEq>
class NonEquilbriumIndices : public EquilibriumIndices
{
public:
    /*!
     * \brief Index for the temperature of the wetting phase in a vector of primary
     *        variables.
     */
    static constexpr unsigned int temperature0Idx = numEq - numEnergyEqFluid - numEnergyEqSolid;

    /*!
     * \brief Index for the temperature of the solid phase in a vector of primary
     *        variables.
     */
    static constexpr unsigned int temperatureSolidIdx = numEq - numEnergyEqSolid;
    /*!
     * \brief Compatibility with non kinetic models
     */
    static constexpr unsigned int temperatureIdx = temperature0Idx;
    /*!
     * \brief Equation index of the energy equation.
     */
    static constexpr unsigned int energyEq0Idx = numEq - numEnergyEqFluid - numEnergyEqSolid;
    /*!
     * \brief Compatibility with non kinetic models
     */
    static constexpr unsigned int energyEqIdx = energyEq0Idx;

    /*!
     * \brief Equation index of the energy equation.
     */
    static constexpr unsigned int energyEqSolidIdx = numEq - numEnergyEqSolid;
};

} // end namespace Dumux

#endif
