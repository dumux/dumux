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
 * \brief Defines the indices required for the 2p1cni model.
 *
 *Important note: The 2p1c model requires the use of the non-isothermal extension found in dumux/implicit/nonisothermal
 */
#ifndef DUMUX_2P1C_INDICES_HH
#define DUMUX_2P1C_INDICES_HH

namespace Dumux
{

/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitIndices
 * \brief Enumerates the formulations which the two-phase n-component model accepts.
 */
struct TwoPOneCFormulation
{
    enum {
        pnsw,
        pwsn
        };
};

/*!
 * \ingroup TwoPOneCModel
 * \ingroup ImplicitIndices
 * \brief The indices for the isothermal 2p1cni model.
 *
 * \tparam formulation The formulation, only pgSwSn
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
class TwoPOneCIndices
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< index of the wetting liquid phase
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< index of the gas phase

    // present phases (-> 'pseudo' primary variable)
    static const int twoPhases = 1; //!< All three phases are present
    static const int wPhaseOnly = 2; //!< Only the water phase is present
    static const int nPhaseOnly = 3; //!< Only gas phase is present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for phase pressure in a solution vector
    static const int switch1Idx = PVOffset + 1; //!< Index of saturation or temperature

    // equation indices
    static const int conti0EqIdx = PVOffset    + 0; //wCompIdx; //!< Index of the mass conservation equation for the water component
    static const int energyEqIdx = PVOffset + 1;//! The index for energy in equation vectors.
};

}

#endif
