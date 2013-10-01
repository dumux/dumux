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
 * \brief Defines the indices required for the 2p2c fully implicit model.
 */
#ifndef DUMUX_2P2C_INDICES_HH
#define DUMUX_2P2C_INDICES_HH

#include "2p2cproperties.hh"

namespace Dumux
{
// \{

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitIndices
 * \brief Enumerates the formulations which the 2p2c model accepts.
 */
struct TwoPTwoCFormulation
{
    static const int pwsn = 0; //!< pw and sn as primary variables
    DUNE_DEPRECATED_MSG("use pwsn (uncapitalized 'S') instead") 
    static const int pwSn = pwsn; //!< \deprecated

    static const int pnsw = 1; //!< pn and sw as primary variables
    DUNE_DEPRECATED_MSG("use pnsw (uncapitalized 'S') instead") 
    static const int pnSw = pnsw; //!< \deprecated
};

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitIndices
 * \brief The indices for the isothermal 2p2c model.
 *
 * \tparam formulation The formulation, either pwsn or pnsw.
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag,
          int formulation = TwoPTwoCFormulation::pwsn,
          int PVOffset = 0>
class TwoPTwoCIndices
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the wetting phase
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the non-wetting phase

    // Component indices
    static const int wCompIdx = FluidSystem::wCompIdx; //!< Index of the primary component of the wetting phase
    static const int nCompIdx = FluidSystem::nCompIdx; //!< Index of the primary component of the non-wetting phase

    // present phases (-> 'pseudo' primary variable)
    static const int wPhaseOnly = 1; //!< Only the wetting phase is present
    static const int nPhaseOnly = 0; //!< Only the non-wetting phase is present
    static const int bothPhases = 2; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on the formulation) in a solution vector
    static const int switchIdx = PVOffset + 1; //!< Index of either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int pwIdx = pressureIdx; //!< Index for wetting phase pressure in a solution vector
    static const int snOrXIdx = switchIdx; //!< Index of either the saturation of the non-wetting phase or the mass fraction secondary component in the only phase

    DUNE_DEPRECATED_MSG("use snOrXIdx (uncapitalized 's') instead") 
    static const int SnOrXIdx = snOrXIdx; //!< \deprecated

    // equation indices
    static const int conti0EqIdx = PVOffset; //!< Index of the mass conservation equation for the first component
    static const int contiWEqIdx = conti0EqIdx + wCompIdx; //!< Index of the mass conservation equation for the primary component of the wetting phase
    static const int contiNEqIdx = conti0EqIdx + nCompIdx; //!< Index of the mass conservation equation for the primary component of the non-wetting phase
};

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitIndices
 * \brief The indices for the isothermal 2p2c model in the pn-sw
 *        formulation.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class TwoPTwoCIndices<TypeTag, TwoPTwoCFormulation::pnsw, PVOffset>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the wetting phase
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the non-wetting phase

    // Component indices
    static const int wCompIdx = FluidSystem::wCompIdx; //!< Index of the primary component of the wetting phase
    static const int nCompIdx = FluidSystem::nCompIdx; //!< Index of the primary component of the non-wetting phase

    // present phases (-> 'pseudo' primary variable)
    static const int wPhaseOnly = 1; //!< Only the wetting phase is present
    static const int nPhaseOnly = 2; //!< Only the non-wetting phase is present
    static const int bothPhases = 3; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on the formulation) in a solution vector
    static const int switchIdx = PVOffset + 1; //!< Index of either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int pnIdx = pressureIdx; //!< Index for non-wetting phase pressure in a solution vector
    static const int swOrXIdx = switchIdx; //!< Index of either the saturation of the liquid phase or the mass fraction of the secondary component in the only phase

    DUNE_DEPRECATED_MSG("use swOrXIdx (uncapitalized 's') instead") 
    static const int SwOrXIdx = swOrXIdx; //!< \deprecated

    // Equation indices
    static const int conti0EqIdx = PVOffset; //!< Index of the mass conservation equation for the first component
    static const int contiWEqIdx = conti0EqIdx + wCompIdx; //!< Index of the mass conservation equation for the primary component of the wetting phase
    static const int contiNEqIdx = conti0EqIdx + nCompIdx; //!< Index of the mass conservation equation for the primary component of the non-wetting phase
};

// \}

}

#endif
