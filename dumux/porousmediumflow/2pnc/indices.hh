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
 * \ingroup TwoPNCModel
 * \brief Defines the indices required for the two-phase n-component
 *        fully implicit model.
 */
#ifndef DUMUX_2PNC_INDICES_HH
#define DUMUX_2PNC_INDICES_HH

#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNCModel
 * \brief Enumerates the formulations which the two-phase n-component model accepts.
 */
struct TwoPNCFormulation
{
    enum {
        pnsw,
        pwsn
        };
};

/*!
 * \ingroup TwoPNCModel
 * \brief The indices for the isothermal two-phase n-component model.
 *
 * \tparam TypeTag The problem Type Tag
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
class TwoPNCIndices
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx;    //!< index of the wetting phase
    static const int nPhaseIdx = FluidSystem::nPhaseIdx;    //!< index of the non-wetting phase

    // Component indices
    static const int wCompIdx = FluidSystem::wCompIdx;      //!< index of the primary component of the wetting phase
    static const int nCompIdx = FluidSystem::nCompIdx;      //!< index of the primary component of the non-wetting phase

    // present phases (-> 'pseudo' primary variable)
    static const int wPhaseOnly = 1; //!< Only the non-wetting phase is present
    static const int nPhaseOnly = 2; //!< Only the wetting phase is present
    static const int bothPhases = 3; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0;    //!< index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx = PVOffset + 1;      //!< index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    // equation indices
    static const int conti0EqIdx = PVOffset + 0;                        //!< Reference index for mass conservation equations.
    static const int contiWEqIdx = conti0EqIdx + FluidSystem::wCompIdx; //!< index of the mass conservation equation for the wetting phase major component
    static const int contiNEqIdx = conti0EqIdx + FluidSystem::nCompIdx; //!< index of the mass conservation equation for the non-wetting phase major component
};

// \}

}

#endif
