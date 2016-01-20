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
 * \brief The indices for the kinetic mass transfer module of the
 *        compositional multi-phase model.
 */
#ifndef DUMUX_MPNC_MASS_INDICES_KINETIC_HH
#define DUMUX_MPNC_MASS_INDICES_KINETIC_HH

#include <dumux/porousmediumflow/mpnc/implicit/indices.hh>

namespace Dumux
{
/*!
 * \brief The indices required for conservation of mass.
 *
 * This is the specialization considering kinetic mass transfer
 * (i.e. not assuming local chemical equilibrium)
 */
template <int PVOffset, class TypeTag>
class MPNCMassIndices<PVOffset, TypeTag, /*enableKinetic=*/true>
{
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
public:
    /*!
     * \brief This module defines one new primary variable.
     */
    static const unsigned int numPrimaryVars = numPhases*numComponents;

    /*!
     * \brief Index for the mole fraction of the first component in
     *        the first phase in a vector of primary variables.
     *
     * The next numPhases*numComponents indices represent the
     * compositions of all phases:
     *
     *  moleFrac00Idx + 0 = mole fraction of component 0 in phase 0
     *  moleFrac00Idx + 1 = mole fraction of component 1 in phase 0
     *  ...
     *  moleFrac00Idx + N - 1 = mole fraction of component N in phase 0
     *  moleFrac00Idx + N = mole fraction of component 0 in phase 1
     *  ...
     */
    static const unsigned int moleFrac00Idx = PVOffset + 0;

    /*!
     * \brief Equation index of the mass conservation equation for the
     *        first component in the first phase.
     *
     * The next numPhases*numComponents indices represent the
     * continuity equations of all components in all phases:
     *
     *  conti00EqIdx + 0 = continuity of component 0 in phase 0
     *  conti00EqIdx + 1 = continuity of component 1 in phase 0
     *  ...
     *  conti00EqIdx + N - 1 = continuity of component N in phase 0
     *  conti00EqIdx + N = continuity of component 0 in phase 1
     *  ...
     */
    static const unsigned int conti0EqIdx = PVOffset + 0;
};

}

#endif
