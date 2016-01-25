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
 * \brief The indices for the mass flow part of the compositional
 *        multi-phase model.
 */
#ifndef DUMUX_MPNC_MASS_INDICES_HH
#define DUMUX_MPNC_MASS_INDICES_HH

#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitIndices
 * \brief The indices required for conservation of mass.
 *
 * This is the specialization for the case without kinetic mass
 * transfer (i.e. assuming chemical equilibrium)
 */
template <int PVOffset,
          class TypeTag,
          bool enableKinetic /*=false*/>
class MPNCMassIndices
{
    static_assert(!enableKinetic,
                  "No kinetic mass transfer module included, "
                  "but kinetic mass transfer enabled.");

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

public:
    /*!
     * \brief This module defines one new primary variable.
     */
    static const unsigned int numPrimaryVars = numComponents;

    /*!
     * \brief Index for the fugacity of the first component in the
     *        first phase in a vector of primary variables.
     *
     * The next numComponents indices represent the remaining
     * fugacities:
     *
     *  fug0Idx + 0 = fugacity of component 0
     *  fug0Idx + 1 = fugacity of component 1
     *  ...
     *  fug0Idx + N - 1 = fugacity of component N
     */
    static const unsigned int fug0Idx = PVOffset + 0;

    /*!
     * \brief Equation index of the mass conservation equation for the
     *        first component.
     *
     * The next numComponents indices represent the equations of all
     * components in all phases:
     *
     *  conti00EqIdx + 0 = continuity of component 0
     *  conti00EqIdx + 1 = continuity of component 1
     *  ...
     *  conti00EqIdx + N - 1 = continuity of component N
     */
    static const unsigned int conti0EqIdx = PVOffset + 0;
};

}

#endif
