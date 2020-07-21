// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup MineralizationModel
 * \brief Defines the properties required for compositional porous medium flow
 *        models considering mineralization processes of one or more of the
 *        components.
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consist only of a storage and a source term:
 * \f$\frac{\partial ( \varrho_\lambda \phi_\lambda )} {\partial t} = q_\lambda\f$
 */

#ifndef DUMUX_MINERALIZATION_MODEL_HH
#define DUMUX_MINERALIZATION_MODEL_HH

#include <string>

namespace Dumux {

/*!
 * \ingroup MineralizationModel
 * \brief Specifies a number properties of
 *        models that consider mineralization proceses.
 *
 * \ţparam NonMinTraits traits class of the underlying model
 *                      not considering mineralization.
 * \tparam numPS number of solid phases to be considered.
 * \tparam numInertSP number of inert solid phases to be considered.
 */
template<class NonMinTraits, int numSC, int numInertSC>
struct MineralizationModelTraits : public NonMinTraits
{
    //! the number of mineral phases
    static constexpr int numSolidComps() { return numSC; }
     //! the number of inert mineral phases
    static constexpr int numInertSolidComps() { return numInertSC; }
    //! we additionally solve one equation per precipitating mineral phase
    static constexpr int numEq() { return NonMinTraits::numEq() + numSC - numInertSC; }
};
} // end namespace Dumux

#endif
