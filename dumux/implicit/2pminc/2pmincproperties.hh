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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TwoPTwoCModel
 *
 * \brief Defines the properties required for the isothermal two-phase,
 * two-component BOX model.
 */
#ifndef DUMUX_2PMINC_PROPERTIES_HH
#define DUMUX_2PMINC_PROPERTIES_HH

#include <dumux/implicit/2p/2pproperties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the two-phase minc model with BOX discretization
NEW_TYPE_TAG(BoxTwoPMinc, INHERITS_FROM(BoxTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumContinua);   //!< Number of continua
/*! How the nested continua are interacting
 *
 * 0: constant volume fractions
 * 1: equally distanced nested continua
 * 2: DFM generated volume fractions
 *
 */
NEW_PROP_TAG(ProblemInteractingContinuaType);
}
}

#endif
