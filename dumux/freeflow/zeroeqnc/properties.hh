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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup BoxZeroEqncModel
 *
 * \file
 *
 * \brief Defines the properties required for the compositional ZeroEq box model.
 */

#ifndef DUMUX_ZEROEQNC_PROPERTIES_HH
#define DUMUX_ZEROEQNC_PROPERTIES_HH

#include <dumux/freeflow/stokesnc/properties.hh>
#include <dumux/freeflow/zeroeq/properties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the ZeroEqnc Problem
NEW_TYPE_TAG(BoxZeroEqnc, INHERITS_FROM(BoxZeroEq, BoxStokesnc));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(ZeroEqEddyDiffusivityModel); //!< Returns which Eddy Diffusivity Model is used
NEW_PROP_TAG(ZeroEqTurbulentSchmidtNumber); //!< Returns the used turbulent Schmidt number (Reynolds analogy)
}
}

#endif // DUMUX_ZEROEQNC_PROPERTIES_HH
