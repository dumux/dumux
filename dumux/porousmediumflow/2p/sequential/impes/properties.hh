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
 * \ingroup SequentialTwoPModel
 * \brief Properties related to the sequential IMPES algorithms.
 */
#ifndef DUMUX_IMPES2P_PROPERTIES_HH
#define DUMUX_IMPES2P_PROPERTIES_HH

#include <dumux/porousmediumflow/sequential/impetproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>

namespace Dumux {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! TypeTag for the two-phase IMPES scheme
// Create new type tags
namespace TTag {
struct IMPESTwoP { using InheritsFrom = std::tuple<SequentialTwoP, IMPET>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
} // end namespace Properties
} // end namespace Dumux

#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>

namespace Dumux {
namespace Properties {
} // end namespace Properties
} // end namespace Dumux

#endif
