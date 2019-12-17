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
 * \brief Defines the properties required for (immiscible) twophase sequential models.
 */
#ifndef DUMUX_MIMETICPROPERTIES2P_SEQUENTIAL_HH
#define DUMUX_MIMETICPROPERTIES2P_SEQUENTIAL_HH

//Dumux-includes
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/sequential/mimetic/properties.hh>
namespace Dumux {

////////////////////////////////
// Forward declarations
////////////////////////////////

////////////////////////////////
// Properties
////////////////////////////////
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the two-phase problems
// Create new type tags
namespace TTag {
struct MimeticPressureTwoP { using InheritsFrom = std::tuple<Mimetic, PressureTwoP>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
} // end namespace Properties
} // end namespace Dumux

#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressure.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/mimetic.hh>

namespace Dumux {
namespace Properties {
template<class TypeTag>
struct PressureModel<TypeTag, TTag::MimeticPressureTwoP> { using type = MimeticPressure2P<TypeTag>; };
template<class TypeTag>
struct LocalStiffness<TypeTag, TTag::MimeticPressureTwoP> { using type = MimeticTwoPLocalStiffness<TypeTag>; };
} // end namespace Properties
} // end namespace Dumux

#endif
