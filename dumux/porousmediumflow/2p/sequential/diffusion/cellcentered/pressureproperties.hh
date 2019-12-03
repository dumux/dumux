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
 * \brief Defines the properties required for finite volume pressure models in a two-phase sequential model.
 */
#ifndef DUMUX_FVPRESSUREPORPERTIES2P_SEQUENTIAL_HH
#define DUMUX_FVPRESSUREPORPERTIES2P_SEQUENTIAL_HH

//Dumux-includes
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>

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

//! The type tag for two-phase problems using a standard finite volume model
// Create new type tags
namespace TTag {
struct FVPressureTwoP { using InheritsFrom = std::tuple<PressureTwoP>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

} // end namespace Properties
} // end namespace Dumux

#include "velocity.hh"
#include "pressure.hh"

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
//! Set velocity reconstruction implementation standard cell centered finite volume schemes as default
template<class TypeTag>
struct Velocity<TypeTag, TTag:: FVPressureTwoP> { using type = FVVelocity2P<TypeTag> ; };
//! Set finite volume implementation of the two-phase pressure equation as default pressure model
template<class TypeTag>
struct PressureModel<TypeTag, TTag::FVPressureTwoP> { using type = FVPressure2P<TypeTag>; };
//! Allow assembling algorithm for the pressure matrix to assemble only from one side of a cell-cell interface
template<class TypeTag>
struct VisitFacesOnlyOnce<TypeTag, TTag::FVPressureTwoP> { static constexpr bool value = true; };

} // end namespace Properties
} // end namespace Dumux

#endif
