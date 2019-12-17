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
#ifndef DUMUX_FVPRESSUREPORPERTIES2P_ADAPTIVE_SEQUENTIAL_HH
#define DUMUX_FVPRESSUREPORPERTIES2P_ADAPTIVE_SEQUENTIAL_HH

//Dumux-includes
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

//! The type tag for two-phase problems using a grid-adaptive finite volume model
// Create new type tags
namespace TTag {
struct FVPressureTwoPAdaptive { using InheritsFrom = std::tuple<PressureTwoP>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

} // end namespace Properties
} // end namespace Dumux

#include "velocityadaptive.hh"
#include "pressureadaptive.hh"

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Set velocity reconstruction implementation for grid-adaptive cell centered finite volume schemes as default
template<class TypeTag>
struct Velocity<TypeTag, TTag:: FVPressureTwoPAdaptive> { using type = FVVelocity2PAdaptive<TypeTag> ; };
//! Set finite volume implementation of the two-phase pressure equation which allows hanging nodes as default pressure model
template<class TypeTag>
struct PressureModel<TypeTag, TTag::FVPressureTwoPAdaptive> { using type = FVPressure2PAdaptive<TypeTag>; };
//! Allow assembling algorithm for the pressure matrix to assemble only from one side of a cell-cell interface
template<class TypeTag>
struct VisitFacesOnlyOnce<TypeTag, TTag::FVPressureTwoPAdaptive> { static constexpr bool value = true; };
} // end namespace Properties
} // end namespace Dumnux

#endif
