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
 * \brief Properties for the MPFA-O method.
 */
#ifndef DUMUX_FVMPFAO2DPROPERTIES2P_HH
#define DUMUX_FVMPFAO2DPROPERTIES2P_HH

// dumux environment
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>

namespace Dumux {
namespace Properties {
// Create new type tags
namespace TTag {
struct FvMpfaO2dPressureTwoP { using InheritsFrom = std::tuple<MPFAProperties, PressureTwoP>; };
} // end namespace TTag
} // end namespace Properties
} // end namespace Dumux

#include "2dpressurevelocity.hh"
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/velocityintransport.hh>

namespace Dumux {
namespace Properties {
template<class TypeTag>
struct PressureModel<TypeTag, TTag::FvMpfaO2dPressureTwoP> { using type = FvMpfaO2dPressureVelocity2p<TypeTag>; };
//! Set velocity reconstruction implementation standard cell centered finite volume schemes as default
template<class TypeTag>
struct Velocity<TypeTag, TTag:: FvMpfaO2dPressureTwoP> { using type = FvMpfaVelocityInTransport<TypeTag> ; };
} // end namespace Properties
} // end namespace Dumux
#endif
