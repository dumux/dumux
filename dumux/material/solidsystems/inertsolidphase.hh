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
 * \ingroup SolidSystems
 * \brief A solid phase consisting of a single inert solid component.
 */
#ifndef DUMUX_SOLIDSYSTEMS_INERT_SOLID_PHASE_HH
#define DUMUX_SOLIDSYSTEMS_INERT_SOLID_PHASE_HH

#include <string>
#include <dune/common/exceptions.hh>
#include <dumux/material/solidsystems/1csolid.hh>
#warning "This header is deprecated and will be removed after release 3.2. Use solidsystem.hh and SolidSystems::OneCSolid."

namespace Dumux {
namespace SolidSystems {

/*!
 * \ingroup SolidSystems
 * \brief A solid phase consisting of a single inert solid component.
 * \note a solid is considered inert if it can't dissolve in a liquid and
 *       and can't increase its mass by precipitation from a fluid phase.
 */
template <class Scalar, class ComponentT>
using InertSolidPhase = OneCSolid<Scalar, ComponentT, /*isInert=*/true>;

} // end namespace SolidSystems
} // end namespace Dumux

#endif
