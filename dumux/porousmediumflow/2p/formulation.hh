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
 * \ingroup TwoPModel
 * \brief Defines an enumeration for the formulations accepted by the two-phase model.
 */

#ifndef DUMUX_2P_FORMULATION_INDICES_HH
#define DUMUX_2P_FORMULATION_INDICES_HH

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Enumerates the formulations which the two-phase model accepts.
 */
enum class TwoPFormulation
{
    p0s1, //!< first phase pressure and second phase saturation as primary variables
    p1s0  //!< first phase saturation and second phase pressure as primary variables
};

} // end namespace Dumux

#endif
