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
 * \ingroup Assembly
 * \brief An enum class to define the colors of elements and vertices required
 *        for partial Jacobian reassembly.
 */
#ifndef DUMUX_ENTITY_COLOR_HH
#define DUMUX_ENTITY_COLOR_HH

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The colors of elements and vertices required for partial
 *        Jacobian reassembly.
 */
enum class EntityColor {
    //! distance from last linearization is above the tolerance
    red,

    //! neighboring entity is red
    yellow,

    /*!
     * A yellow entity that has only non-green neighbor elements.
     *
     * This means that its relative error is below the tolerance,
     * but its defect can be linearized without any additional
     * cost.
     */
    orange,

    //! does not need to be reassembled
    green
};

} // end namespace Dumux

#endif
