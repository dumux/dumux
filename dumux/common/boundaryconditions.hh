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
 * \ingroup Common
 * \brief Definition of boundary condition types, extend if necessary
 * \todo can this be removed for the sake of boundarytypes.hh?
 */
#ifndef DUMUX_BOUNDARYCONDITIONS_HH
#define DUMUX_BOUNDARYCONDITIONS_HH

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Define a class containing boundary condition flags
 */
struct BoundaryConditions
{
    /** \brief These values are ordered according to precedence */
    enum Flags {
        outflow = 0, //!< An outflow boundary
        neumann = 1, //!< Neumann boundary
        process = 2, //!< Processor boundary
        dirichlet = 3 //!< Dirichlet boundary
    };
};

} // end namespace Dumux

#endif
