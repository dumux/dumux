// $Id$
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_BOUNDARYCONDITIONS_HH
#define DUMUX_BOUNDARYCONDITIONS_HH

/**
* @file
* @brief  Definition of boundary condition types, extend if necessary
* @author Peter Bastian
*/
namespace Dumux
{
/*!
 * \ingroup IMPETbc
 */
/**
* @brief Define a class containing boundary condition flags
*
*/

//! base Class that defines boundary condition flags
struct BoundaryConditions
{
    /** \brief These values are ordered according to precedence */
    enum Flags {
        couplingOutflow = -2, //!< An outflow boundary for coupled models
        couplingInflow = -1, //!< An inflow boundary for coupled models
        outflow = 0, //!< An outflow boundary
        neumann = 1, //!< Neumann boundary
        process = 2, //!< Processor boundary
        dirichlet = 3 //!< Dirichlet boundary
    };
};

/** @} */
}
#endif
