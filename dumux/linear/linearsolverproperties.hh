// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Properties
 * \ingroup Linear
 * \file
 *
 * \brief Defines a type tag and some fundamental properties for
 *        linear solvers
 */
#ifndef DUMUX_LINEAR_SOLVER_PROPERTIES_HH
#define DUMUX_LINEAR_SOLVER_PROPERTIES_HH

#include <dumux/common/basicproperties.hh>

namespace Dumux
{
namespace Properties
{
//! Linear solver type tag for all models.
NEW_TYPE_TAG(LinearSolverTypeTag);

///////////////////////////////////
// Property tag declarations:
///////////////////////////////////

//! Block level depth for the preconditioner
// Set this to more than one if the matrix to solve is nested multiple times
// e.g. for Dune::MultiTypeBlockMatrix'es.
NEW_PROP_TAG(LinearSolverPreconditionerBlockLevel);

//! Size of the matrix/vector blocks
/*!
 * The number of different types of equations which build the system of equations to solve
 * can differ from the number of equations given by the mathematical/physical model (e.g. IMPES).
 * Thus, the block size does not have to be equal to NumEq.
 * (Especially important for the SuperLU solver!)
 */
NEW_PROP_TAG(LinearSolverBlockSize);

///////////////////////////////////
// Default values for properties:
///////////////////////////////////

//! set the block level to 1, suitable for e.g. a simple Dune::BCRSMatrix.
SET_INT_PROP(LinearSolverTypeTag, LinearSolverPreconditionerBlockLevel, 1);

//! set the block size to 1 as default
SET_INT_PROP(LinearSolverTypeTag, LinearSolverBlockSize, 1);

} // namespace Properties
} // namespace Dumux

#endif
