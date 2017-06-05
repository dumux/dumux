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

//! The type of the linear solver to be used
NEW_PROP_TAG(LinearSolver);

/*!
 * \brief Specifies the verbosity of the linear solver
 *
 * By default it is 0, i.e. it doesn't print anything. Setting this
 * property to 1 prints aggregated convergence rates, 2 prints the
 * convergence rate of every iteration of the scheme.
 */
NEW_PROP_TAG(LinearSolverVerbosity);

//! target reduction of the initial residual
NEW_PROP_TAG(LinearSolverResidualReduction);

//! maximum number of iterations of solver
NEW_PROP_TAG(LinearSolverMaxIterations);

//! relaxation parameter for the preconditioner
NEW_PROP_TAG(LinearSolverPreconditionerRelaxation);

//! number of preconditioner iterations per solver iteration
NEW_PROP_TAG(LinearSolverPreconditionerIterations);

//! Verbosity level of the preconditioner
NEW_PROP_TAG(LinearSolverPreconditionerVerbosity);

//! Block level depth for the preconditioner
// Set this to more than one if the matrix to solve is nested multiple times
// e.g. for Dune::MultiTypeBlockMatrix'es.
NEW_PROP_TAG(LinearSolverPreconditionerBlockLevel);

//! restart parameter for GMRes
NEW_PROP_TAG(LinearSolverGMResRestart);

//! Size of the matrix/vector blocks
/*!
 * The number of different types of equations which build the system of equations to solve
 * can differ from the number of equations given by the mathematical/physical model (e.g. IMPES).
 * Thus, the block size does not have to be equal to NumEq.
 * (Especially important for the SuperLU solver!)
 */
NEW_PROP_TAG(LinearSolverBlockSize);

SET_INT_PROP(LinearSolverTypeTag, LinearSolverVerbosity, 0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_SCALAR_PROP(LinearSolverTypeTag, LinearSolverPreconditionerRelaxation, 1.0);

//! set the preconditioner iterations to 1 by default
SET_INT_PROP(LinearSolverTypeTag, LinearSolverPreconditionerIterations, 1);

//! set the preconditioner verbosity level to 0 by default
SET_INT_PROP(LinearSolverTypeTag, LinearSolverPreconditionerVerbosity, 0);

//! set the block level to 1, suitable for e.g. a simple Dune::BCRSMatrix.
SET_INT_PROP(LinearSolverTypeTag, LinearSolverPreconditionerBlockLevel, 1);

//! set the GMRes restart parameter to 10 by default
SET_INT_PROP(LinearSolverTypeTag, LinearSolverGMResRestart, 10);

} // namespace Properties
} // namespace Dumux

#endif
