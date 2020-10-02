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
 * \ingroup Linear
 * \brief The specialized Dumux macro and tag for the ISTL registry to choose the
 *        solver and preconditioner at runtime
 */
#ifndef DUMUX_LINEAR_ISTL_SOLVER_REGISTRY_HH
#define DUMUX_LINEAR_ISTL_SOLVER_REGISTRY_HH

#include <dune/istl/common/registry.hh>

/*!
 * \brief Register a Dumux preconditioner
 *
 * Use this macro in namespace Dumux.
 * Example:
 * DUMUX_REGISTER_PRECONDITIONER("mypreconditioner", Dumux::MultiTypeBlockMatrixPreconditionerTag, Dune::defaultPreconditionerBlockLevelCreator<Dumux::MyPreconditioner, 1>());
 * Expicitly specifying the namespaces is required.
 * Set parameter Preconditioner.Type to "mypreconditioner" to use it through the factory.
 *
 * In the macro implementation, the final static_assert forces implementers
 * to put a semicolon after every DUMUX_REGISTER_PRECONDITIONER macro call (cf. example)
 * and avoids a compiler warning for an empty line semicolon at the same time
 */
#define DUMUX_REGISTER_PRECONDITIONER(name, tag, ...)                 \
} namespace Dune {                                               \
DUNE_REGISTRY_PUT(tag, name, __VA_ARGS__);  \
} namespace Dumux { \
static_assert(true, "Require semicolon after macro call")

namespace Dumux {
namespace {
struct MultiTypeBlockMatrixPreconditionerTag {};
struct MultiTypeBlockMatrixDirectSolverTag {};
} // end namespace
} // end namespace Dumux

#endif
