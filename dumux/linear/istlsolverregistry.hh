// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * Explicitly specifying the namespaces is required.
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
