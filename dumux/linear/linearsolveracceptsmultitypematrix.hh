// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Trait checking if linear solvers accept Dune::MultiTypeBlockMatrix or we need to convert the matrix
 */
#ifndef DUMUX_LINEAR_SOLVER_ACCEPTS_MULTITYPEMATRIX_HH
#define DUMUX_LINEAR_SOLVER_ACCEPTS_MULTITYPEMATRIX_HH

#include <dumux/linear/seqsolverbackend.hh>
#ifndef DUMUX_SUPPRESS_LINEAR_SOLVER_ACCEPTS_MULTITYPEMATRIX_WARNING
#warning "This header is deprecated and will be removed after release 3.7."
#endif

// suppress all secondary deprecation warning from this deprecated file
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

namespace Dumux {

//! The default
template<class LinearSolver>
struct linearSolverAcceptsMultiTypeMatrix : public std::true_type {};

//! Solvers that don't accept multi-type matrices
//! Those are all with ILU preconditioner that doesn't support the additional block level
//! And the direct solvers that have BCRS Matrix hardcoded

class ILUnBiCGSTABBackend;

template<>
struct linearSolverAcceptsMultiTypeMatrix<ILUnBiCGSTABBackend> : public std::false_type {};

class ILUnCGBackend;

template<>
struct linearSolverAcceptsMultiTypeMatrix<ILUnCGBackend> : public std::false_type {};

class ILU0BiCGSTABBackend;

template<>
struct linearSolverAcceptsMultiTypeMatrix<ILU0BiCGSTABBackend> : public std::false_type {};

class ILU0CGBackend;

template<>
struct linearSolverAcceptsMultiTypeMatrix<ILU0CGBackend> : public std::false_type {};

class ILU0RestartedGMResBackend;

template<>
struct linearSolverAcceptsMultiTypeMatrix<ILU0RestartedGMResBackend> : public std::false_type {};

class ILUnRestartedGMResBackend;

template<>
struct linearSolverAcceptsMultiTypeMatrix<ILUnRestartedGMResBackend> : public std::false_type {};

#if HAVE_SUPERLU
class SuperLUBackend;
template<>
struct linearSolverAcceptsMultiTypeMatrix<SuperLUBackend> : public std::false_type {};
#endif // HAVE_SUPERLU

#if HAVE_UMFPACK
class UMFPackBackend;
template<>
struct linearSolverAcceptsMultiTypeMatrix<UMFPackBackend> : public std::false_type {};
#endif // HAVE_UMFPACK

} // end namespace Dumux

#pragma GCC diagnostic pop

#endif
