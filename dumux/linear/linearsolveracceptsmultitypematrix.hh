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
 * \file
 * \ingroup Linear
 * \brief Trait checking if linear solvers accept Dune::MultiTypeBlockMatrix or we need to convert the matrix
 */
#ifndef DUMUX_LINEAR_SOLVER_ACCEPTS_MULTITYPEMATRIX_HH
#define DUMUX_LINEAR_SOLVER_ACCEPTS_MULTITYPEMATRIX_HH

#include <dumux/linear/seqsolverbackend.hh>

namespace Dumux {

//! The default
template<typename TypeTag, typename LinearSolver>
struct LinearSolverAcceptsMultiTypeMatrixImpl
{ static constexpr bool value = true; };

/*!
 * \ingroup Linear
 * \brief Trait checking if linear solvers accept Dune::MultiTypeBlockMatrix or we need to convert the matrix
 */
template <typename TypeTag>
using LinearSolverAcceptsMultiTypeMatrix = LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, typename GET_PROP_TYPE(TypeTag, LinearSolver)>;

//! Solvers that don't accept multi-type matrices
//! Those are all with ILU preconditioner that doesn't support the additional block level
//! And the direct solvers that have BCRS Matrix hardcoded

template<typename TypeTag>
struct LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, ILUnBiCGSTABBackend>
{ static constexpr bool value = false; };

template<typename TypeTag>
struct LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, ILUnCGBackend>
{ static constexpr bool value = false; };

template<typename TypeTag>
struct LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, ILU0BiCGSTABBackend>
{ static constexpr bool value = false; };

template<typename TypeTag>
struct LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, ILU0CGBackend>
{ static constexpr bool value = false; };

template<typename TypeTag>
struct LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, ILU0RestartedGMResBackend>
{ static constexpr bool value = false; };

template<typename TypeTag>
struct LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, ILUnRestartedGMResBackend>
{ static constexpr bool value = false; };

#if HAVE_SUPERLU
template<typename TypeTag>
struct LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, SuperLUBackend>
{ static constexpr bool value = false; };
#endif // HAVE_SUPERLU

#if HAVE_UMFPACK
template<typename TypeTag>
struct LinearSolverAcceptsMultiTypeMatrixImpl<TypeTag, UMFPackBackend>
{ static constexpr bool value = false; };
#endif // HAVE_UMFPACK

} // end namespace Dumux

#endif
