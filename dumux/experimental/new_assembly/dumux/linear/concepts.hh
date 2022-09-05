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
 * \brief Concepts related to linear solvers
 */
#ifndef DUMUX_LINEAR_CONCEPTS_HH
#define DUMUX_LINEAR_CONCEPTS_HH

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/matrixpattern.hh>

namespace Dumux::Concepts {

template<typename T>
concept LinearAlgebraBackend = requires (const T& t) {
    typename T::Vector;
    typename T::Matrix;
    typename T::IndexStrategy;

    { t.makeIndexStrategy(std::size_t{}) };
    { t.setSparsityPattern(std::declval<typename T::Matrix&>(), MatrixPattern{}) };
    { t.setSize(std::declval<typename T::Vector&>(), std::size_t{}) };
};

template<typename Solver,
         typename Domain,
         typename Range,
         typename LinearOperator>
concept LinearSolver
    = requires(Solver& solver,
               const Solver& constSolver,
               Domain& domain,
               const Range& range,
               const LinearOperator& op) {
        { solver.setLinearOperator(op) };
        { constSolver.solve(domain, range) };
        { constSolver.scalarProduct(range, range) } -> Concepts::Arithmetic;
};

} // namespace Dumux::Concepts

#endif
