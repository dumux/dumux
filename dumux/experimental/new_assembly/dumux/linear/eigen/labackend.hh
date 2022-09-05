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
 * \brief Default linear algebra backend for Eigen linear solvers.
 */
#ifndef DUMUX_EXPERIMENTAL_LINEAR_EIGEN_LA_BACKEND_HH
#define DUMUX_EXPERIMENTAL_LINEAR_EIGEN_LA_BACKEND_HH

#include <config.h>

#if HAVE_EIGEN3

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/indexstrategies.hh>

// include helper functions for linear algebra backends
#include <dumux/experimental/new_assembly/dumux/linear/labackend.hh>

// include traits for eigen linear systems
#include <dumux/experimental/new_assembly/dumux/linear/eigen/systemtraits.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

struct EigenTriplet
{
    std::size_t rowIndex;
    std::size_t colIndex;

    std::size_t row() const { return rowIndex; }
    std::size_t col() const { return colIndex; }
    double value() const { return 0.0; }
};

std::vector<EigenTriplet> asEigenTriplets(const MatrixPattern& pattern)
{
    std::vector<EigenTriplet> triplets;
    triplets.reserve(pattern.nnz());
    for (const auto& index : indices(pattern))
        triplets.emplace_back(EigenTriplet{index.rowIndex, index.colIndex});
    return triplets;
}

} // namespace Detail
#endif // DOXYGEN

//! Set the sparsity pattern of an Eigen sparse matrix
template<typename T, int O, typename S>
void setSparsityPattern(Eigen::SparseMatrix<T, O, S>& m, const MatrixPattern& pattern)
{
    m.resize(pattern.numRows(), pattern.numCols());
    const auto triplets = Detail::asEigenTriplets(pattern);
    m.setFromTriplets(triplets.begin(), triplets.end());
}


//! Backend for Eigen linear algebra types
template<Concepts::Arithmetic Scalar, int n>
struct EigenLinearAlgebraBackend
{
    using Matrix = Eigen::SparseMatrix<Scalar>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using IndexStrategy = FlatIndexStrategy<2>;

    static IndexStrategy makeIndexStrategy(std::size_t size)
    { return {BlockSizes{size, n}}; }

    static void setSparsityPattern(Matrix& m, const MatrixPattern& pattern)
    { Dumux::setSparsityPattern(m, asFlatPattern<n>(pattern)); }

    static void setSize(Vector& v, std::size_t size)
    { v.resize(size); }

};

} // namespace Dumux

#endif // HAVE_EIGEN3
#endif
