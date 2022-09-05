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
 * \brief Traits related to solving linear systems with Eigen.
 */
#ifndef DUMUX_EXPERIMENTAL_LINEAR_EIGEN_SYSTEM_TRAITS_HH
#define DUMUX_EXPERIMENTAL_LINEAR_EIGEN_SYSTEM_TRAITS_HH

#include <config.h>

#if HAVE_EIGEN3

#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>

#include <dumux/experimental/new_assembly/dumux/linear/systemtraits.hh>
#include <dumux/experimental/new_assembly/dumux/linear/eigen/detail.hh>

namespace Dumux::LinearSystem::Traits {

namespace Eigen_ {

template<typename Matrix, Concepts::Arithmetic V>
struct MatFill
{
    static void fill(Matrix& m, const V& value)
    { m.fill(value); }
};

template<typename SparseMatrix, Concepts::Arithmetic V>
struct SparseMatFill
{
    static void fill(SparseMatrix& m, const V& value)
    { m.coeffs() = value; }
};

template<typename Y, typename X, Concepts::Arithmetic A>
struct Axpy
{
    static void axpy(Y& y, const A& a, const X& x)
    { y += a*x; }
};

template<typename Vector, Concepts::MultiIndex I> requires(I::size() == 1)
struct VecAccess
{
    static decltype(auto) get(const Vector& v, const I& index)
    { return v(index.template getAs<std::size_t, 0>()); }

    template<typename V>
    static void set(Vector& v, const I& index, const V& value)
    { v(index.template getAs<std::size_t, 0>()) = value; }
};

template<typename Matrix, Concepts::MultiIndex I, Concepts::MultiIndex J> requires(
    I::size() == 1 and J::size() == 1)
struct MatAccess
{
    static decltype(auto) get(const Matrix& m, const I& i, const J& j)
    {
        return m(i.template getAs<std::size_t, 0>(),
                 j.template getAs<std::size_t, 0>());
    }

    template<typename V>
    static void set(Matrix& m, const I& i, const J& j, const V& value)
    {
        m(i.template getAs<std::size_t, 0>(),
          j.template getAs<std::size_t, 0>()) = value;
    }
};

template<typename Matrix, Concepts::MultiIndex I, Concepts::MultiIndex J> requires(
    I::size() == 1 and J::size() == 1)
struct SparseMatAccess
{
    static decltype(auto) get(const Matrix& m, const I& i, const J& j)
    {
        return m.coeff(i.template getAs<std::size_t, 0>(),
                       j.template getAs<std::size_t, 0>());
    }

    template<typename V>
    static void set(Matrix& m, const I& i, const J& j, const V& value)
    {
         m.coeffRef(i.template getAs<std::size_t, 0>(),
                    j.template getAs<std::size_t, 0>()) = value;
    }
};

} // namespace Eigen_

template<typename T, int r, int c, int o, int rc, int cc>
struct Scalar<Eigen::Matrix<T, r, c, o, rc, cc>>
: public std::type_identity<T> {};

template<typename T, int O, typename S>
struct Scalar<Eigen::SparseMatrix<T, O, S>>
: public std::type_identity<T> {};

template<typename T, int r, int c, int O, int rc, int cc, Concepts::Arithmetic V>
struct Fill<Eigen::Matrix<T, r, c, O, rc, cc>, V>
: public Eigen_::MatFill<Eigen::Matrix<T, r, c, O, rc, cc>, V> {};

template<typename T, int O, typename S, Concepts::Arithmetic V>
struct Fill<Eigen::SparseMatrix<T, O, S>, V>
: public Eigen_::SparseMatFill<Eigen::SparseMatrix<T, O, S>, V> {};

template<Detail::EigenVector V, Detail::EigenVector X, Concepts::Arithmetic A>
struct Axpy<V, X, A> : public Eigen_::Axpy<V, X, A> {};

template<Detail::EigenMatrix V, Detail::EigenMatrix X, Concepts::Arithmetic A>
struct Axpy<V, X, A> : public Eigen_::Axpy<V, X, A> {};

template<Detail::EigenVector V, Concepts::MultiIndex I>
struct VectorAccess<V, I> : public Eigen_::VecAccess<V, I> {};

template<typename T, int r, int c, int O, int rc, int cc, Concepts::MultiIndex I, Concepts::MultiIndex J>
struct MatrixAccess<Eigen::Matrix<T, r, c, O, rc, cc>, I, J>
: public Eigen_::MatAccess<Eigen::Matrix<T, r, c, O, rc, cc>, I, J>
{};

template<typename T, int O, typename S, Concepts::MultiIndex I, Concepts::MultiIndex J>
struct MatrixAccess<Eigen::SparseMatrix<T, O, S>, I, J>
: public Eigen_::SparseMatAccess<Eigen::SparseMatrix<T, O, S>, I, J>
{};

} // namespace Dumux::LinearSystem::Traits

#endif // HAVE_EIGEN3
#endif
