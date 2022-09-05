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
 * \ingroup LinearSystem
 * \brief Concepts & functions related to operations on matrices/vectors.
 */
#ifndef DUMUX_COMMON_LINEAR_SYSTEM_HH
#define DUMUX_COMMON_LINEAR_SYSTEM_HH

#include <utility>
#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>
#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>
#include <dumux/experimental/new_assembly/dumux/linear/systemtraits.hh>

namespace Dumux::LinearSystem {

//! Concept for types that are registered as linear system types
template<typename X>
concept Interoperable = isComplete<Traits::Scalar<std::decay_t<X>>>;


//! Concept for types for which the axpy operation is registered
template<typename Y, typename X, typename A>
concept ScaleAddable
    = isComplete<Traits::Axpy<Y, X, A>>
    and requires(Y& y, const A& a, const X& x) {
        { Traits::Axpy<Y, X, A>::axpy(y, a, x) };
    };


//! Concept for types for which the fill operation is registered
template<typename T, typename V>
concept Fillable
    = isComplete<Traits::Fill<T, V>>
    and requires(T& t, const V& v) {
        { Traits::Fill<T, V>::fill(t, v) };
    };


//! Concept for types for which the vector-access operation is registered
template<typename T, typename I>
concept VectorAccessible
    = isComplete<Traits::VectorAccess<T, I>>
    and requires(const T& t, const I& i) {
        { Traits::VectorAccess<T, I>::get(t, i) };
    };


//! Concept for types for which the vector-write-access operation is registered
template<typename T, typename I, typename Value>
concept VectorSettable
    = isComplete<Traits::VectorAccess<T, I>>
    and requires(T& t, const I& i, const Value& value) {
        { Traits::VectorAccess<T, I>::set(t, i, value) };
    };


//! Concept for types for which the matrix-access operation is registered
template<typename T, typename I, typename J>
concept MatrixAccessible
    = isComplete<Traits::MatrixAccess<T, I, J>>
    and requires(const T& t, const I& i, const J& j) {
        { Traits::MatrixAccess<T, I, J>::get(t, i, j) };
    };


//! Concept for types for which the matrix-write-access operation is registered
template<typename T, typename I, typename J, typename Value>
concept MatrixSettable
    = isComplete<Traits::MatrixAccess<T, I, J>>
    and requires(T& t, const I& i, const J& j, const Value& value) {
        { Traits::MatrixAccess<T, I, J>::set(t, i, j, value) };
    };


//! Convenience alias for the scalar type used by a matrix/vector
template<typename T> requires isComplete<Traits::Scalar<T>>
using ScalarType = typename Traits::Scalar<T>::type;


//! Convenience alias for the type returned from vector access
template<typename T, Concepts::MultiIndex I> requires VectorAccessible<T, I>
using VectorReferenceType = decltype(
    Traits::VectorAccess<T, I>::get(std::declval<const T&>(), std::declval<const I&>())
);


//! Convenience alias for the type returned from matrix access
template<typename T, Concepts::MultiIndex I, Concepts::MultiIndex J> requires MatrixAccessible<T, I, J>
using MatrixReferenceType = decltype(
    Traits::MatrixAccess<T, I, J>::get(std::declval<const T&>(), std::declval<const I&>(), std::declval<const J&>())
);


/*!
 * \ingroup LinearSystem
 * \brief Fill a vector or matrix with a scalar
 */
template<typename T, Concepts::Arithmetic V> requires(
    Fillable<T, V>)
inline constexpr void fill(T& t, const V& v)
{ Traits::Fill<T, V>::fill(t, v); }


/*!
 * \ingroup LinearSytem
 * \brief Add two vectors/matrices
 */
template<typename T, typename V> requires(
    ScaleAddable<T, V, ScalarType<T>>)
inline constexpr void add(T& t, const V& v)
{ Traits::Axpy<T, V, ScalarType<T>>::axpy(t, 1.0, v); }


/*!
 * \ingroup LinearSytem
 * \brief Subtract two vectors/matrices
 */
template<typename T, typename V> requires(
    ScaleAddable<T, V, ScalarType<T>>)
inline constexpr void subtract(T& t, const V& v)
{ Traits::Axpy<T, V, ScalarType<T>>::axpy(t, -1.0, v); }


/*!
 * \ingroup LinearSytem
 * \brief Scale a vector/matrix and add it to another one
 */
template<typename Y, typename X, Concepts::Arithmetic A> requires(
    ScaleAddable<Y, X, A>)
inline constexpr void axpy(Y& y, const A& a, const X& x)
{ Traits::Axpy<Y, X, A>::axpy(y, a, x); }


/*!
 * \ingroup LinearSytem
 * \brief Set the value in a vector
 */
template<typename T, Concepts::MultiIndex I, typename Value> requires(
    VectorSettable<T, I, Value>)
inline constexpr void set(T& t, const I& i, const Value& value)
{ Traits::VectorAccess<T, I>::set(t, i, value); }


/*!
 * \ingroup LinearSytem
 * \brief Get a value from a vector
 */
template<typename T, Concepts::MultiIndex I> requires(
    VectorAccessible<T, I>)
inline constexpr decltype(auto) get(const T& t, const I& i)
{ return Traits::VectorAccess<T, I>::get(t, i); }


/*!
 * \ingroup LinearSytem
 * \brief Add to a value of a vector
 */
template<typename T, Concepts::MultiIndex I, typename Value> requires(
    Concepts::Addable<VectorReferenceType<T, I>, Value>)
inline constexpr void add(T& t, const I& i, const Value& value)
{ set(t, i, get(t, i) + value); }


/*!
 * \ingroup LinearSytem
 * \brief Set the value in a matrix
 */
template<typename T, Concepts::MultiIndex I, Concepts::MultiIndex J, typename Value> requires(
    MatrixSettable<T, I, J, Value>)
inline constexpr void set(T& t, const I& i, const J& j, const Value& value)
{ Traits::MatrixAccess<T, I, J>::set(t, i, j, value); }


/*!
 * \ingroup LinearSytem
 * \brief Get a value from a matrix
 */
template<typename T, Concepts::MultiIndex I, Concepts::MultiIndex J> requires(
    MatrixAccessible<T, I, J>)
inline constexpr decltype(auto) get(const T& t, const I& i, const J& j)
{ return Traits::MatrixAccess<T, I, J>::get(t, i, j); }


/*!
 * \ingroup LinearSytem
 * \brief Add to a value of a matrix
 */
template<typename T, Concepts::MultiIndex I, Concepts::MultiIndex J, typename Value> requires(
    Concepts::Addable<MatrixReferenceType<T, I, J>, Value>)
inline constexpr void add(T& t, const I& i, const J& j, const Value& value)
{ set(t, i, j, get(t, i, j) + value); }

} // namespace Dumux::LinearSystem

#endif
