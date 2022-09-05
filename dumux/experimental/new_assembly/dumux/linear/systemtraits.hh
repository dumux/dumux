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
 * \brief Traits related to operations on matrices/vectors.
 */
#ifndef DUMUX_LINEAR_SYSTEM_TRAITS_HH
#define DUMUX_LINEAR_SYSTEM_TRAITS_HH

#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>

namespace Dumux::LinearSystem::Traits {

//! Extract the scalar type a vector/matrix operates on
template<typename T>
struct Scalar;

//! Function to fill a vector/matrix with a scalar.
//! Specializations are expected to provide the static `fill(T&, const V&)` function
template<typename T, Concepts::Arithmetic V>
struct Fill;

//! Function to add a vector of type `X`, scaled by a value of type `A`, to a vector of type `Y`
//! Specializations are expected to provide the static `axpy(Y&, const A&, const X&)` function
template<typename Y, typename X, Concepts::Arithmetic A>
struct Axpy;

//! Function to access or modify values in a vector at the position given by the multi index `I`
//! Specializations are expected to provide the static `get(T&, const I&)` and the static
//! `set(T&, const I&, const A&)` functions, where `A` is the type to be set at the position
template<typename T, Concepts::MultiIndex I>
struct VectorAccess;

//! Function to access or modify values in a matrix at the position given by the multi indices `I` and `J`
//! Specializations are expected to provide the static `get(T&, const I&, const J&)` and the static
//! `set(T&, const I&, const J&, const A&)` functions, where `A` is the type to be set at the position
template<typename T, Concepts::MultiIndex I, Concepts::MultiIndex J>
struct MatrixAccess;


// default implementations for arithmetic types
template<Concepts::Arithmetic T>
struct Scalar<T> : public std::type_identity<T> {};

template<Concepts::Arithmetic T, Concepts::Arithmetic V>
struct Fill<T, V>
{
    static constexpr void fill(T& t, const V& v)
    { t = v; }
};

template<Concepts::Arithmetic Y, Concepts::Arithmetic X, Concepts::Arithmetic A>
struct Axpy<Y, X, A>
{
    static constexpr void fill(Y& y, const A& a, const X& x)
    { y += a*x; }
};

} // namespace Dumux::LinearSystem::Traits

// include predefined traits
#include <dumux/experimental/new_assembly/dumux/linear/dune/systemtraits.hh>
#include <dumux/experimental/new_assembly/dumux/linear/eigen/systemtraits.hh>

#endif
