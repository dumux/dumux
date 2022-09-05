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
 * \brief Traits related to Dune matrix/vector types.
 */
#ifndef DUMUX_LINEAR_DUNE_SYSTEM_TRAITS_HH
#define DUMUX_LINEAR_DUNE_SYSTEM_TRAITS_HH

#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>
#include <dumux/experimental/new_assembly/dumux/linear/systemtraits.hh>

#include "detail.hh"
#include "access.hh"

namespace Dumux {
namespace LinearSystem::Traits {

#ifndef DOXYGEN

template<typename T, Concepts::Arithmetic V>
struct DuneFill
{
    static void fill(T& t, const V& v)
    { t = v; }
};

template<typename Y, typename X, Concepts::Arithmetic A>
struct DuneAxpy
{
    static void axpy(Y& y, const A& a, const X& x)
    { y.axpy(a, x); }
};

#endif // DOXYGEN

template<typename T> requires(Detail::DuneVector<T> or Detail::DuneMatrix<T>)
struct Scalar<T>
: public std::type_identity<typename T::field_type>
{};

template<typename T, Concepts::Arithmetic V> requires(
    Detail::DuneVector<T> or Detail::DuneMatrix<T>)
struct Fill<T, V> : public DuneFill<T, V> {};

template<typename Y, typename X, Concepts::Arithmetic A> requires(
    Detail::DuneVector<Y> or Detail::DuneMatrix<Y>)
struct Axpy<Y, X, A> : public DuneAxpy<Y, X, A> {};

template<typename T, Concepts::MultiIndex I> requires(Detail::DuneVector<T>)
struct VectorAccess<T, I>
{
    static decltype(auto) get(const T& vector, const I& multiIndex)
    { return DuneBlockVectorAccess::get(vector, multiIndex); }

    template<typename Value>
    static void set(T& vector, const I& multiIndex, const Value& v)
    { DuneBlockVectorAccess::get(vector, multiIndex) = v; }
};

template<typename T, Concepts::MultiIndex I, Concepts::MultiIndex J> requires(Detail::DuneMatrix<T>)
struct MatrixAccess<T, I, J>
{
    static decltype(auto) get(const T& matrix, const I& i, const J& j)
    { return DuneBlockMatrixAccess::get(matrix, i, j); }

    template<typename Value>
    static void set(T& matrix, const I& i, const J& j, const Value& v)
    { DuneBlockMatrixAccess::get(matrix, i, j) = v; }
};

} // namespace LinearSystem::Traits
} // namespace Dumux

#endif
