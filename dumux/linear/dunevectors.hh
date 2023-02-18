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
 * \brief Helper to extract native Dune vector types from particular Dumux types
 */
#ifndef DUMUX_LINEAR_DUNE_VECTORS_HH
#define DUMUX_LINEAR_DUNE_VECTORS_HH

#include <type_traits>
#include <dune/common/typetraits.hh>
#include <dune/common/std/type_traits.hh>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

// implementation details specific to this header
namespace Dumux::Detail::DuneVectors {

template<class T, std::enable_if_t<Dune::IsNumber<std::decay_t<T>>::value, int> = 0>
constexpr std::size_t blockSize() { return 1; }

template<class T, std::enable_if_t<!Dune::IsNumber<std::decay_t<T>>::value, int> = 0>
constexpr std::size_t blockSize() { return std::decay_t<T>::size(); }

template<class C> using StateDetector = decltype(std::declval<C>()[0].state());

} // end namespace Dumux::Detail::DuneVectors

// implementation details used elsewhere in Dumux
namespace Dumux::Detail {

template<class V, bool hasState>
struct NativeDuneVectorTypeImpl
{ using type = V; };

// only support single-nested block vectors for now
template<class V>
struct NativeDuneVectorTypeImpl<V, true>
{
    using Scalar = std::decay_t<decltype(std::declval<V>()[0][0])>;
    static constexpr auto blockSize = Detail::DuneVectors::blockSize<decltype(std::declval<V>()[0])>();
    using BlockType = Dune::FieldVector<Scalar, blockSize>;
    using type = Dune::BlockVector<BlockType>;
};

// get native dune type from dumux vector (possibly with state)
template<class V>
struct NativeDuneVectorType
{
    using type = typename NativeDuneVectorTypeImpl<
        V, Dune::Std::is_detected<Detail::DuneVectors::StateDetector, V>{}
    >::type;
};

// specialization for multitype vector
template<class... Args>
struct NativeDuneVectorType<Dune::MultiTypeBlockVector<Args...>>
{
    using type = Dune::MultiTypeBlockVector<
        typename NativeDuneVectorType<Args>::type...
    >;
};

} // end namespace Dumux::Detail

#endif
