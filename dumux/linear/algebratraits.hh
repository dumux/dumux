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
 * \brief Define traits for linear algebra.
 */
#ifndef DUMUX_LINEAR_ALGEBRA_TRAITS_HH
#define DUMUX_LINEAR_ALGEBRA_TRAITS_HH

#include <utility>
#include <type_traits>

#include <dune/istl/bvector.hh>
#include <dune/common/std/type_traits.hh>

#include<dumux/common/typetraits/vector.hh>

namespace Dumux {

namespace Detail {

template<class T, std::enable_if_t<Dune::IsNumber<std::decay_t<T>>::value, int> = 0>
constexpr std::size_t blockSize() { return 1; }

template<class T, std::enable_if_t<!Dune::IsNumber<std::decay_t<T>>::value, int> = 0>
constexpr std::size_t blockSize() { return std::decay_t<T>::size(); }

} // end namespace Detail

template<class M, class V>
struct LinearAlgebraTraits
{
    using Matrix = M;
    using Vector = V;
};

template<class Assembler, bool isMultiTypeBlockVector>
struct LATraitsFromAssemblerImpl
{
private:
    using VectorPossiblyWithState = typename Assembler::ResidualType;
    using Scalar = std::decay_t<decltype(std::declval<VectorPossiblyWithState>()[0][0])>;
    static constexpr auto blockSize = Detail::blockSize<decltype(std::declval<VectorPossiblyWithState>()[0])>();
    using BlockType = Dune::FieldVector<Scalar, blockSize>;
public:
    using Vector = Dune::BlockVector<BlockType>;
    using Matrix = typename Assembler::JacobianMatrix;
};

template<class Assembler>
struct LATraitsFromAssemblerImpl<Assembler, true>
{
    using Vector = typename Assembler::ResidualType;
    using Matrix = typename Assembler::JacobianMatrix;
};

template<class Assembler>
using LinearAlgebraTraitsFromAssembler = LATraitsFromAssemblerImpl<Assembler,
                                                                   isMultiTypeBlockVector<typename Assembler::ResidualType>::value>;

} // end namespace Dumux

#endif
