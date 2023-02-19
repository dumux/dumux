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
 * \brief Define traits for linear algebra backends
 */
#ifndef DUMUX_LINEAR_ALGEBRA_TRAITS_HH
#define DUMUX_LINEAR_ALGEBRA_TRAITS_HH

#include <utility>
#include <type_traits>

#include <dumux/common/typetraits/vector.hh>
#include <dumux/linear/matrixconverter.hh>

namespace Dumux::Detail::LATraits {

template<class Assembler, bool isMultiType = false>
struct LATraitsFromAssemblerImpl
{
    using Vector = typename Assembler::ResidualType;
    using Matrix = typename Assembler::JacobianMatrix;
    using SingleTypeVector = Vector;
    using SingleTypeMatrix = Matrix;
};

template<class Assembler>
struct LATraitsFromAssemblerImpl<Assembler, true>
{
    using Vector = typename Assembler::ResidualType;
    using Matrix = typename Assembler::JacobianMatrix;
    using SingleTypeVector = decltype(VectorConverter<Vector>::multiTypeToBlockVector(std::declval<Vector>()));
    using SingleTypeMatrix = decltype(MatrixConverter<Matrix>::multiTypeToBCRSMatrix(std::declval<Matrix>()));
};

} // end namespace Dumux::Detail::LATraits

namespace Dumux {

/*
 * \ingroup Linear
 * \brief Traits providing linear algebra types (vector, matrix)
 */
template<class M, class V, class STM = M, class STV = V>
struct LinearAlgebraTraits
{
    using Matrix = M;
    using Vector = V;
    using SingleTypeMatrix = STM;
    using SingleTypeVector = STV;
};

/*
 * \ingroup Linear
 * \brief Helper to extract linear algebra types from an assembler
 */
template<class Assembler>
using LinearAlgebraTraitsFromAssembler =
    Detail::LATraits::LATraitsFromAssemblerImpl<
        Assembler, isMultiTypeBlockVector<typename Assembler::ResidualType>::value
    >;

} // end namespace Dumux

#endif
