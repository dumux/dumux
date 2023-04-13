// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
