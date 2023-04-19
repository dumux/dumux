// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits to be used with matrix types
 */
#ifndef DUMUX_TYPETRAITS_MATRIX_HH
#define DUMUX_TYPETRAITS_MATRIX_HH

#include <type_traits>

// Forward declare to avoid includes
namespace Dune {
template <class Block, class Allocator>
class BCRSMatrix;

template <class FirstRow, class ... Args>
class MultiTypeBlockMatrix;
} // end namespace Dune

namespace Dumux {

//! Helper type to determine whether a given type is a Dune::BCRSMatrix
template<class T>
struct isBCRSMatrix : public std::false_type {};

template<class B, class A>
struct isBCRSMatrix<Dune::BCRSMatrix<B, A>> : public std::true_type {};

//! Helper type to determine whether a given type is a Dune::MultiTypeBlockMatrix
template<class... Args>
struct isMultiTypeBlockMatrix : public std::false_type {};

template<class... Args>
struct isMultiTypeBlockMatrix<Dune::MultiTypeBlockMatrix<Args...>> : public std::true_type {};

} // end namespace Dumux

#endif
