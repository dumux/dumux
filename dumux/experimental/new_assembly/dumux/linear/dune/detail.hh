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
#ifndef DUMUX_LINEAR_DUNE_DETAIL_HH
#define DUMUX_LINEAR_DUNE_DETAIL_HH
#ifndef DOXYGEN

#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>

// forward declarations
namespace Dune {

template<typename T, int size>
class FieldVector;

template<typename T, typename Allocator>
class DynamicVector;

template<typename BlockType, typename Allocator>
class BlockVector;

template<typename... Blocks>
class MultiTypeBlockVector;

template<typename T, int rows, int cols>
class FieldMatrix;

template<typename T>
class DynamicMatrix;

template<typename BlockType, typename Allocator>
class BCRSMatrix;

template<typename FirstRow, typename... Rows>
class MultiTypeBlockMatrix;

} // namespace Dune

namespace Dumux::LinearSystem::Detail {

template<typename T>
struct IsDuneVector : public std::false_type {};

template<typename T, int size>
struct IsDuneVector<Dune::FieldVector<T, size>> : public std::true_type {};

template<typename T, typename A>
struct IsDuneVector<Dune::DynamicVector<T, A>> : public std::true_type {};

template<typename BlockType, typename Allocator>
struct IsDuneVector<Dune::BlockVector<BlockType, Allocator>> : public std::true_type {};

template<typename... Blocks>
struct IsDuneVector<Dune::MultiTypeBlockVector<Blocks...>> : public std::true_type {};


template<typename T>
struct IsDuneMatrix : public std::false_type {};

template<typename T, int rows, int cols>
struct IsDuneMatrix<Dune::FieldMatrix<T, rows, cols>> : public std::true_type {};

template<typename T>
struct IsDuneMatrix<Dune::DynamicMatrix<T>> : public std::true_type {};

template<typename BlockType, typename Allocator>
struct IsDuneMatrix<Dune::BCRSMatrix<BlockType, Allocator>> : public std::true_type {};

template<typename FirstRow, typename... Rows>
struct IsDuneMatrix<Dune::MultiTypeBlockMatrix<FirstRow, Rows...>> : public std::true_type {};


template<typename T>
struct IsDynamic : public std::false_type {};

template<typename T, typename A>
struct IsDynamic<Dune::DynamicVector<T, A>> : public std::true_type {};

template<typename T, typename A>
struct IsDynamic<Dune::BlockVector<T, A>> : public std::true_type {};

template<typename T>
struct IsDynamic<Dune::DynamicMatrix<T>> : public std::true_type {};

template<typename T, typename A>
struct IsDynamic<Dune::BCRSMatrix<T, A>> : public std::true_type {};


template<typename T>
concept DuneVector = IsDuneVector<T>::value;

template<typename T>
concept DuneMatrix = IsDuneMatrix<T>::value;

template<typename T>
concept DuneBlockVector = DuneVector<T> and DuneVector<typename T::value_type>;

template<typename T>
concept DuneBlockMatrix = DuneMatrix<T> and DuneMatrix<typename T::block_type>;

template<typename T>
concept DynamicDuneVector = DuneVector<T> and IsDynamic<T>::value;

template<typename T>
concept DynamicDuneMatrix = DuneMatrix<T> and IsDynamic<T>::value;

} // namespace Dumux::LinearSystem::Detail

#endif // DOXYGEN

#endif
