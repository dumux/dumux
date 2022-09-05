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
 * \ingroup Common
 * \brief Default type for storing data.
 */
#ifndef DUMUX_COMMON_STORAGAE_HH
#define DUMUX_COMMON_STORAGAE_HH

#include <vector>
#include <type_traits>

#include <dune/common/reservedvector.hh>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>

namespace Dumux {

template<typename T, std::size_t s>
using ReservedStorage = Dune::ReservedVector<T, s>;


#ifndef DOXYGEN
namespace Detail {

template<typename T, Concepts::Size Size, Concepts::Size auto size>
struct DefaultStorage;

template<typename T, std::integral Size, std::integral auto size>
struct DefaultStorage<T, Size, size>
: public std::type_identity<ReservedStorage<T, size>>
{};

template<typename T, DynamicSize size>
struct DefaultStorage<T, DynamicSize, size>
: public std::type_identity<std::vector<T>>
{};

} // namespace Detail
#endif // DOXYGEN

template<typename T, Concepts::Size auto size>
using DefaultStorage = Detail::DefaultStorage<T, std::decay_t<decltype(size)>, size>::type;

} // namespace Dumux

#endif
