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
 * \copydoc Dumux::DefaultStorage
 */
#ifndef DUMUX_COMMON_STORAGE_HH
#define DUMUX_COMMON_STORAGE_HH

#include <vector>
#include <ranges>
#include <cassert>
#include <type_traits>

#include <dune/common/reservedvector.hh>

#include <dumux/experimental/new_assembly/dumux/common/size.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<typename T, std::integral auto size>
using ReservedStorage = Dune::ReservedVector<T, size>;

template<typename T>
struct IsReserved : public std::false_type {};

template<typename T, std::integral auto size>
struct IsReserved<ReservedStorage<T, size>> : public std::true_type {};

template<typename T, typename Size, auto size>
struct DefaultStorage;

template<typename T, std::integral Size, std::integral auto size>
struct DefaultStorage<T, Size, size>
: public std::type_identity<Dune::ReservedVector<T, size>>
{};

template<typename T, Size::Dynamic size>
struct DefaultStorage<T, Size::Dynamic, size>
: public std::type_identity<std::vector<T>>
{};

} // namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup Common
 * \brief Default container type for storing data dynamically.
 *        Depending on the given template argument for the (maximum) container
 *        size (statically known or fully dynamic size), the underlying container
 *        may use statically pre-allocated memory, while still exposing the
 *        interface of a dynamic (std::)vector. This avoids dynamic memory allocation,
 *        but can only be used if an upper bound on the required memory is know.
 * \todo We may want to use std::pmr::vector with a statically-sized
 *       memory_ressource such that we don't run into segfaults when using
 *       more space than was reserved with `s` (as happens with Dune::ReservedVector).
 * \note We inherit privately from the underlying vector type such that we can
 *       completely control the exposed public interface with using directives.
 */
template<typename T, Concepts::Size auto s>
class DefaultStorage
: private Detail::DefaultStorage<T, std::decay_t<decltype(s)>, s>::type
{
    using ParentType = Detail::DefaultStorage<T, std::decay_t<decltype(s)>, s>::type;
    static constexpr bool isReserved = Detail::IsReserved<ParentType>::value;

public:
    using ParentType::ParentType;
    using typename ParentType::iterator;
    using typename ParentType::const_iterator;

    using ParentType::size;
    using ParentType::empty;

    using ParentType::clear;
    using ParentType::push_back;
    using ParentType::emplace_back;

    using ParentType::operator[];
    using ParentType::begin;
    using ParentType::end;

    using ParentType::resize;
    void resize(std::integral auto newSize, const T& value)
    {
        if constexpr (isReserved)
        {
            ParentType::resize(newSize);
            std::ranges::fill(*this, value);
        }
        else
            ParentType::resize(newSize, value);
    }

    void reserve(std::integral auto newSize)
    {
        if constexpr (isReserved)
            assert(newSize <= this->capacity());
        else
            ParentType::reserve(newSize);
    }

    void assign(std::integral auto newSize, const T& value)
    {
        if constexpr (isReserved)
        {
            this->resize(newSize);
            std::ranges::fill(*this, value);
        }
        else
            this->assign(newSize, value);
    }

    void shrink_to_fit()
    {
        if constexpr (!isReserved)
            ParentType::shrink_to_fit();
    }
};

} // namespace Dumux

#endif
