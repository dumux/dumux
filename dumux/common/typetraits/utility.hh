// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup TypeTraits
 * \brief Utilities for template meta programming
 */
#ifndef DUMUX_COMMON_TYPETRAITS_UTILITY_HH
#define DUMUX_COMMON_TYPETRAITS_UTILITY_HH

namespace Dumux {

/*
 * \ingroup TypeTraits
 * \brief create a variadic template from indexed types
 * \tparam V a variadic template that we want to create
 * \tparam T an indexed type (type that gets an index as template parameter)
 * \tparam U the list of indices
 */
template <template<typename... Args> class Variadic, template<std::size_t> class Indexed, class U>
struct makeFromIndexedType;

template <template<typename... Args> class Variadic, template<std::size_t> class Indexed, std::size_t... IndexSeq>
struct makeFromIndexedType<Variadic, Indexed, std::index_sequence<IndexSeq...>>
{
    using type = Variadic<Indexed<IndexSeq>...>;
};

} // end namespace Dumux

#endif
