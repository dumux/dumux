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
 * \ingroup TypeTraits
 */
#ifndef DUMUX_TYPE_TRAITS_HH
#define DUMUX_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/common/typetraits.hh>

namespace Dumux {
    /*!
     * \brief Template which always yields a false value
     * \tparam T Some type.
     */
    template<typename T>
    struct AlwaysFalse : public std::false_type {};

    /*! \brief We define our own is_indexable type in order
     *         to avoid several version checks throughout dumux.
     *         This should be deleted when the deprecation phase is over.
     */
    template<typename T, typename I = std::size_t>
    using IsIndexable = typename Dune::IsIndexable<T, I>;

} // end namespace Dumux
#endif
