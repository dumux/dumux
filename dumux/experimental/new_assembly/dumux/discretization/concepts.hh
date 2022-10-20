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
 * \ingroup Discretization
 * \brief Basic discretization-related concepts
 */
#ifndef DUMUX_DISCRETIZATION_CONCEPTS_HH
#define DUMUX_DISCRETIZATION_CONCEPTS_HH

#include <concepts>
#include <utility>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>

namespace Dumux::Concepts {

template<typename T>
concept GridView = requires(const T& t) {
    { int{T::dimension} };
    { int{T::dimensionworld} };

    typename T::Intersection;
    typename T::template Codim<0>;
    typename T::template Codim<0>::Entity;
    typename T::template Codim<T::dimension>;
    typename T::template Codim<T::dimension>::Entity;

    { elements(t) } -> RangeOf<const typename T::template Codim<0>::Entity>;
    { vertices(t) } -> RangeOf<const typename T::template Codim<T::dimension>::Entity>;
    { intersections(t, std::declval<typename T::template Codim<0>::Entity>()) } -> RangeOf<const typename T::Intersection>;
};

template<typename T>
concept NetworkGridView = GridView<T> and (int{T::dimension} < int{T::dimensionworld});

} // end namespace Dumux::Concepts

#endif
