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
 * \brief Discretization-related concepts.
 */
#ifndef DUMUX_DISCRETIZATION_CONCEPTS_HH
#define DUMUX_DISCRETIZATION_CONCEPTS_HH

#include <type_traits>
#include <concepts>
#include <ranges>

#include <dumux/discretization/method.hh>

namespace Dumux::Concepts {

template<typename T>
concept Viewable = requires(T& t) {
    typename T::LocalView;
    { localView(t) } -> std::convertible_to<typename T::LocalView>;
};

template<typename T>
concept GridGeometry = Viewable<std::add_const_t<T>> and requires {
    T::discMethod;
};

template<typename T>
concept GridGeometryLocalView = requires (T& t) {
    typename T::GridGeometry;
    GridGeometry<typename T::GridGeometry>;

    { t.bind(std::declval<const typename T::GridGeometry::GridView::template Codim<0>::Entity&>() ) };
};

template<typename T>
concept FVGridGeometryLocalView = GridGeometryLocalView<T> and requires (const T& t) {
    typename T::SubControlVolume;
    typename T::SubControlVolumeFace;
    { t.numScv() } -> std::integral;
    { t.numScvf() } -> std::integral;
    { scvs(t) } -> std::ranges::range;
    { scvfs(t) } -> std::ranges::range;
    std::same_as<std::ranges::range_value_t<decltype(scvs(t))>, typename T::SubControlVolume>;
    std::same_as<std::ranges::range_value_t<decltype(scvfs(t))>, typename T::SubControlVolumeFace>;
} and requires (const T& t,
                const typename T::SubControlVolume& scv,
                const typename T::SubControlVolumeFace& scvf) {
    { t.dofIndex(scv) } -> std::integral;
};

template<typename T>
concept FVGridGeometry = GridGeometry<T> and requires (const T& t) {
    typename T::SubControlVolume;
    typename T::SubControlVolumeFace;
    FVGridGeometryLocalView<typename T::LocalView>;
};

template<typename T>
concept CCTpfaGridGeometry = FVGridGeometry<T> and T::discMethod == DiscretizationMethods::cctpfa;

template<typename T>
concept CCTpfaGridGeometryLocalView = FVGridGeometryLocalView<T> and CCTpfaGridGeometry<typename T::GridGeometry>;

} // namespace Dumux::Concepts

#endif
