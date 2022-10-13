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

template<typename T, typename GridView>
concept ElementMapper = requires(const T& t, const typename GridView::template Codim<0>::Entity& e) {
    { t.index(e) } -> std::integral;
};

template<typename T, typename GridView>
concept UpdatableElementMapper = ElementMapper<T, GridView> and requires(T& t, const GridView& gv) {
    { t.update(gv) };
};

template<typename T>
concept CCGridGeometryLocalView = requires(T& localView, const T& constLocalView) {
    typename T::GridGeometry;
    typename T::Element;
    typename T::SubControlVolume;
    typename T::SubControlVolumeFace;

    { constLocalView.gridGeometry() } -> std::same_as<const typename T::GridGeometry&>;
    { constLocalView.element() } -> std::same_as<const typename T::Element&>;

    { localView.bind(std::declval<const typename T::Element&>()) };
    { localView.bindElement(std::declval<const typename T::Element&>()) };

    { scvs(constLocalView) } -> RangeOf<typename T::SubControlVolume>;
    { scvfs(constLocalView) } -> RangeOf<typename T::SubControlVolumeFace>;

    { neighborScvs(constLocalView) } -> RangeOf<typename T::SubControlVolume>;
    { neighborScvfs(constLocalView) } -> RangeOf<typename T::SubControlVolumeFace>;

    { constLocalView.onBoundary(std::declval<const typename T::SubControlVolumeFace&>()) } -> std::convertible_to<bool>;
    { constLocalView.insideScv(std::declval<const typename T::SubControlVolumeFace&>()) } -> std::convertible_to<typename T::SubControlVolume>;
    { constLocalView.outsideScv(std::declval<const typename T::SubControlVolumeFace&>()) } -> std::convertible_to<typename T::SubControlVolume>;
    { constLocalView.outsideScv(std::declval<const typename T::SubControlVolumeFace&>(), unsigned{}) } -> std::convertible_to<typename T::SubControlVolume>;
    { constLocalView.numOutsideNeighbors(std::declval<const typename T::SubControlVolumeFace&>()) } -> std::integral;

    { constLocalView.geometry(std::declval<const typename T::SubControlVolume&>()) };
    { constLocalView.geometry(std::declval<const typename T::SubControlVolumeFace&>()) };
};

} // end namespace Dumux::Concepts

#endif
