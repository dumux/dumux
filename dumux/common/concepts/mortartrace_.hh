// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concepts for the trace of a subdomain on a mortar domain
 */
#ifndef DUMUX_CONCEPTS_MORTAR_TRACE__HH
#define DUMUX_CONCEPTS_MORTAR_TRACE__HH

#include <concepts>
#include <cstddef>
#include <memory>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/common/concepts/entityset_.hh>

namespace Dumux::Concept {

/*!
 * \ingroup Concepts
 * \brief The trace of a subdomain on a mortar domain: a set of trace cells with geometry,
 *        indexed contiguously from zero as given by its entity set.
 *
 * Data exchanged between a subdomain and a mortar lives on the trace cells, and the
 * projections between the two are assembled over the intersections of the trace cells with
 * the mortar cells. What a trace cell is depends on the subdomain: a boundary facet of a
 * bulk grid, or the cross-section of a network's boundary entity with the interface.
 */
template<class T>
concept MortarTrace = requires(const T& trace)
{
    typename T::EntitySet;
    requires GeometricEntitySet<typename T::EntitySet>;
    { trace.size() } -> std::convertible_to<std::size_t>;
    { trace.entitySet() } -> std::convertible_to<std::shared_ptr<const typename T::EntitySet>>;
    { trace.boundingBoxTree() } -> std::convertible_to<const BoundingBoxTree<typename T::EntitySet>&>;
};

/*!
 * \ingroup Concepts
 * \brief A mortar trace whose cells form a grid, so that conforming function spaces can be
 *        defined on it.
 */
template<class T>
concept GridMortarTrace = MortarTrace<T> && requires(const T& trace)
{
    typename T::GridView;
    { trace.gridView() } -> std::convertible_to<typename T::GridView>;
    requires std::same_as<
        typename T::EntitySet::Entity,
        typename T::GridView::template Codim<0>::Entity
    >;
};

} // end namespace Dumux::Concept

#endif
