// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concept for finite sets of entities with geometry
 */
#ifndef DUMUX_CONCEPTS_ENTITYSET__HH
#define DUMUX_CONCEPTS_ENTITYSET__HH

#include <concepts>
#include <cstddef>
#include <iterator>

namespace Dumux::Concept {

/*!
 * \ingroup Concepts
 * \brief A finite, contiguously indexed set of entities that carry a geometry, as consumed
 *        by bounding box trees, intersection entity sets and bases defined per entity.
 *
 * The entities may be those of a grid view or free geometries; the set is what makes
 * them addressable by an index in `[0, size())` and back.
 */
template<class ES>
concept GeometricEntitySet = requires(const ES& set, const typename ES::Entity& entity, std::size_t i)
{
    typename ES::Entity;
    typename ES::Entity::Geometry;
    typename ES::ctype;
    { ES::dimensionworld } -> std::convertible_to<int>;
    { set.size() } -> std::convertible_to<std::size_t>;
    { set.index(entity) } -> std::convertible_to<std::size_t>;
    { set.entity(i) } -> std::convertible_to<typename ES::Entity>;
    { entity.geometry() } -> std::convertible_to<typename ES::Entity::Geometry>;
    { set.begin() } -> std::input_iterator;
    { set.end() } -> std::sentinel_for<decltype(set.begin())>;
};

} // end namespace Dumux::Concept

#endif
