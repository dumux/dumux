// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief A map from indices to entities using grid entity seeds
 */
#ifndef DUMUX_ENTITY_INDEX_MAP_HH
#define DUMUX_ENTITY_INDEX_MAP_HH

#include <vector>
#include <utility>
#include <dune/geometry/dimension.hh>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief A map from indices to entities using grid entity seeds
 */
template <class GridView, int codim = 0>
class EntityMap
{
public:
    using Grid = typename GridView::Traits::Grid;
    using Entity = typename Grid::template Codim<codim>::Entity;
    using EntitySeed = typename Grid::template Codim<codim>::EntitySeed;

    //! constructor moving a ready seed list in here
    EntityMap(const Grid& grid, std::vector<EntitySeed>&& seeds)
    : grid_(grid)
    , seeds_(std::move(seeds))
    {}

    //! constructor with all entities of codim
    template<class Mapper>
    EntityMap(const Grid& grid, const Mapper& mapper)
    : grid_(grid)
    {
        update(mapper);
    }

    //! update the map after the grid changed
    void update(std::vector<EntitySeed>&& seeds)
    { seeds_.swap(std::move(seeds)); }

    //! update the map after the grid changed
    template<class Mapper>
    void update(const Mapper& mapper)
    {
        const auto& gv = grid_.leafGridView();
        seeds_.resize(gv.size(codim));
        for (const auto& entity : entities(gv, Dune::Codim<codim>()))
            seeds_[mapper.index(entity)] = entity.seed();
    }

    //! get an element from an index i
    Entity operator[](std::size_t i) const
    { return grid_.entity(seeds_[i]); }

    //! get the size of the map
    std::size_t size() const
    { return seeds_.size(); }

private:
    const Grid& grid_;
    std::vector<EntitySeed> seeds_;
};

template<class GridView>
using ElementMap = EntityMap<GridView, 0>;

} // end namespace Dumux

#endif
