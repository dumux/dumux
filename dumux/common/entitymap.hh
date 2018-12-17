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
 * \brief A map from indices to entities using grid entity seeds
 */
#ifndef DUMUX_ENTITY_INDEX_MAP_HH
#define DUMUX_ENTITY_INDEX_MAP_HH

#include <vector>
#include <utility>
#include <dune/geometry/dimension.hh>

namespace Dumux {

/*!
 * \ingroup Common
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

    //! constructor with all entites of codim
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
