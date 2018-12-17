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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Geometry
 * \brief An interface for a set of geometric entities
 * \note This can be used e.g. to contruct a bounding box volume hierarchy of a grid
 * It defines the minimum requirement for such a set
 */
#ifndef DUMUX_GRIDVIEW_GEOMETRIC_ENTITY_SET_HH
#define DUMUX_GRIDVIEW_GEOMETRIC_ENTITY_SET_HH

#include <memory>
#include <dune/grid/common/mcmgmapper.hh>
#include <dumux/common/entitymap.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief An interface for a set of geometric entities
 * \note This can be used e.g. to contruct a bounding box volume hierarchy of a grid
 * It defines the minimum requirement for such a set
 */
template <class GridView, int codim = 0, class Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>>
class GridViewGeometricEntitySet
{
public:
    using Entity = typename GridView::template Codim<codim>::Entity;

    GridViewGeometricEntitySet(const GridView& gridView)
    : gridView_(gridView)
    , mapper_(gridView, Dune::mcmgLayout(Dune::Codim<codim>()))
    , entityMap_(std::make_shared<EntityMap<GridView, codim>>(gridView.grid(), mapper_))
    {}

    GridViewGeometricEntitySet(const GridView& gridView,
                               const Mapper& mapper,
                               std::shared_ptr<const EntityMap<GridView, codim>> entityMap)
    : gridView_(gridView)
    , mapper_(mapper)
    , entityMap_(entityMap)
    {}

    /*!
     * \brief The world dimension of the entity set
     */
    enum { dimensionworld = GridView::dimensionworld };

    /*!
     * \brief the coordinate type
     */
    using ctype = typename GridView::ctype;

    /*!
     * \brief the number of entities in this set
     */
    decltype(auto) size() const
    { return gridView_.size(codim); }

    /*!
     * \brief begin iterator to enable range-based for iteration
     */
    decltype(auto) begin() const
    { return entities(gridView_, Dune::Codim<codim>()).begin(); }

    /*!
     * \brief end iterator to enable range-based for iteration
     */
    decltype(auto) end() const
    { return entities(gridView_, Dune::Codim<codim>()).end(); }

    /*!
     * \brief get an entities index
     */
    std::size_t index(const Entity& e) const
    { return mapper_.index(e); }

    /*!
     * \brief get an entity from an index
     */
    Entity entity(std::size_t index) const
    { return (*entityMap_)[index]; }

private:
    GridView gridView_;
    Mapper mapper_;
    std::shared_ptr<const EntityMap<GridView, codim>> entityMap_;

};

} // end namespace Dumux

#endif
