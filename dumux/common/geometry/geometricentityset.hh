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
#ifndef DUMUX_GEOMETRIC_ENTITY_SET_HH
#define DUMUX_GEOMETRIC_ENTITY_SET_HH

#include <memory>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dumux/common/entitymap.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief An interface for a set of geometric entities based on a GridView
 * \note This can be used e.g. to contruct a bounding box volume hierarchy of a grid
 * It defines the minimum requirement for such a set
 */
template <class GridView, int codim = 0, class Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>>
class GridViewGeometricEntitySet
{
    using EntityMap = Dumux::EntityMap<GridView, codim>;
public:
    using Entity = typename GridView::template Codim<codim>::Entity;

    GridViewGeometricEntitySet(const GridView& gridView)
    : GridViewGeometricEntitySet(gridView, Mapper(gridView, Dune::mcmgLayout(Dune::Codim<codim>())))
    {}

    GridViewGeometricEntitySet(const GridView& gridView, const Mapper& mapper)
    : gridView_(gridView)
    , mapper_(mapper)
    , entityMap_(std::make_shared<EntityMap>(gridView.grid(), mapper_))
    {}

    GridViewGeometricEntitySet(const GridView& gridView,
                               const Mapper& mapper,
                               std::shared_ptr<const EntityMap> entityMap)
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
    std::shared_ptr<const EntityMap> entityMap_;

};

/*!
 * \ingroup Geometry
 * \brief An interface for a set of geometric entities
 * \note This can be used e.g. to contruct a bounding box volume hierarchy of a grid
 * It defines the minimum requirement for such a set
 */
template<class GeoType>
class GeometriesEntitySet
{
    /*!
     * \brief Helper class to wrap the Dune entity
     */
    class EntityWrapper
    {
    public:
        using GeometryType = GeoType;

        /*!
         * \brief Copy constructor
         */
        EntityWrapper(const GeometryType& geo, const std::size_t seed) : geo_(geo), seed_(seed) {}

        /*!
         * \brief Move constructor
         */
        EntityWrapper(GeometryType&& geo, const std::size_t seed) : geo_(std::move(geo)), seed_(seed) {}

        /*!
         * \brief Returns the geometry
         */
        const GeometryType& geometry() const
        { return geo_; }

        /*!
         * \brief Returns the seed (index) of the geometry
         */
        std::size_t seed() const
        { return seed_; }

    private:
        GeometryType geo_;
        std::size_t seed_;
    };

    using GeometryType = typename EntityWrapper::GeometryType;

public:
    using Entity = EntityWrapper;

    /*!
     * \brief Copy constructor for single geometry
     */
    GeometriesEntitySet(const GeometryType& geometry)
    {
        const std::size_t seed = 0;
        entities_.emplace_back(geometry, seed);
    }

    /*!
     * \brief Move constructor for single geometry
     */
    GeometriesEntitySet(GeometryType&& geometry)
    {
        const std::size_t seed = 0;
        entities_.emplace_back(std::move(geometry), seed);
    }

    /*!
     * \brief Constructor for initializer_list
     */
    GeometriesEntitySet(std::initializer_list<GeometryType>&& geometries)
    {
        std::size_t seed = 0;
        // note: std::initializer_list::begin() returns const T*,
        // thus no moving will be performed and only the copy ctor of
        // EntityWrapper can be called
        for (auto&& g : geometries)
            entities_.emplace_back(g, seed++);
    }

    /*!
     * \brief Copy constructor for a vector of geometries
     */
    GeometriesEntitySet(const std::vector<GeometryType>& geometries)
    {
        std::size_t seed = 0;
        for (auto&& g : geometries)
            entities_.emplace_back(g, seed++);
    }

    /*!
     * \brief Move constructor for a vector of geometries
     */
    GeometriesEntitySet(std::vector<GeometryType>&& geometries)
    {
        std::size_t seed = 0;
        for (auto&& g : geometries)
            entities_.emplace_back(std::move(g), seed++);
    }

    /*!
     * \brief The world dimension of the entity set
     */
    enum { dimensionworld = GeometryType::coorddimension };

    /*!
     * \brief the coordinate type
     */
    using ctype = typename GeometryType::ctype;

    /*!
     * \brief the number of entities in this set
     */
    decltype(auto) size() const
    { return entities_.size(); }

    /*!
     * \brief begin iterator to enable range-based for iteration
     */
    decltype(auto) begin() const
    { return entities_.begin(); }

    /*!
     * \brief end iterator to enable range-based for iteration
     */
    decltype(auto) end() const
    { return entities_.end(); }

    /*!
     * \brief get an entities index
     */
    template<class Entity>
    std::size_t index(const Entity& e) const
    { return e.seed(); }

    /*!
     * \brief get an entity from an index
     */
    Entity entity(std::size_t index) const
    { return entities_[index]; }

private:
    std::vector<Entity> entities_;
};

} // end namespace Dumux

#endif
