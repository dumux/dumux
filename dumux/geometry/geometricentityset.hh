// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief An interface for a set of geometric entities
 * \note This can be used e.g. to construct a bounding box volume hierarchy of a grid
 * It defines the minimum requirement for such a set
 */
#ifndef DUMUX_GEOMETRY_GEOMETRIC_ENTITY_SET_HH
#define DUMUX_GEOMETRY_GEOMETRIC_ENTITY_SET_HH

#include <memory>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dumux/common/entitymap.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief An interface for a set of geometric entities based on a GridView
 * \note This can be used e.g. to construct a bounding box volume hierarchy of a grid
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
    { assert(index < entityMap_->size()); return (*entityMap_)[index]; }

private:
    GridView gridView_;
    Mapper mapper_;
    std::shared_ptr<const EntityMap> entityMap_;
};

/*!
 * \ingroup Geometry
 * \brief An interface for a set of geometric entities
 * \note This can be used e.g. to construct a bounding box volume hierarchy of a grid
 * It defines the minimum requirement for such a set
 */
template<class GeoType>
class GeometriesEntitySet
{
    /*!
     * \brief Wrapper to turn a geometry into a geometric entity
     */
    class EntityWrapper
    {
    public:
        using Geometry = GeoType;

        /*!
         * \brief Constructor
         */
        EntityWrapper(const Geometry& geo, const std::size_t index) : geo_(geo), index_(index) {}

        /*!
         * \brief Constructor
         */
        EntityWrapper(Geometry&& geo, const std::size_t index) : geo_(std::move(geo)), index_(index) {}

        /*!
         * \brief Returns the geometry
         */
        const Geometry& geometry() const
        { return geo_; }

        /*!
         * \brief Returns the index of the geometry
         */
        std::size_t index() const
        { return index_; }

    private:
        Geometry geo_;
        std::size_t index_;
    };

public:
    using Entity = EntityWrapper;

    /*!
     * \brief Constructor for initializer_list
     */
    GeometriesEntitySet(std::initializer_list<typename Entity::Geometry>&& geometries)
    {
        std::size_t index = 0;
        // note: std::initializer_list::begin() returns const T*,
        // thus no moving will be performed and only the copying ctor of
        // EntityWrapper can be called
        for (auto&& g : geometries)
            entities_.emplace_back(g, index++);
    }

    /*!
     * \brief Constructor for a vector of geometries
     */
    GeometriesEntitySet(const std::vector<typename Entity::Geometry>& geometries)
    {
        std::size_t index = 0;
        for (auto&& g : geometries)
            entities_.emplace_back(g, index++);
    }

    /*!
     * \brief Constructor for a vector of geometries
     */
    GeometriesEntitySet(std::vector<typename Entity::Geometry>&& geometries)
    {
        std::size_t index = 0;
        for (auto&& g : geometries)
            entities_.emplace_back(std::move(g), index++);
    }

    /*!
     * \brief The world dimension of the entity set
     */
    enum { dimensionworld = Entity::Geometry::coorddimension };

    /*!
     * \brief the coordinate type
     */
    using ctype = typename Entity::Geometry::ctype;

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
    { return e.index(); }

    /*!
     * \brief get an entity from an index
     */
    const Entity& entity(std::size_t index) const
    { assert(index < entities_.size()); return entities_[index]; }

private:
    std::vector<Entity> entities_;
};

} // end namespace Dumux

#endif
