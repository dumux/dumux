// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief A basic implementation of a grid geometry with some common interfaces
 */
#ifndef DUMUX_DISCRETIZATION_BASIC_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_BASIC_GRID_GEOMETRY_HH

#include <memory>
#include <utility>
#include <type_traits>

#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/common/entitymap.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief An implementation of a grid geometry with some basic features
 * \tparam GV the grid view type
 * \tparam EM the type of the element mapper
 * \tparam VM the type of the vertex mapper
 */
template<class GV, class EM, class VM>
class BasicGridGeometry
{
    using ElementMap = EntityMap<GV, 0>;
    using ElementSet = GridViewGeometricEntitySet<GV, 0, EM>;
    using BoundingBoxTree = Dumux::BoundingBoxTree<ElementSet>;

    static constexpr int dim = GV::dimension;
    static constexpr int dimWorld = GV::dimensionworld;

    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using Element = typename GV::template Codim<0>::Entity;

public:
    //! export the grid type
    using Grid = typename GV::Grid;
    //! export the grid view type
    using GridView = GV;
    //! export the global coordinate type
    using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
    //! export the element mapper type
    using ElementMapper = EM;
    //! export the vertex mapper type
    using VertexMapper = VM;

    /*!
     * \ingroup Discretization
     * \brief Constructor computes the bounding box of the entire domain, for e.g. setting boundary conditions
     * \param gridView the grid view on which to construct the grid geometry
     */
    BasicGridGeometry(const GridView& gridView)
    : gridView_(gridView)
    , elementMapper_(makeElementMapper_(gridView))
    , vertexMapper_(makeVertexMapper_(gridView))
    , bBoxMin_(std::numeric_limits<double>::max())
    , bBoxMax_(-std::numeric_limits<double>::max())
    {
        computeGlobalBoundingBox_();
        update_();
    }

    /*!
     * \brief Update internal state after grid changed
     */
    void update(const GridView& gridView)
    {
        gridView_ = gridView;
        update_();
    }

    /*!
     * \brief Update internal state after grid changed
     */
    void update(GridView&& gridView)
    {
        gridView_ = std::move(gridView);
        update_();
    }

    /*!
     * \brief Return the gridView this grid geometry object lives on
     */
    const GridView& gridView() const
    { return gridView_; }

    /*!
     * \brief Returns the mapper for vertices to indices for constant grids.
     */
    const VertexMapper &vertexMapper() const
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices for constant grids.
     */
    const ElementMapper &elementMapper() const
    { return elementMapper_; }

    /*!
     * \brief Returns the mapper for vertices to indices for possibly adaptive grids.
     */
    VertexMapper &vertexMapper()
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices for possibly adaptive grids.
     */
    ElementMapper &elementMapper()
    { return elementMapper_; }

    /*!
     * \brief Returns the bounding box tree of the grid
     */
    const BoundingBoxTree& boundingBoxTree() const
    { return *boundingBoxTree_; }

    /*!
     * \brief Returns the element index to element map
     */
    const ElementMap& elementMap() const
    { return *elementMap_; }

    /*!
     * \brief Get an element from a global element index
     */
    Element element(GridIndexType eIdx) const
    { return elementMap()[eIdx]; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the smallest values.
     */
    const GlobalCoordinate &bBoxMin() const
    { return bBoxMin_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalCoordinate &bBoxMax() const
    { return bBoxMax_; }

private:

    //! Return an instance of the element mapper
    ElementMapper makeElementMapper_(const GridView& gridView) const
    {
        if constexpr (std::is_constructible<ElementMapper, GridView, Dune::MCMGLayout>())
            return ElementMapper(gridView, Dune::mcmgElementLayout());
        else
            return ElementMapper(gridView);
    }

    //! Return an instance of the vertex mapper
    VertexMapper makeVertexMapper_(const GridView& gridView) const
    {
        if constexpr (std::is_constructible<VertexMapper, GridView, Dune::MCMGLayout>())
            return VertexMapper(gridView, Dune::mcmgVertexLayout());
        else
            return VertexMapper(gridView);
    }

    //! Compute the bounding box of the entire domain, for e.g. setting boundary conditions
    void computeGlobalBoundingBox_()
    {
        // calculate the bounding box of the local partition of the grid view
        for (const auto& vertex : vertices(gridView_))
        {
            for (int i=0; i<dimWorld; i++)
            {
                using std::min;
                using std::max;
                bBoxMin_[i] = min(bBoxMin_[i], vertex.geometry().corner(0)[i]);
                bBoxMax_[i] = max(bBoxMax_[i], vertex.geometry().corner(0)[i]);
            }
        }

        // communicate to get the bounding box of the whole domain
        if (gridView_.comm().size() > 1)
        {
            for (int i = 0; i < dimWorld; ++i)
            {
                bBoxMin_[i] = gridView_.comm().min(bBoxMin_[i]);
                bBoxMax_[i] = gridView_.comm().max(bBoxMax_[i]);
            }
        }
    }

    void update_()
    {
        // Update the mappers
        elementMapper_.update(gridView_);
        vertexMapper_.update(gridView_);

        // Compute the bounding box of the entire domain, for e.g. setting boundary conditions
        computeGlobalBoundingBox_();

        // update element map and bounding box tree
        // always building these comes at a memory overhead but improved
        // performance and thread-safe element level access (e.g. during assembly)
        // for all simulation that use these features
        elementMap_ = std::make_shared<ElementMap>(gridView_.grid(), elementMapper_);
        boundingBoxTree_ = std::make_unique<BoundingBoxTree>(
            std::make_shared<ElementSet>(gridView_, elementMapper(), elementMap_)
        );
    }

    //! the process grid view
    GridView gridView_;

    //! entity mappers
    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    //! the bounding box tree of the grid view for efficient element intersections
    std::unique_ptr<const BoundingBoxTree> boundingBoxTree_;

    //! a map from element index to elements
    std::shared_ptr<const ElementMap> elementMap_;

    //! the bounding box of the whole domain
    GlobalCoordinate bBoxMin_;
    GlobalCoordinate bBoxMax_;
};

} // end namespace Dumux

#endif
