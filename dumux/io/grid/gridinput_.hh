// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Grids
 * \brief Grid factory adapter for reading grid data
 */
#ifndef DUMUX_IO_GRID_GRID_FACTORY__HH
#define DUMUX_IO_GRID_GRID_FACTORY__HH

#include <memory>
#include <array>
#include <ranges>
#include <algorithm>
#include <functional>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/capabilities.hh>

#include <dumux/geometry/intersectingentities.hh>

namespace Dumux::Detail::GridData  {

template<class Grid>
struct GridInput
{
    using Intersection = typename Grid::LeafIntersection;
public:
    GridInput(std::shared_ptr<Grid> grid, std::shared_ptr<Dune::GridFactory<Grid>> factory)
    : grid_(std::move(grid))
    , factory_(std::move(factory))
    {
        // assume that all entities are on rank 0
        numElements_ = grid_->levelGridView(0).size(0);
        numVertices_ = grid_->levelGridView(0).size(Grid::dimension);

        if (grid_->comm().size() > 1)
        {
            numElements_ = grid_->comm().sum(numElements_);
            numVertices_ = grid_->comm().sum(numVertices_);
        }
    }

    bool wasInserted(const Intersection& intersection) const
    { return factory_->wasInserted(intersection); }

    template<class Entity>
    auto insertionIndex(const Entity& e) const
    { return factory_->insertionIndex(e); }

    std::size_t numElements() const
    { return numElements_; }

    std::size_t numVertices() const
    { return numVertices_; }

private:
    std::shared_ptr<Grid> grid_;
    std::shared_ptr<Dune::GridFactory<Grid>> factory_;
    std::size_t numElements_, numVertices_;
};

template<class Grid>
concept CartesianGrid = Dune::Capabilities::isCartesian<Grid>::v;

template<CartesianGrid Grid>
struct GridInput<Grid>
{
    using Element = typename Grid::template Codim<0>::Entity;
    using Vertex = typename Grid::template Codim<Grid::dimension>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Intersection = typename Grid::LeafIntersection;
public:

    template<class ImageGrid>
    GridInput(std::shared_ptr<Grid> grid, std::shared_ptr<ImageGrid> factory)
    : grid_(std::move(grid))
    {
        cells_ = vertices_ = factory->cells();
        for (std::size_t i = 0; i < Grid::dimension; ++i)
            vertices_[i] += 1;

        lowerLeft_ = lowerLeftVertices_ = factory->lowerLeft();
        upperRight_ = upperRightVertices_ = factory->upperRight();

        const auto spacing = factory->spacing();
        for (std::size_t i = 0; i < Grid::dimension; ++i)
        {
            lowerLeftVertices_[i] -= 0.5 * spacing[i];
            upperRightVertices_[i] += 0.5 * spacing[i];
        }

        numElements_ = std::accumulate(cells_.begin(), cells_.end(), 1, std::multiplies<int>());
        numVertices_ = std::accumulate(vertices_.begin(), vertices_.end(), 1, std::multiplies<int>());
    }

    bool wasInserted(const Intersection& intersection) const
    { return false; }

    unsigned int insertionIndex(const Element& e) const
    {
        const auto c = e.geometry().center();
        return intersectingEntityCartesianGrid(c, lowerLeft_, upperRight_, cells_);
    }

    unsigned int insertionIndex(const Vertex& v) const
    {
        const auto c = v.geometry().corner(0);
        return intersectingEntityCartesianGrid(c, lowerLeftVertices_, upperRightVertices_, vertices_);
    }

    std::size_t numElements() const { return numElements_; }
    std::size_t numVertices() const { return numVertices_; }

    const GlobalPosition& lowerLeft() const { return lowerLeft_; }
    const GlobalPosition& upperRight() const { return upperRight_; }

    const GlobalPosition& lowerLeftVertices() const { return lowerLeftVertices_; }
    const GlobalPosition& upperRightVertices() const { return upperRightVertices_; }

    const std::array<int, Grid::dimension>& cells() const { return cells_; }
    const std::array<int, Grid::dimension>& vertices() const { return vertices_; }

private:
    std::array<int, Grid::dimension> cells_, vertices_;
    GlobalPosition lowerLeft_, lowerLeftVertices_;
    GlobalPosition upperRight_, upperRightVertices_;

    std::shared_ptr<Grid> grid_;
    std::size_t numElements_, numVertices_;
};

}  // end namespace Dumux::Detail::GridIO

#endif
