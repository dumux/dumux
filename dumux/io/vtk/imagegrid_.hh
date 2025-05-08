// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Grid factory adapter for structured grids
 */
#ifndef DUMUX_IO_VTK_IMAGE_GRID__HH
#define DUMUX_IO_VTK_IMAGE_GRID__HH

#include <memory>
#include <array>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/capabilities.hh>
#include <gridformat/gridformat.hpp>

namespace Dumux::Detail::VTKReader  {

template<class ctype, std::size_t dim>
class ImageGrid
{
public:
    ImageGrid(std::shared_ptr<GridFormat::Reader> grid)
    : grid_(std::move(grid))
    {}

    Dune::FieldVector<ctype, dim> lowerLeft() const
    {
        const auto origin = grid_->origin();
        Dune::FieldVector<ctype, dim> lowerLeft;
        for (std::size_t i = 0; i < dim; ++i)
            lowerLeft[i] = origin[i];
        return lowerLeft;
    }

    Dune::FieldVector<ctype, dim> upperRight() const
    {
        const auto spacing = grid_->spacing();
        const auto extents = grid_->extents();
        auto upperRight = lowerLeft();
        for (std::size_t i = 0; i < dim; ++i)
            upperRight[i] += extents[i] * spacing[i];
        return upperRight;
    }

    std::array<int, dim> cells() const
    {
        std::array<int, dim> cells;
        const auto extents = grid_->extents();
        for (std::size_t i = 0; i < dim; ++i)
            cells[i] = extents[i];
        return cells;
    }

    std::array<double, dim> spacing() const
    {
        std::array<double, dim> spacing;
        const auto spacingArray = grid_->spacing();
        for (std::size_t i = 0; i < dim; ++i)
            spacing[i] = spacingArray[i];
        return spacing;
    }
private:
    std::shared_ptr<GridFormat::Reader> grid_;
};

}  // end namespace Dumux::Detail::VTKReader

#endif
