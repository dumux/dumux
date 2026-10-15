// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVEModel
 * \brief Manages the maps between the coarse and fine level.
 */

#ifndef DUMUX_TWOPVE_COLUMN_MAPPING_HH
#define DUMUX_TWOPVE_COLUMN_MAPPING_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVEModel
 * \brief Maps horizontal coarse cells to vertically stacked fine cells using horizontal projection
 *
 * The fine-level elements of each column are sorted by increasing vertical coordinate.
 */
template<class CoarseGridGeometry, class FineGridGeometry, class Scalar>
class TwoPVEColumnMapping
{
    using GridView = typename FineGridGeometry::GridView;
    using CoarseGridView = typename CoarseGridGeometry::GridView;
    static_assert(CoarseGridView::dimension + 1 == GridView::dimension);
    static_assert(CoarseGridView::dimensionworld + 1 == GridView::dimensionworld);
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using Column = std::vector<Element>;

    TwoPVEColumnMapping(std::shared_ptr<const CoarseGridGeometry> gridGeometryCoarse,
                        std::shared_ptr<const FineGridGeometry> gridGeometryFine)
    {
        build_(*gridGeometryCoarse, *gridGeometryFine);
    }

    /*!
     * \brief Returns the column of fine-level elements given a coarse-level index
     *
     * \param idxCoarse coarse-level index of the grid column
     */
    const Column& column(std::size_t idxCoarse) const
    {
        return coarseToFine_.at(idxCoarse);
    }

    /*!
     * \brief Returns the coarse-level column index the given fine-level index belongs to
     *
     * \param idxFine index of a fine-level elements
     */
    std::size_t coarseIndex(std::size_t idxFine) const
    {
        return fineToCoarse_.at(idxFine);
    }

    /*!
     * \brief Returns the number of coarse-level columns
     */
    std::size_t numberOfColumns() const noexcept
    {
        return coarseToFine_.size();
    }

private:

    /*!
     * \brief Builds and stores the coarse-to-fine and fine-to-coarse maps
     *
     * \param gridGeometryCoarse coarse-level grid geometry
     * \param gridGeometryFine   fine-level grid geometry
     */
    void build_(const CoarseGridGeometry& gridGeometryCoarse,
                const FineGridGeometry& gridGeometryFine)
    {
        const std::size_t numberOfCoarseElements = gridGeometryCoarse.elementMapper().size();
        const std::size_t numberOfFineElements = gridGeometryFine.elementMapper().size();

        if (numberOfCoarseElements == 0)
            DUNE_THROW(Dune::InvalidStateException, "Cannot construct columns for an empty coarse grid");

        coarseToFine_.clear();
        coarseToFine_.resize(numberOfCoarseElements);

        fineToCoarse_.clear();
        fineToCoarse_.resize(numberOfFineElements);

        for (const auto& elementFine : elements(gridGeometryFine.gridView()))
        {
            const auto idxFine = gridGeometryFine.elementMapper().index(elementFine);
            const auto fineCenter = elementFine.geometry().center();

            Dune::FieldVector<Scalar, CoarseGridView::dimensionworld> horizontalCenter;
            for (int axis = 0; axis < CoarseGridView::dimensionworld; ++axis)
                horizontalCenter[axis] = fineCenter[axis];
            const auto coarseCandidates = intersectingEntities(horizontalCenter, gridGeometryCoarse.boundingBoxTree());
            if (coarseCandidates.empty())
                DUNE_THROW(Dune::InvalidStateException, "Fine element " << idxFine << " is not contained in a coarse element");
            if (coarseCandidates.size() != 1)
                DUNE_THROW(Dune::InvalidStateException, "Fine element " << idxFine << " intersects more than one coarse element");
            const auto idxCoarse = coarseCandidates.front();
            const auto coarseGeometry = gridGeometryCoarse.element(idxCoarse).geometry();
            for (int axis = 0; axis < CoarseGridView::dimensionworld; ++axis)
            {
                const auto coarseBounds = axisBounds_(coarseGeometry, axis);
                const auto fineBounds = axisBounds_(elementFine.geometry(), axis);
                const Scalar tolerance = 1e-10*(coarseBounds[1] - coarseBounds[0]);
                using std::abs;
                if (abs(coarseBounds[0] - fineBounds[0]) > tolerance || abs(coarseBounds[1] - fineBounds[1]) > tolerance)
                    DUNE_THROW(Dune::InvalidStateException, "Coarse and fine grids must have matching horizontal cell bounds");
            }

            coarseToFine_[idxCoarse].push_back(elementFine);
            fineToCoarse_[idxFine] = idxCoarse;
        }

        sortColumns_();
        sanityCheckColumn_(numberOfCoarseElements, numberOfFineElements);
    }

    //! Bounds of an axis-aligned cell along one coordinate.
    template<class Geometry>
    static std::array<Scalar, 2> axisBounds_(const Geometry& geometry, int axis)
    {
        std::array<Scalar, 2> bounds{geometry.corner(0)[axis], geometry.corner(0)[axis]};
        for (int corner = 1; corner < geometry.corners(); ++corner)
        {
            bounds[0] = std::min(bounds[0], geometry.corner(corner)[axis]);
            bounds[1] = std::max(bounds[1], geometry.corner(corner)[axis]);
        }
        return bounds;
    }

    //! Sort fine elements by increasing height within each column.
    void sortColumns_()
    {
        static constexpr int verticalAxis = GridView::dimensionworld - 1;
        for (auto& column : coarseToFine_)
        {
            std::sort(column.begin(), column.end(), [](const Element& a, const Element& b)
            {
                return a.geometry().center()[verticalAxis] < b.geometry().center()[verticalAxis];
            });
        }
    }

    /*!
     * \brief Checks a few conditions for each column to make sure they are in valid state
     *
     * \param numberOfCoarseElements number of coarse-level elements in grid
     * \param numberOfFineElements   number of fine-level elements in grid
     */
    void sanityCheckColumn_(const std::size_t numberOfCoarseElements,
                            const std::size_t numberOfFineElements)
    {
        static constexpr int verticalAxis = GridView::dimensionworld - 1;
        if (numberOfFineElements == 0 || numberOfFineElements % numberOfCoarseElements != 0)
            DUNE_THROW(Dune::InvalidStateException, "Each coarse element must map to a nonempty column of equal size");
        const auto expectedColumnSize = numberOfFineElements / numberOfCoarseElements;

        for (std::size_t coarseIdx = 0; coarseIdx < coarseToFine_.size(); ++coarseIdx)
        {
            const auto& column = coarseToFine_[coarseIdx];

            // check size
            if (column.size() != expectedColumnSize)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected number of fine elements in column " << coarseIdx);

            // check ordering
            for (std::size_t i = 1; i < column.size(); ++i)
            {
                const Scalar previousZ = column[i - 1].geometry().center()[verticalAxis];
                const Scalar currentZ = column[i].geometry().center()[verticalAxis];

                if (!(previousZ < currentZ))
                    DUNE_THROW(Dune::InvalidStateException, "Fine elements do not have strictly increasing vertical coordinates in column " << coarseIdx);
            }
        }
    }

    std::vector<Column> coarseToFine_;
    std::vector<std::size_t> fineToCoarse_;
};

} // end namespace Dumux

#endif
