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
#include <cstddef>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVEModel
 * \brief Maps each coarse-level column to the fine-level elements it contains and each fine-level element to its column
 *
 * The fine-level elements of each column are sorted by increasing vertical coordinate.
 */
template<class GridGeometry, class Scalar>
class TwoPVEColumnMapping
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using Column = std::vector<Element>;

    TwoPVEColumnMapping(std::shared_ptr<const GridGeometry> gridGeometryCoarse,
                        std::shared_ptr<const GridGeometry> gridGeometryFine)
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
    void build_(const GridGeometry& gridGeometryCoarse,
                const GridGeometry& gridGeometryFine)
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

            const auto coarseCandidates = intersectingEntities(fineCenter, gridGeometryCoarse.boundingBoxTree());
            if (coarseCandidates.empty())
                DUNE_THROW(Dune::InvalidStateException, "Fine element " << idxFine << " is not contained in a coarse element");
            if (coarseCandidates.size() != 1)
                DUNE_THROW(Dune::InvalidStateException, "Fine element " << idxFine << " intersects more than one coarse element");
            const auto idxCoarse = coarseCandidates.front();

            coarseToFine_[idxCoarse].push_back(elementFine);
            fineToCoarse_[idxFine] = idxCoarse;
        }

        sortColumns_();
        sanityCheckColumn_(numberOfCoarseElements, numberOfFineElements);
    }

    /*!
     * \brief Sorts the fine-element index within the columns from lower to highest z-coordinate. This sorting is naturally given by the indexing scheme of YASP grid.
     */
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
