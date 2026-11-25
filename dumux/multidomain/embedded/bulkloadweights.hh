// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Load weights of the bulk elements for partitioning a bulk grid with an embedded lower-dimensional grid
 */
#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_BULK_LOAD_WEIGHTS_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_BULK_LOAD_WEIGHTS_HH

#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include <dune/grid/common/partitionset.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Load weights of the elements of the coarsest level of a bulk grid that include the work of an
 *        embedded coupling
 *
 * A leaf element weighs one plus pieceWeight times the number of pieces of lower-dimensional elements
 * inside it, which is where the coupling places its point sources and samples. An element of the coarsest
 * level weighs the sum of the weights of its interior leaf elements, since grids distribute whole trees of
 * elements. How much work a piece causes relative to a bulk element depends on the coupling mode and the
 * discretizations.
 *
 * \param bulkGridView the bulk grid view of this process
 * \param lowDimGridView the whole lower-dimensional grid
 * \param pieceWeight the weight of a piece relative to a bulk element
 * \return the weight per element of the coarsest level of the bulk grid, indexed by its level index; zero
 *         for elements without interior leaf elements on this process
 */
template<class BulkGridView, class LowDimGridView>
std::vector<int> bulkLoadWeights(const BulkGridView& bulkGridView, const LowDimGridView& lowDimGridView, double pieceWeight)
{
    using BulkSet = GridViewGeometricEntitySet<BulkGridView, 0>;
    using LowDimSet = GridViewGeometricEntitySet<LowDimGridView, 0>;
    const auto bulkSet = std::make_shared<BulkSet>(bulkGridView);

    std::vector<std::size_t> pieces(bulkSet->size(), 0);
    if (bulkSet->size() > 0 && lowDimGridView.size(0) > 0)
    {
        const BoundingBoxTree<BulkSet> bulkTree(bulkSet);
        const BoundingBoxTree<LowDimSet> lowDimTree(std::make_shared<LowDimSet>(lowDimGridView));
        for (const auto& intersection : intersectingEntities(bulkTree, lowDimTree))
            ++pieces[intersection.first()];
    }

    const auto& coarseIndexSet = bulkGridView.grid().levelGridView(0).indexSet();
    std::vector<double> weights(coarseIndexSet.size(0), 0.0);
    for (const auto& element : elements(bulkGridView, Dune::Partitions::interior))
    {
        auto coarse = element;
        while (coarse.hasFather())
            coarse = coarse.father();
        weights[coarseIndexSet.index(coarse)] += 1.0 + pieceWeight*pieces[bulkSet->index(element)];
    }

    std::vector<int> roundedWeights(weights.size());
    for (std::size_t i = 0; i < weights.size(); ++i)
        roundedWeights[i] = static_cast<int>(std::lround(weights[i]));
    return roundedWeights;
}

} // end namespace Dumux

#endif
