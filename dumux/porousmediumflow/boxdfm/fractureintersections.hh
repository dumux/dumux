// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDFMModel
 * \copybrief Dumux::BoxDfmFractureIntersections
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_FRACTURE_INTERSECTIONS_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_FRACTURE_INTERSECTIONS_HH

#include <type_traits>
#include <algorithm>
#include <optional>
#include <functional>
#include <utility>

#include <dune/common/indices.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/entitymap.hh>
#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief Exposes a visitor over all bulk grid intersections that coincide with fractures.
 *        This implementation uses the facet coupling grid manager under the hood.
 */
template<std::size_t fractureGridId, class... Grids>
class BoxDfmFractureIntersections
{
    using GridManager = FacetCouplingGridManager<Grids...>;
    using GridData = typename GridManager::GridData;

    static_assert(fractureGridId > 0);
    static_assert(fractureGridId < GridManager::numGrids);
    static constexpr auto bulkGridId = fractureGridId - 1;

    using Embeddings = typename GridManager::Embeddings;
    using GridAdapter = CodimOneGridAdapter<Embeddings, bulkGridId, fractureGridId>;

    using FractureGrid = typename GridManager::Grid<fractureGridId>;
    using FractureGridElement = typename FractureGrid::template Codim<0>::Entity;
    using FractureGridVertex = typename FractureGrid::template Codim<FractureGrid::dimension>::Entity;

    using BulkGrid = typename GridManager::Grid<bulkGridId>;
    using BulkGridElement = typename BulkGrid::template Codim<0>::Entity;
    using BulkGridReferenceElement = std::decay_t<decltype(referenceElement(std::declval<const BulkGridElement&>()))>;
    using BulkGridIntersection = typename BulkGrid::LeafGridView::Intersection;
    using BulkElementMap = EntityMap<typename BulkGrid::LeafGridView>;
    using BulkGridIndexType = typename BulkGrid::LeafGridView::IndexSet::IndexType;

public:
    using BarrierMarker = std::function<bool(const FractureGridElement&, const Embeddings&, std::optional<int>)>;

    //! Structure to hold information on a bulk grid intersection coinciding with a fracture
    struct BulkFractureIntersection
    {
        const BulkGridElement& element;
        const typename BulkGridElement::Geometry& elementGeometry;
        const BulkGridReferenceElement& referenceElement;

        const BulkGridIntersection& intersection;
        const typename BulkGridIntersection::Geometry& intersectionGeometry;
        const std::vector<BulkGridIndexType>& intersectionVertexIndices;
    };

    //! Construction from a grid manager and the grid id of the fracture grid.
    BoxDfmFractureIntersections(Dune::index_constant<fractureGridId>,
                                const GridManager& gridManager,
                                BarrierMarker&& isBarrier = [] (auto&&...) { return false; })
    : gridManager_(gridManager)
    , embeddings_{gridManager_.getEmbeddings()}
    , adapter_{embeddings_}
    , bulkElementMap_{makeBulkElementMap_()}
    , gridData_{}
    , isBarrierFunction_{std::move(isBarrier)}
    , barriersTakePrecendence_{getParam<bool>("BoxDFM.BarriersTakePrecedence", false)}
    {
        const auto& bulkGridView = gridManager_.template grid<bulkGridId>().leafGridView();
        bulkInsertionToElementIndex_.resize(bulkGridView.size(0));
        std::size_t i = 0;
        for (const auto& element : elements(bulkGridView))
            bulkInsertionToElementIndex_.at(embeddings_->template insertionIndex<bulkGridId>(element)) = i++;

        if (gridManager_.hasGridData())
            gridData_ = gridManager_.getGridData();
    }

    //! Visit all bulk grid intersections that coincide with barrier
    template<class Visitor>
    void visitBarriers(const Visitor& visitor) const
    {
        static_assert(
            std::is_invocable_v<Visitor, const FractureGridElement&>,
            "Visitor must be invocable with const FractureGridElement&"
        );
        for (const auto& element : elements(fractureGridView_()))
            if (isBarrier_(element))
                visitor(element);
    }

    //! Visit all bulk grid intersections that coincide with barriers
    template<class Visitor>
    void visitBarrierIntersections(const Visitor visitor) const
    { visit_(visitor, [&] (const auto& e) { return isBarrier_(e); }); }


    //! Visit all bulk grid intersections that coincide with fractures
    template<class Visitor>
    void visit(const Visitor& visitor) const
    { visit_(visitor, [&] (const auto& e) { return !isBarrier_(e); }); }

    //! Enrich the given vertex mapper according to the defined barriers
    template<class VertexMapper>
    void enrich(VertexMapper& mapper) const
    {
        mapper.enrich(fractureGridView_(),
                      adapter_,
                      [&] (const auto& e) { return isBarrier_(e); });
    }

private:
    auto fractureGridView_() const
    { return gridManager_.template grid<fractureGridId>().leafGridView(); }

    //! Visit all bulk grid intersections that coincide with fractures
    template<class Visitor, class Filter>
    void visit_(const Visitor& visitor, const Filter& filter) const
    {
        static_assert(
            std::is_invocable_v<Visitor, const BulkFractureIntersection&>,
            "Visitor must be invocable with const BulkFractureIntersection&"
        );

        // predeclare storage to avoid dynamic allocations
        std::vector<BulkGridIndexType> fractureVertexIndices;
        std::vector<BulkGridIndexType> isVertexIndices;
        std::vector<BulkGridIndexType> sortedIsVertexIndices;
        for (const auto& element : elements(fractureGridView_()))
        {
            if (!filter(element))
                continue;

            const auto numElementCorners = element.subEntities(FractureGrid::dimension);
            fractureVertexIndices.resize(numElementCorners);
            for (int i = 0; i < numElementCorners; ++i)
                fractureVertexIndices[i] = bulkGridVertexIndex_(element.template subEntity<FractureGrid::dimension>(i));
            std::sort(fractureVertexIndices.begin(), fractureVertexIndices.end());

            for (int k = 0; k < numAdjoinedBulkElements_(element); ++k)
            {
                const auto bulkElement = adjoinedBulkElement_(element, k);
                const auto bulkElementGeometry = bulkElement.geometry();
                const auto bulkRefElement = referenceElement(bulkElementGeometry);

                bool intersectionFound = false;
                const auto& bulkGridView = gridManager_.template grid<bulkGridId>().leafGridView();
                for (const auto& bulkIntersection : intersections(bulkGridView, bulkElement))
                {
                    const auto bulkIntersectionGeometry = bulkIntersection.geometry();
                    if (bulkIntersectionGeometry.corners() != numElementCorners)
                        continue;

                    isVertexIndices.resize(numElementCorners);
                    for (unsigned int vIdxLocal = 0; vIdxLocal < numElementCorners; ++vIdxLocal)
                        isVertexIndices[vIdxLocal] = adapter_.bulkGridVertexMapper().subIndex(
                            bulkElement,
                            bulkRefElement.subEntity(bulkIntersection.indexInInside(), 1, vIdxLocal, BulkGrid::dimension),
                            BulkGrid::dimension
                        );
                    sortedIsVertexIndices = isVertexIndices;
                    std::sort(sortedIsVertexIndices.begin(), sortedIsVertexIndices.end());

                    if (std::equal(fractureVertexIndices.begin(), fractureVertexIndices.end(), sortedIsVertexIndices.begin()))
                    {
                        intersectionFound = true;
                        visitor(BulkFractureIntersection{
                            bulkElement, bulkElementGeometry, bulkRefElement,
                            bulkIntersection, bulkIntersectionGeometry, isVertexIndices
                        });
                        break;
                    }
                }

                if (!intersectionFound)
                    DUNE_THROW(Dune::InvalidStateException, "Could not find bulk intersection coinciding with fracture");
            }
        }
    }

    bool isBarrier_(const FractureGridElement& fractureElement) const
    {
        std::optional<int> marker;
        if (gridData_)
            marker = gridData_.value()->template getElementDomainMarker<fractureGridId>(fractureElement);
        return isBarrierFunction_(fractureElement, *embeddings_, marker);
    }

    unsigned int numAdjoinedBulkElements_(const FractureGridElement& fractureElement) const
    { return adapter_.numEmbedments(fractureElement); }

    BulkGridElement adjoinedBulkElement_(const FractureGridElement& fractureElement, unsigned int i) const
    {
        const auto insertionIndex = embeddings_->template adjoinedEntityIndices<fractureGridId>(fractureElement).at(i);
        return bulkElementMap_[bulkInsertionToElementIndex_.at(insertionIndex)];
    }

    auto bulkGridVertexIndex_(const FractureGridVertex& fractureVertex) const
    { return adapter_.bulkGridVertexIndex(fractureVertex); }

    BulkElementMap makeBulkElementMap_() const
    {
        const auto& bulkGrid = gridManager_.template grid<bulkGridId>();
        const auto& bulkGridView = bulkGrid.leafGridView();
        std::vector<typename BulkGrid::template Codim<0>::EntitySeed> seeds;
        seeds.reserve(bulkGridView.size(0));
        for (const auto& element : elements(bulkGridView))
            seeds.emplace_back(element.seed());
        return {bulkGrid, std::move(seeds)};
    }

    const GridManager& gridManager_;
    std::shared_ptr<const Embeddings> embeddings_;
    GridAdapter adapter_;
    BulkElementMap bulkElementMap_;
    std::vector<std::size_t> bulkInsertionToElementIndex_;
    std::optional<std::shared_ptr<const GridData>> gridData_;
    BarrierMarker isBarrierFunction_;
    bool barriersTakePrecendence_;
};

} // namespace Dumux

#endif
