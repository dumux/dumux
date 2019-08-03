// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup FacetCoupling
 * \copydoc Dumux::EnrichedGridManager
 */
#ifndef DUMUX_FACETCOUPLING_ENRICHED_GRID_MANAGER
#define DUMUX_FACETCOUPLING_ENRICHED_GRID_MANAGER

#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/io/grid/griddata.hh>
#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Grid manager that creates grids in which the faces coinciding
 *        with a lower-dimensional facet grid are doubled. This way,
 *        interior boundaries are interpreted as actual boundaries in the grid.
 */
template<class GridType>
class EnrichedGridManager
{
    using GridFactory = Dune::GridFactory<GridType>;
    using GV = typename GridType::LeafGridView;
    using Index = typename IndexTraits<GV>::GridIndex;
    using ctype = typename GridType::ctype;

    static constexpr int dim = GridType::dimension;
    using RefElements = Dune::ReferenceElements<ctype, dim>;

public:
    using Grid = GridType;
    using GridData = Dumux::GridData<Grid>;

    /*!
     * \brief Create the grid from a given grid view and vertex dof mapper.
     * \note This expects a vertex dof mapper that fulfills the EnrichedVertexDofMapper
     *       interface. This means that the vertices lying on the lower-dimensional
     *       facet grid already have multiple degrees of freedom assigned to them.
     *
     * \param gridView The grid view
     * \param elementMapper Maps elements to indices
     * \param vertexDofMapper Enriched vertex dof mapper
     * \param grid data (optional). If argument is given, a new grid data
     *        data on the enriched grid will be constructed and stored
     */
    template<class GridView, class ElementMapper, class VertexDofMapper>
    void init(const GridView& gridView,
              const ElementMapper& elementMapper,
              const VertexDofMapper& vertexDofMapper,
              std::shared_ptr<const Dumux::GridData<typename GridView::Grid>> gridData = nullptr)
    {
        static_assert(int(GridView::dimension) == dim, "Provided grid must have the same dimension");

        elementIndexMap_.resize(gridView.size(0));
        vertexIndexMap_.resize(gridView.size(GridView::dimension));

        std::vector<Dune::GeometryType> elemGeometryTypes(gridView.size(0));
        std::vector<std::vector<Index>> elemCornerIndices(gridView.size(0));
        std::vector<bool> dofHandled(vertexDofMapper.size(), false);
        std::vector<Index> dofToInsertionIdx(vertexDofMapper.size());

        std::vector<int> elementMarkers;
        std::vector<int> boundaryMarkers;
        std::vector<std::vector<Index>> boundarySegments;
        if (gridData)
            elementMarkers.resize(gridView.size(0));

        std::size_t elementCount = 0;
        std::size_t vertexCount = 0;

        // insert vertices
        GridFactory gridFactory;
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = elementMapper.index(element);

            elemGeometryTypes[elementCount] = element.geometry().type();
            elementIndexMap_[eIdx] = elementCount;
            if (gridData)
                elementMarkers[elementCount] = gridData->getElementDomainMarker(element);

            auto& corners = elemCornerIndices[eIdx];
            for (unsigned int i = 0; i < element.subEntities(dim); ++i)
            {
                const auto vIdx = vertexDofMapper.vertexIndex(element, i, dim);
                const auto dofIdx = vertexDofMapper.subIndex(element, i, dim);

                if (dofHandled[dofIdx])
                    corners.push_back(dofToInsertionIdx[dofIdx]);
                else
                {
                    const auto& v = element.template subEntity<dim>(i);
                    gridFactory.insertVertex(v.geometry().center());

                    corners.push_back(vertexCount);
                    vertexIndexMap_[vIdx].push_back(vertexCount);
                    dofToInsertionIdx[dofIdx] = vertexCount++;
                    dofHandled[dofIdx] = true;
                }
            }

            // check for boundary segments
            if (gridData)
            {
                const auto refElem = RefElements::general(elemGeometryTypes[elementCount]);

                for (const auto& is : intersections(gridView, element))
                {
                    if (is.boundary() && gridData->wasInserted(is))
                    {
                        std::vector<Index> isCorners;
                        for (unsigned int c = 0; c < is.geometry().corners(); ++c)
                        {
                            const auto vIdxLocal = refElem.subEntity(is.indexInInside(), 1, c, dim);
                            const auto dofIdx = vertexDofMapper.subIndex(element, vIdxLocal, dim);

                            assert(dofHandled[dofIdx]);
                            isCorners.push_back(dofToInsertionIdx[dofIdx]);
                        }

                        boundarySegments.push_back(isCorners);
                        boundaryMarkers.push_back(gridData->getBoundaryDomainMarker(is));
                    }
                }
            }

            // go to next element
            elementCount++;
        }

        // insert elements
        for (std::size_t i = 0; i < elementCount; ++i)
            gridFactory.insertElement(elemGeometryTypes[i], elemCornerIndices[i]);

        // insert boundary segments
        if (gridData)
            for (const auto& bs : boundarySegments)
                gridFactory.insertBoundarySegment(bs);

        gridPtr_ = std::shared_ptr<Grid>(gridFactory.createGrid());
        gridFactory_ = std::make_shared<GridFactory>(gridFactory);
        if (gridData)
            gridDataPtr_ = std::make_shared<GridData>(gridPtr_,
                                                      gridFactory_,
                                                      std::move(elementMarkers),
                                                      std::move(boundaryMarkers));
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    { return *gridPtr_; }

    /*!
     * \brief Call loadBalance() function of the grid.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
        {
            // if we have gmsh parameters we have to manually load balance the data
            if (gridDataPtr_)
            {
                // element and face markers are communicated during load balance
                auto dh = gridDataPtr_->createGmshDataHandle();
                gridPtr_->loadBalance(dh.interface());
                gridPtr_->communicate(dh.interface(), Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
            }
            else
                gridPtr_->loadBalance();
        }
    }

    /*!
     * \brief Return shared pointer to grid data object.
     */
    std::shared_ptr<GridData> getGridData() const
    {
        if (!gridDataPtr_)
            DUNE_THROW(Dune::IOError, "No grid data available");
        return gridDataPtr_;
    }

    /*!
     * \brief Return the insertion index of the element
     *        with index eIdx in the original grid
     */
    Index elementInsertionIndex(std::size_t eIdx)
    { return elementIndexMap_[eIdx]; }

    /*!
     * \brief Return the insertion indices of all vertices
     *        that were inserted for all dofs at the
     *        original vertex with index vIdx
     */
    Index vertexInsertionIndices(std::size_t vIdx)
    { return vertexIndexMap_[vIdx]; }

private:
    std::shared_ptr<Grid> gridPtr_;
    std::shared_ptr<GridData> gridDataPtr_;
    std::shared_ptr<GridFactory> gridFactory_;

    std::vector<Index> elementIndexMap_;
    std::vector<std::vector<Index>> vertexIndexMap_;
};

} // namespace Dumux

#endif
