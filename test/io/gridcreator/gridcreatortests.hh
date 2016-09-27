/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief Some tests for the grid creator
 * \note This can/should be extended
 */

#ifndef DUMUX_GRIDCREATOR_TESTS_HH
#define DUMUX_GRIDCREATOR_TESTS_HH

#include <dune/common/version.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dumux/io/gridcreator.hh>

namespace Dumux {

template<class Grid>
class GridCreatorTests
{
    using GridView = typename Grid::LeafGridView;
    using Scalar = double;
    static const int dim = Grid::dimension;
    using GridCreator = typename Dumux::GridCreatorImpl<Grid, DiscretizationMethods::None>;
    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

public:

    static void testBoundaryDomainMarkers()
    {
        // initialize the grid
        initialize_();

        // Read the boundary markers and convert them to vertex flags (e.g. for use in a box method)
        // Write a map from vertex position to boundaryMarker
        std::vector<int> boundaryMarker, rank;
        getBoundaryDomainMarkers_(boundaryMarker, rank);

        // construct a vtk output writer and attach the boundaryMakers
        Dune::VTKSequenceWriter<typename Grid::LeafGridView> vtkWriter(GridCreator::grid().leafGridView(), "bifurcation", ".", "");
        vtkWriter.addVertexData(boundaryMarker, "boundaryMarker");
        vtkWriter.addCellData(rank, "rank");
        vtkWriter.write(0);

        // refine grid once. Due to parametrized boundaries this will result in a grid closer to the orginal geometry.
        GridCreator::grid().globalRefine(1);
        getBoundaryDomainMarkers_(boundaryMarker, rank);
        vtkWriter.write(1);
    }

    static void testElementDomainMarkers(const std::string& type = "gmsh",
                                         const std::string& vtkFileName = "test")
    {
        // initialize the grid
        initialize_();

        // Read the boundary markers and convert them to vertex flags (e.g. for use in a box method)
        // Write a map from vertex position to boundaryMarker
        std::vector<int> elementMarker, rank;
        getElementDomainMarkers_(elementMarker, rank, type);

        // construct a vtk output writer and attach the boundaryMakers
        Dune::VTKSequenceWriter<typename Grid::LeafGridView> vtkWriter(GridCreator::grid().leafGridView(), vtkFileName, ".", "");
        vtkWriter.addCellData(elementMarker, "elementMarker");
        vtkWriter.addCellData(rank, "rank");
        vtkWriter.write(0);

        // refine grid once. Check if element markers still fit
        GridCreator::grid().globalRefine(1);
        getElementDomainMarkers_(elementMarker, rank, type);
        vtkWriter.write(1);
    }

private:

    static void getElementDomainMarkers_(std::vector<int>& elementMarker,
                                         std::vector<int>& rank,
                                         const std::string& type)
    {
        const auto& gridView = GridCreator::grid().leafGridView();

        Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper(gridView, Dune::mcmgElementLayout());
        elementMarker.clear();
        elementMarker.resize(gridView.size(0));
        rank.clear();
        rank.resize(gridView.size(0));

        for(const auto& element : elements(gridView))
        {
            auto eIdx = elementMapper.index(element);
            rank[eIdx] = gridView.comm().rank();

            if (type == "gmsh")
                elementMarker[eIdx] = GridCreator::getElementDomainMarker(element);
            else if (type == "dgf")
                elementMarker[eIdx] = GridCreator::parameters(element)[0];
            else
                DUNE_THROW(Dune::InvalidStateException, "No parameter type for " << type);
        }
    }

    static void getBoundaryDomainMarkers_(std::vector<int>& boundaryMarker,
                                          std::vector<int>& rank)
    {
        const auto& gridView = GridCreator::grid().leafGridView();

        Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper(gridView, Dune::mcmgElementLayout());
        Dune::MultipleCodimMultipleGeomTypeMapper<GridView> vertexMapper(gridView, Dune::mcmgVertexLayout());

        boundaryMarker.clear();
        boundaryMarker.resize(gridView.size(dim));
        rank.clear();
        rank.resize(gridView.size(0));

        for(const auto& element : elements(gridView))
        {
            auto eIdx = elementMapper.index(element);
            rank[eIdx] = gridView.comm().rank();
            for(const auto& intersection : intersections(gridView, element))
            {
                if(!intersection.boundary())
                    continue;

                const auto refElement = ReferenceElements::general(element.geometry().type());
                const auto numVertices = refElement.size(intersection.indexInInside(), 1, dim);
                // loop over vertices of the intersection facet
                for(std::size_t vIdx = 0; vIdx < numVertices; vIdx++)
                {
                    // get local vertex index with respect to the element
                    auto vIdxLocal = refElement.subEntity(intersection.indexInInside(), 1, vIdx, dim);
                    auto vIdxGlobal = vertexMapper.subIndex(element, vIdxLocal, dim);

                    // make sure we always take the lowest non-zero marker (problem dependent!)
                    if (boundaryMarker[vIdxGlobal] == 0)
                        boundaryMarker[vIdxGlobal] = GridCreator::getBoundaryDomainMarker(intersection.boundarySegmentIndex());
                    else
                    {
                        if (boundaryMarker[vIdxGlobal] > GridCreator::getBoundaryDomainMarker(intersection.boundarySegmentIndex()))
                            boundaryMarker[vIdxGlobal] = GridCreator::getBoundaryDomainMarker(intersection.boundarySegmentIndex());
                    }
                }
            }
        }
    }

    static void initialize_()
    {
        // Make the grid
        GridCreator::makeGrid();

        // Load balancing if parallel
        GridCreator::loadBalance();
    }
};

} // end namespace Dumux

#endif
