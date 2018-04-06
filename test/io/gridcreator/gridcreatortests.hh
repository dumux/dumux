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
#include <dumux/discretization/methods.hh>

namespace Dumux {

template<class Grid>
class GridCreatorTests
{
    using GridView = typename Grid::LeafGridView;
    using Scalar = double;
    static const int dim = Grid::dimension;
    using GridCreator = typename Dumux::GridCreatorImpl<Grid, DiscretizationMethod::none>;
    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

public:

    static void testBoundaryAndElementMarkers(const std::string& type = "gmsh",
                                              const std::string& vtkFileName = "test")
    {
        // initialize the grid
        initialize_();

        // read the boundary and element markers as well as the rank
        std::vector<int> boundaryMarker, elementMarker, rank;
        getBoundaryMarkers_(boundaryMarker);
        getElementMarkers_(elementMarker, type);
        getRank_(rank);

        // construct a vtk output writer and attach the boundary markers
        Dune::VTKSequenceWriter<typename Grid::LeafGridView> vtkWriter(GridCreator::grid().leafGridView(), vtkFileName, ".", "");
        vtkWriter.addVertexData(boundaryMarker, "boundaryMarker");
        vtkWriter.addCellData(elementMarker, "elementMarker");
        vtkWriter.addCellData(rank, "rank");
        vtkWriter.write(0);

        // refine grid once and write out the markers again
        GridCreator::grid().globalRefine(1);
        getBoundaryMarkers_(boundaryMarker);
        getElementMarkers_(elementMarker, type);
        getRank_(rank);
        vtkWriter.write(1);
    }

    static void testElementMarkers(const std::string& type = "gmsh",
                                   const std::string& vtkFileName = "test")
    {
        // initialize the grid
        initialize_();

        // read the element markers and the rank
        std::vector<int> elementMarker, rank;
        getElementMarkers_(elementMarker, type);
        getRank_(rank);

        // construct a vtk output writer and attach the element markers
        Dune::VTKSequenceWriter<typename Grid::LeafGridView> vtkWriter(GridCreator::grid().leafGridView(), vtkFileName, ".", "");
        vtkWriter.addCellData(elementMarker, "elementMarker");
        vtkWriter.addCellData(rank, "rank");
        vtkWriter.write(0);

        // refine grid once and write out the markers again
        GridCreator::grid().globalRefine(1);
        getElementMarkers_(elementMarker, type);
        getRank_(rank);
        vtkWriter.write(1);
    }

private:

    static void getRank_(std::vector<int>& rank)
    {
        const auto& gridView = GridCreator::grid().leafGridView();

        rank.clear();
        rank.resize(gridView.size(0));

        for(const auto& element : elements(gridView))
        {
            auto eIdx = gridView.indexSet().index(element);
            rank[eIdx] = gridView.comm().rank();
        }
    }

    static void getElementMarkers_(std::vector<int>& elementMarker,
                                   const std::string& type)
    {
        const auto& gridView = GridCreator::grid().leafGridView();

        elementMarker.clear();
        elementMarker.resize(gridView.size(0));

        for(const auto& element : elements(gridView))
        {
            auto eIdx = gridView.indexSet().index(element);

            if (type == "gmsh")
                elementMarker[eIdx] = GridCreator::getElementDomainMarker(element);
            else if (type == "dgf")
                elementMarker[eIdx] = GridCreator::parameters(element)[0];
            else
                DUNE_THROW(Dune::InvalidStateException, "No parameters for type " << type);
        }
    }

    static void getBoundaryMarkers_(std::vector<int>& boundaryMarker)
    {
        const auto& gridView = GridCreator::grid().leafGridView();

        boundaryMarker.clear();
        boundaryMarker.resize(gridView.size(dim), 0);

        for(const auto& element : elements(gridView))
        {
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
                    auto vIdxGlobal = gridView.indexSet().subIndex(element, vIdxLocal, dim);

                    // make sure we always take the lowest non-zero marker (problem dependent!)
                    if (boundaryMarker[vIdxGlobal] == 0)
                        boundaryMarker[vIdxGlobal] = GridCreator::getBoundaryDomainMarker(intersection);
                    else
                    {
                        if (boundaryMarker[vIdxGlobal] > GridCreator::getBoundaryDomainMarker(intersection))
                            boundaryMarker[vIdxGlobal] = GridCreator::getBoundaryDomainMarker(intersection);
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
