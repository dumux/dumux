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
 *
 * \brief Some tests for the grid creator
 * \note This can/should be extended
 */

#ifndef DUMUX_TEST_IO_GRIDMANAGER_TESTS_HH
#define DUMUX_TEST_IO_GRIDMANAGER_TESTS_HH

#include <dune/common/version.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \brief data handle for parallel communication which takes
 *        the minimum non-zero values that are attached to vertices
 */
template <class GridView>
class VertexHandleNonZeroMin
    : public Dune::CommDataHandleIF< VertexHandleNonZeroMin<GridView>, int >
{
    using Container = std::vector<int>;

public:
    VertexHandleNonZeroMin(Container &container, const GridView &gridView)
    : gridView_(gridView)
    , container_(container)
    {}

    bool contains(int dim, int codim) const
    {
        // only communicate vertices
        return codim == dim;
    }

    bool fixedSize(int dim, int codim) const
    {
        // for each vertex we communicate a single field vector which
        // has a fixed size
        return true;
    }

    template<class EntityType>
    size_t size (const EntityType &e) const
    {
        // communicate a field type per entity
        return 1;
    }

    template<class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp &buff, const EntityType &e) const
    {
        int vIdx = gridView_.indexSet().index(e);
        buff.write(container_[vIdx]);
    }

    template<class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        int vIdx = gridView_.indexSet().index(e);
        int tmp;
        buff.read(tmp);
        using std::min;
        if (tmp > 0)
        {
            if (container_[vIdx] == 0)
                container_[vIdx] = tmp;
            else
                container_[vIdx] = min(container_[vIdx], tmp);
        }
    }

private:
    const GridView gridView_;
    Container &container_;
};

template<class Grid>
class GridManagerTests
{
    using GridView = typename Grid::LeafGridView;
    using Scalar = double;
    static const int dim = Grid::dimension;
    using GridManager = typename Dumux::GridManager<Grid>;

public:

    static void testBoundaryAndElementMarkers(const std::string& type = "gmsh",
                                              const std::string& vtkFileName = "test",
                                              bool refine = true)
    {
        // make the grid manager and initialize the grid
        GridManager gridManager;
        gridManager.init();
        auto gridData = gridManager.getGridData();
        const auto& leafGridView = gridManager.grid().leafGridView();

        // read the boundary and element markers as well as the rank
        std::vector<int> boundaryMarker, elementMarker, rank;
        getBoundaryMarkers_(leafGridView, gridData, boundaryMarker);
        getElementMarkers_(leafGridView, gridData, elementMarker, type);
        getRank_(leafGridView, rank);

        // construct a vtk output writer and attach the boundary markers
        Dune::VTKSequenceWriter<typename Grid::LeafGridView> vtkWriter(leafGridView, vtkFileName, ".", "");
        vtkWriter.addVertexData(boundaryMarker, "boundaryMarker");
        vtkWriter.addCellData(elementMarker, "elementMarker");
        vtkWriter.addCellData(rank, "rank");
        vtkWriter.write(0);

        if (refine)
        {
            // refine grid once and write out the markers again
            gridManager.grid().globalRefine(1);
            getBoundaryMarkers_(leafGridView, gridData, boundaryMarker);
            getElementMarkers_(leafGridView, gridData, elementMarker, type);
            getRank_(leafGridView, rank);
            vtkWriter.write(1);
        }
    }

    static void testElementMarkers(const std::string& type = "gmsh",
                                   const std::string& vtkFileName = "test",
                                   bool refine = true)
    {
        // make the grid manager and initialize the grid
        GridManager gridManager;
        gridManager.init();
        auto gridData = gridManager.getGridData();

        // read the element markers and the rank
        std::vector<int> elementMarker, rank;
        getElementMarkers_(gridManager.grid().leafGridView(), gridData, elementMarker, type);
        getRank_(gridManager.grid().leafGridView(), rank);

        // construct a vtk output writer and attach the element markers
        Dune::VTKSequenceWriter<typename Grid::LeafGridView> vtkWriter(gridManager.grid().leafGridView(), vtkFileName, ".", "");
        vtkWriter.addCellData(elementMarker, "elementMarker");
        vtkWriter.addCellData(rank, "rank");
        vtkWriter.write(0);

        if (refine)
        {
            // refine grid once and write out the markers again
            gridManager.grid().globalRefine(1);
            getElementMarkers_(gridManager.grid().leafGridView(), gridData, elementMarker, type);
            getRank_(gridManager.grid().leafGridView(), rank);
            vtkWriter.write(1);
        }
    }

    static void testVertexMarkers(const std::string& type = "dgf",
                                  const std::string& vtkFileName = "test")
    {
        // make the grid manager and initialize the grid
        GridManager gridManager;
        gridManager.init();
        auto gridData = gridManager.getGridData();

        // read the element markers and the rank
        std::vector<int> vertexMarker;
        getVertexMarkers_(gridManager.grid().levelGridView(0), gridData, vertexMarker, type);

        // construct a vtk output writer and attach the element markers
        Dune::VTKSequenceWriter<typename Grid::LevelGridView> vtkWriter(gridManager.grid().levelGridView(0), vtkFileName, ".", "");
        vtkWriter.addVertexData(vertexMarker, "vertexData");
        vtkWriter.write(0);
    }

private:

    static void getRank_(const GridView& gridView, std::vector<int>& rank)
    {
        rank.clear();
        rank.resize(gridView.size(0));

        for(const auto& element : elements(gridView))
        {
            auto eIdx = gridView.indexSet().index(element);
            rank[eIdx] = gridView.comm().rank();
        }
    }

    template<class GridData>
    static void getElementMarkers_(const GridView& gridView,
                                   const GridData& gridData,
                                   std::vector<int>& elementMarker,
                                   const std::string& type)
    {
        elementMarker.clear();
        elementMarker.resize(gridView.size(0));

        for(const auto& element : elements(gridView))
        {
            const auto eIdx = gridView.indexSet().index(element);
            if (type == "gmsh")
                elementMarker[eIdx] = gridData->getElementDomainMarker(element);
            else if (type == "dgf")
                elementMarker[eIdx] = gridData->parameters(element)[0];
            else
                DUNE_THROW(Dune::InvalidStateException, "No parameters for type " << type);
        }
    }

    template<class LevelGridView, class GridData>
    static void getVertexMarkers_(const LevelGridView& gridView,
                                  const GridData& gridData,
                                  std::vector<int>& vertexMarker,
                                  const std::string& type)
    {
        if (type != "dgf")
            DUNE_THROW(Dune::InvalidStateException, "Vertex marker only exist for dgf grids.");

        vertexMarker.clear();
        vertexMarker.resize(gridView.size(Grid::dimension));

        for (const auto& vertex : vertices(gridView.grid().levelGridView(0)))
        {
            const auto vIdx = gridView.indexSet().index(vertex);
            vertexMarker[vIdx] = gridData->parameters(vertex)[0];
        }
    }

    template<class GridData>
    static void getBoundaryMarkers_(const GridView& gridView,
                                    const GridData& gridData,
                                    std::vector<int>& boundaryMarker)
    {
        boundaryMarker.clear();
        boundaryMarker.resize(gridView.size(dim), 0);

        for(const auto& element : elements(gridView))
        {
            for(const auto& intersection : intersections(gridView, element))
            {
                if (!intersection.boundary())
                    continue;

                const auto refElement = referenceElement(element.geometry());
                const auto numVertices = refElement.size(intersection.indexInInside(), 1, dim);
                // loop over vertices of the intersection facet
                for(std::size_t vIdx = 0; vIdx < numVertices; vIdx++)
                {
                    // get local vertex index with respect to the element
                    auto vIdxLocal = refElement.subEntity(intersection.indexInInside(), 1, vIdx, dim);
                    auto vIdxGlobal = gridView.indexSet().subIndex(element, vIdxLocal, dim);

                    // make sure we always take the lowest non-zero marker (problem dependent!)
                    if (boundaryMarker[vIdxGlobal] == 0)
                        boundaryMarker[vIdxGlobal] = gridData->getBoundaryDomainMarker(intersection);
                    else
                    {
                        if (boundaryMarker[vIdxGlobal] > gridData->getBoundaryDomainMarker(intersection))
                            boundaryMarker[vIdxGlobal] = gridData->getBoundaryDomainMarker(intersection);
                    }
                }
            }
        }

        // In a parallel setting, it is possible that not all boundary vertices
        // will be reached on each process by the procedure above that loops
        // over the boundary intersections and their vertices. This mandates the
        // following synchronization.
        if (gridView.comm().size() > 1)
        {
            VertexHandleNonZeroMin<GridView> dataHandle(boundaryMarker, gridView);
            gridView.communicate(dataHandle,
                                 Dune::InteriorBorder_All_Interface,
                                 Dune::ForwardCommunication);
        }
    }
};

} // end namespace Dumux

#endif
