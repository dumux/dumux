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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup FacetCoupling
 * \copydoc Dumux::FacetCouplingGmshReader
 */
#ifndef DUMUX_FACETCOUPLING_GMSH_READER_HH
#define DUMUX_FACETCOUPLING_GMSH_READER_HH

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <unordered_map>

#include <dune/common/timer.hh>
#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Reads gmsh files where (n-1)-dimensional grids are defined on the faces
 *        or edges of n-dimensional grids.
 *
 * \note Lower-dimensional entities appearing in the grid file are interpreted
 *       either as parts of a lower-dimensional grid living on the sub-entities of
 *       the grids with higher dimension, or as boundary segments. Per default, we
 *       consider all entities as part of the lower-dimensional grids. If you want
 *       to specify boundary segments as well, provide a threshold physical entity
 *       index. All entities with physical entity indices below this threshold will
 *       then be interpreted as boundary segments. Use respective physical entity
 *       indexing in your grid file in that case.
 *
 * \tparam BulkGrid The type of the highest-dimensional grid in the hierachy
 * \tparam numGrids The number of grids to be considered in the hierarchy
 */
template <class BulkGrid, int numGrids>
class FacetCouplingGmshReader
{
    // extract some necessary info from bulk grid
    static constexpr int bulkDim = BulkGrid::dimension;
    static constexpr int bulkDimWorld = BulkGrid::dimensionworld;
    using ctype = typename BulkGrid::ctype;
    using GridIndexType = typename IndexTraits< typename BulkGrid::LeafGridView >::GridIndex;
    using GlobalPosition = Dune::FieldVector<ctype, bulkDimWorld>;

    // determine minimum dimension for which a grid is created
    static constexpr int minGridDim = BulkGrid::dimension - numGrids + 1;
    static_assert(minGridDim >= 1, "Grids with dim < 1 cannot be read!");

    // structure to store data on an element
    using VertexIndexSet = std::vector<GridIndexType>;
    struct ElementData
    {
        Dune::GeometryType gt;
        VertexIndexSet cornerIndices;
    };

public:
    //! Reads the data from a given mesh file
    //! Use this routine if you don't specify boundary segments in the grid file
    void read(const std::string& fileName, bool verbose = false)
    {
        read(fileName, 0, verbose);
    }

    //! Reads the data from a given mesh file
    void read(const std::string& fileName, std::size_t boundarySegThresh, bool verbose = false)
    {
        Dune::Timer watch;
        if (verbose) std::cout << "Opening " << fileName << std::endl;
        std::ifstream gridFile(fileName);
        if (gridFile.fail())
            DUNE_THROW(Dune::InvalidStateException, "Could not open the given .msh file. Make sure it exists");

        // currently we only support version 2 file format
        std::string line;
        std::getline(gridFile, line);
        if (line.find("$MeshFormat") == std::string::npos)
            DUNE_THROW(Dune::InvalidStateException, "Expected $MeshFormat in the first line of the grid file!");
        std::getline(gridFile, line);
        if (line.find_first_of("2") != 0)
            DUNE_THROW(Dune::InvalidStateException, "Currently only Gmsh mesh file format version 2 is supported!");

        // read file until we get to the list of nodes
        while (line.find("$Nodes") == std::string::npos)
            std::getline(gridFile, line);

        // read all vertices
        std::getline(gridFile, line);
        const auto numVertices = convertString<std::size_t>(line);
        gridVertices_.resize(numVertices);

        std::getline(gridFile, line);
        std::size_t vertexCount = 0;
        while (line.find("$EndNodes")  == std::string::npos)
        {
            // drop first entry in line (vertex index) and read coordinates
            std::istringstream stream(line);
            std::string buf; stream >> buf;
            GlobalPosition v;
            for (auto& coord : v)
            {
                stream >> coord;
                if (stream.fail()) DUNE_THROW(Dune::IOError, "Could not read vertex coordinate");
            }

            // insert vertex into container and move to next line
            gridVertices_[vertexCount++] = v;
            std::getline(gridFile, line);
        }

        // we should always find as many vertices as the mesh file states
        if (vertexCount != numVertices)
            DUNE_THROW(Dune::InvalidStateException, "Couldn't find as many vertices as stated in the .msh file");

        // read file until we get to the list of elements
        while(line.find("$Elements") == std::string::npos)
            std::getline(gridFile, line);

        // read elements
        std::getline(gridFile, line);
        const auto numElements = convertString<std::size_t>(line);

        // keep track of number of elements
        std::array<std::size_t, numGrids> elementCount;
        std::fill(elementCount.begin(), elementCount.end(), 0);

        // Construct maps that map bulk grid vertex
        // indices to lowDim vertex indices. -1 indicates non-initialized status
        std::size_t elemCount = 0;
        std::array<std::size_t, numGrids> gridVertexCount;
        std::array<std::vector<GridIndexType>, numGrids> gridVertexMap;
        std::array<std::vector<bool>, numGrids> idxIsAssigned;
        std::fill(gridVertexCount.begin(), gridVertexCount.end(), 0);
        std::fill(gridVertexMap.begin(), gridVertexMap.end(), std::vector<GridIndexType>(vertexCount));
        std::fill(idxIsAssigned.begin(), idxIsAssigned.end(), std::vector<bool>(vertexCount, false));
        std::getline(gridFile, line);
        while (line.find("$EndElements") == std::string::npos)
        {
            // pass all indices into vector
            std::istringstream stream(line);
            std::string buf;
            std::vector<std::size_t> lineData;
            while (stream >> buf) lineData.push_back(convertString<std::size_t>(buf));
            assert(lineData.size() >= 4 && "Grid format erroneous or unsupported");

            // obtain geometry type
            const auto gt = obtainGeometryType( lineData[1] );
            const std::size_t physicalIndex = lineData[3];
            const auto geoDim = gt.dim();
            const bool isBoundarySeg = geoDim != bulkDim && physicalIndex < boundarySegThresh;
            if (geoDim >= minGridDim-1)
            {
                // insert boundary segment
                if ((isBoundarySeg || geoDim == minGridDim-1))
                {
                    const unsigned int nextLevelGridIdx = bulkDim-geoDim-1;

                    VertexIndexSet corners;
                    auto it = lineData.begin()+2+lineData[2]+1;
                    for (; it != lineData.end(); ++it)
                    {
                        *it -= 1; // gmsh indices start from 1

                        // insert map if vertex is not inserted yet
                        if (!idxIsAssigned[nextLevelGridIdx][*it])
                        {
                            gridVertexMap[nextLevelGridIdx][*it] = gridVertexCount[nextLevelGridIdx]++;
                            idxIsAssigned[nextLevelGridIdx][*it] = true;
                            vertexIndices_[nextLevelGridIdx].push_back(*it);
                        }

                        corners.push_back(gridVertexMap[nextLevelGridIdx][*it]);
                    }

                    // marker = physical entity index
                    boundaryMarkerMaps_[nextLevelGridIdx].push_back(physicalIndex);
                    boundarySegments_[nextLevelGridIdx].push_back(corners);
                }

                // insert element
                else
                {
                    const unsigned int gridIdx = bulkDim-geoDim;

                    VertexIndexSet corners;
                    auto it = lineData.begin()+2+lineData[2]+1;
                    for (; it != lineData.end(); ++it)
                    {
                        *it -= 1; // gmsh indices start from 1

                        // insert map if vertex is not inserted yet
                        if (!idxIsAssigned[gridIdx][*it])
                        {
                            gridVertexMap[gridIdx][*it] = gridVertexCount[gridIdx]++;
                            idxIsAssigned[gridIdx][*it] = true;
                            vertexIndices_[gridIdx].push_back(*it);
                        }

                        corners.push_back(gridVertexMap[gridIdx][*it]);
                    }

                    // add data to embedments/embeddings
                    if (geoDim > minGridDim)
                    {
                        const auto gridElemCount = elementData_[gridIdx].size();
                        const auto& embeddedVIndices = vertexIndices_[gridIdx+1];
                        const auto& embeddedIndicesAssigned = idxIsAssigned[gridIdx+1];

                        VertexIndexSet cornerIndicesGlobal(corners.size());
                        for (unsigned int i = 0; i < corners.size(); ++i)
                            cornerIndicesGlobal[i] = vertexIndices_[gridIdx][corners[i]];
                        addEmbeddings(cornerIndicesGlobal, gridIdx, gridElemCount, embeddedVIndices, embeddedIndicesAssigned);
                    }

                    // ensure dune-specific corner ordering
                    reorder(gt, corners);

                    // insert element data to grid's container
                    elementMarkerMaps_[gridIdx].push_back(physicalIndex);
                    elementData_[gridIdx].emplace_back(ElementData({gt, corners}));
                }
            }

            // get next line
            std::getline(gridFile, line);
            elemCount++;
        }

        // make sure we read all elements
        if (elemCount != numElements)
            DUNE_THROW(Dune::InvalidStateException, "Didn't read as many elements as stated in the .msh file");

        if (verbose)
        {
            std::cout << "Finished reading gmsh file" << std::endl;
            for (std::size_t id = 0; id < numGrids; ++id)
            {
                std::cout << elementData_[id].size() << " "
                          << bulkDim-id << "-dimensional elements comprising of "
                          << gridVertexCount[id] << " vertices";
                if (id < numGrids-1) std::cout << "," << std::endl;
            }
            std::cout << " have been read in " << watch.elapsed() << " seconds." << std::endl;
        }
    }

    //! Returns the vector with all grid vertices (entire hierarchy)
    const std::vector<GlobalPosition>& gridVertices() const
    { return gridVertices_; }

    //! Returns a grid's vertex indices
    VertexIndexSet& vertexIndices(std::size_t id)
    {
        assert(id < numGrids && "Index exceeds number of grids provided");
        return vertexIndices_[id];
    }

    //! Returns the vector of read elements for a grid
    const std::vector<ElementData>& elementData(std::size_t id) const
    {
        assert(id < numGrids && "Index exceeds number of grids provided");
        return elementData_[id];
    }

    //! Returns the vector of read elements for a grid
    const std::vector<VertexIndexSet>& boundarySegmentData(std::size_t id) const
    {
        assert(id < numGrids && "Index exceeds number of grids provided");
        return boundarySegments_[id];
    }

    //! Returns the maps of element markers
    std::vector<int>& elementMarkerMap(std::size_t id)
    {
        assert(id < numGrids && "Index exceeds number of grids provided");
        return elementMarkerMaps_[id];
    }

    //! Returns the maps of domain markers
    std::vector<int>& boundaryMarkerMap(std::size_t id)
    {
        assert(id < numGrids && "Index exceeds number of grids provided");
        return boundaryMarkerMaps_[id];
    }

    //! Returns the maps of the embedded entities
    std::unordered_map< GridIndexType, std::vector<GridIndexType> >& embeddedEntityMap(std::size_t id)
    {
        assert(id < numGrids && "Index exceeds number of grids provided");
        return embeddedEntityMaps_[id];
    }

    //! Returns the maps of the embedments
    std::unordered_map< GridIndexType, std::vector<GridIndexType> >& adjoinedEntityMap(std::size_t id)
    {
        assert(id < numGrids && "Index exceeds number of grids provided");
        return adjoinedEntityMaps_[id];
    }

private:
    //! converts a value contained in a string
    template<class T>
    T convertString(const std::string& string) const
    {
        T value;
        std::istringstream stream(string);
        stream >> value;
        if (stream.fail())
            DUNE_THROW(Dune::InvalidStateException, "Couldn't convert string: " << string << "to type: " << typeid(T).name());
        return value;
    }

    //! obtain Dune::GeometryType from a given gmsh element type
    Dune::GeometryType obtainGeometryType(std::size_t gmshElemType) const
    {
        // TODO: Version check with Dune 2.5!
        switch (gmshElemType)
        {
            case 15: return Dune::GeometryTypes::vertex;        // points
            case 1:  return Dune::GeometryTypes::line;          // lines
            case 2:  return Dune::GeometryTypes::triangle;      // triangle
            case 3:  return Dune::GeometryTypes::quadrilateral; // quadrilateral
            case 4:  return Dune::GeometryTypes::tetrahedron;   // tetrahedron
            case 5:  return Dune::GeometryTypes::hexahedron;    // hexahedron
            default:
                DUNE_THROW(Dune::NotImplemented, "FacetCoupling gmsh reader for gmsh element type " << gmshElemType);
        }
    }

    //! reorders in a dune way a set of given element corners in gmsh ordering
    void reorder(const Dune::GeometryType gt, VertexIndexSet& cornerIndices) const
    {
        // triangles, lines & tetrahedra need no reordering
        if (gt == Dune::GeometryTypes::hexahedron)
        {
            using std::swap;
            assert(cornerIndices.size() == 8);
            swap(cornerIndices[2], cornerIndices[3]);
            swap(cornerIndices[6], cornerIndices[7]);
        }
        else if (gt == Dune::GeometryTypes::quadrilateral)
        {
            using std::swap;
            assert(cornerIndices.size() == 4);
            swap(cornerIndices[2], cornerIndices[3]);
        }
    }

    //! adds embeddings/embedments to the map for a given element
    void addEmbeddings(const VertexIndexSet& globalCornerIndices,
                       unsigned int gridIdx,
                       std::size_t curElemIdx,
                       const std::vector<GridIndexType>& embeddedVIndices,
                       const std::vector<bool>& embededIndexIsAssigned)
    {
        const unsigned int embeddedGridIdx = gridIdx+1;

        // check for embedments only if any of the corners exist in the embedded grid
        for (auto globalCornerIdx : globalCornerIndices)
        {
            if (embededIndexIsAssigned[globalCornerIdx])
            {
                for (std::size_t i = 0; i < elementData_[embeddedGridIdx].size(); ++i)
                {
                    const auto& e = elementData_[embeddedGridIdx][i];

                    auto vertIsContained = [&embeddedVIndices, &globalCornerIndices] (auto eCornerIdx)
                    {
                        return std::find(globalCornerIndices.begin(),
                                         globalCornerIndices.end(),
                                         embeddedVIndices[eCornerIdx]) != globalCornerIndices.end();
                    };

                    // if all corners are contained within this element, it is embedded
                    if ( std::all_of(e.cornerIndices.begin(), e.cornerIndices.end(), vertIsContained) )
                    {
                        embeddedEntityMaps_[gridIdx][curElemIdx].push_back(i);
                        adjoinedEntityMaps_[embeddedGridIdx][i].push_back(curElemIdx);
                    }
                }

                return;
            }
        }
    }

    //! data on grid entities
    std::vector<GlobalPosition> gridVertices_;
    std::array<VertexIndexSet, numGrids> vertexIndices_;
    std::array<std::vector<ElementData>, numGrids> elementData_;
    std::array<std::vector<VertexIndexSet>, numGrids> boundarySegments_;

    //! data on connectivity between the grids
    std::array< std::unordered_map< GridIndexType, std::vector<GridIndexType> >, numGrids > embeddedEntityMaps_;
    std::array< std::unordered_map< GridIndexType, std::vector<GridIndexType> >, numGrids > adjoinedEntityMaps_;

    //! data on domain and boundary markers
    std::array< std::vector<int>, numGrids > elementMarkerMaps_;
    std::array< std::vector<int>, numGrids > boundaryMarkerMaps_;
};

} // end namespace Dumux

#endif
