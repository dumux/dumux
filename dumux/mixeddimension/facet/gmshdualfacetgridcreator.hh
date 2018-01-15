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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Reads gmsh files where fractures are defined as lines in surfaces
 *        or surfaces in volumes and constructs a grid for the bulk part and a
 *        lower-dimensional grid for the lines/surfaces.
 */

#ifndef DUMUX_GMSHDUALFACETGRIDCREATOR_HH
#define DUMUX_GMSHDUALFACETGRIDCREATOR_HH

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/common/math.hh>
#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
}

/*!
 * \brief Reads gmsh files where fractures are defined as lines in surfaces
 *        or surfaces in volumes and constructs a grid for the bulk part and a
 *        lower-dimensional grid for the lines/surfaces.
 */
template <class TypeTag>
class GmshDualFacetGridCreator
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    // obtain the grid types of the sub problems
    using BulkGrid = typename GET_PROP_TYPE(BulkProblemTypeTag, Grid);
    using LowDimGrid = typename GET_PROP_TYPE(LowDimProblemTypeTag, Grid);
    using BulkElement = typename BulkGrid::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGrid::template Codim<0>::Entity;

    // The Bulk Element Mapper
    using BulkElementMapper = typename GET_PROP_TYPE(BulkProblemTypeTag, ElementMapper);

    static constexpr int bulkDim = BulkGrid::dimension;
    static constexpr int lowDimDim = LowDimGrid::dimension;
    static constexpr int dimWorld = BulkGrid::dimensionworld;
    static constexpr int lowDimDimWorld = LowDimGrid::dimensionworld;
    static_assert(dimWorld == lowDimDimWorld, "dimWorld cannot be different for the two sub-domains!");

    static constexpr int invalidIndex = -5;

    using GlobalPosition = Dune::FieldVector<double, dimWorld>;
    using BulkIndexType = typename BulkGrid::LeafGridView::IndexSet::IndexType;
    using LowDimIndexType = typename LowDimGrid::LeafGridView::IndexSet::IndexType;

public:
    /*!
     * \brief Create the Grid.
     *        Specialization for bulk dim = 2 and lowDim dim = 1
     */
    template<int bd = bulkDim, int ld = lowDimDim>
    static typename std::enable_if<bd == 2 && ld == 1>::type
    makeGrid()
    {
        const std::string fileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Grid, File);

        // set up the grid factories
        auto& bulkFactory = bulkFactory_();
        auto& lowDimFactory = lowDimFactory_();

        // open the mesh file
        std::cout << "Opening " << fileName << std::endl;
        std::ifstream mshFile(fileName.c_str());
        if (mshFile.fail())
            DUNE_THROW(Dune::InvalidStateException, "Could not open the given .msh file. Make sure it exists");
        if (getFileExtension_(fileName) != "msh")
            DUNE_THROW(Dune::InvalidStateException, "Please provide a .msh file for this grid creator");

        // read file until we get to the list of nodes
        std::string line;
        std::getline(mshFile, line);
        while (line.compare(0, 6, "$Nodes") != 0)
            std::getline(mshFile, line);

        // get total number of nodes
        std::getline(mshFile, line);
        const auto numVertices = stringToNumber_(line);

        // read vertices line by line
        std::getline(mshFile, line);
        std::vector<GlobalPosition> bulkVertices(numVertices);
        unsigned int vertexCounter = 0;
        while (line.compare(0, 9, "$EndNodes") != 0)
        {
            // read the coordinates of the vertex
            std::string buf;
            std::stringstream stream(line);

            // drop first entry in line
            stream >> buf;
            std::vector<std::string> coords;
            while (stream >> buf)
                coords.push_back(buf);

            // "create" a vertex
            GlobalPosition v;
            for (int i = 0; i < dimWorld; ++i)
                v[i] = stringToNumber_(coords[i]);

            // insert vertex into container and the bulk grid factory
            bulkVertices[vertexCounter] = v;
            bulkFactory.insertVertex(v);
            std::getline(mshFile, line);
            vertexCounter++;
        }

        // we should always find as many vertices as the mesh file states
        assert(vertexCounter == numVertices);

        // read file until we get to the list of elements
        while(line.compare(0, 9, "$Elements") != 0)
            std::getline(mshFile, line);

        // get total number of elements
        std::getline(mshFile, line);
        const auto numElements = stringToNumber_(line);

        // while inserting the low dim vertices, keep track of indices that have
        // been inserted already and element vertex indices in "original" indexing
        std::vector<GlobalPosition> lowDimVertices;
        std::vector<LowDimIndexType> lowDimVertexIndices;
        std::vector<std::vector<BulkIndexType>> lowDimElemVertices;

        lowDimVertices.reserve(numVertices);
        lowDimVertexIndices.reserve(numVertices);
        lowDimElemVertices.reserve(numElements);

        // the coupling data container lowDim -> bulk
        auto& lowDimCouplingMap = lowDimCouplingElements_();

        // the map to the physical entities
        auto& elemToEntityMap = bulkInsertionIdxToPhysicalEntityMap_();
        elemToEntityMap.resize(numElements);

        // read in the elements line by line
        std::getline(mshFile, line);
        BulkIndexType bulkElementCounter = 0;
        LowDimIndexType lowDimElementCounter = 0;

        while (line.compare(0, 12, "$EndElements") != 0)
        {
            // read the data of this element
            std::string buf;
            std::stringstream stream(line);
            std::vector<unsigned int> elemData;
            while (stream >> buf)
                elemData.push_back(stringToNumber_(buf));

            // get number of gmsh tags for this element
            const auto noOfTags = elemData[2];

            // get element type, make container of indexvertices and insert element
            const auto elemType = elemData[1];

            switch (elemType)
            {
                case 15: // points
                {
                    // do nothing for points
                    break;
                }

                case 1: // 2-node line
                {
                    // the geometry type for lines
                    static const Dune::GeometryType gt = [] () { Dune::GeometryType t; t.makeLine(); return t; } ();

                    // get the connected vertex indices
                    std::vector<LowDimIndexType> elemVertexIndices(2);
                    std::vector<BulkIndexType> originalVertexIndices(2);

                    originalVertexIndices[0] = elemData[3 + noOfTags];
                    originalVertexIndices[1] = elemData[3 + noOfTags + 1];

                    elemVertexIndices[0] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[0]);
                    elemVertexIndices[1] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[1]);

                    for (unsigned int i = 0; i < 2; ++i)
                    {
                        if (elemVertexIndices[i] == invalidIndex)
                        {
                            const auto vIdx = originalVertexIndices[i];
                            const auto v = bulkVertices[vIdx-1];

                            lowDimVertices.push_back(v);
                            lowDimVertexIndices.push_back(vIdx);
                            lowDimFactory.insertVertex(v);

                            elemVertexIndices[i] = lowDimVertices.size()-1;
                        }
                    }

                    // store the low dim elements' vertices
                    lowDimElemVertices.push_back(originalVertexIndices);

                    // insert element into the factory
                    lowDimFactory.insertElement(gt, elemVertexIndices);

                    // increase low dim element counter
                    lowDimElementCounter++;
                    break;
                }

                case 2: // 3-node triangle
                {
                    // the geometry type for triangles
                    static const Dune::GeometryType gt = Dune::GeometryType(Dune::GeometryType::simplex, 2);

                    // we know that the low dim elements have been read already
                    if (lowDimCouplingMap.empty())
                        lowDimCouplingMap.resize(lowDimElementCounter);

                    // get the vertex indices of this bulk element
                    std::vector<BulkIndexType> elemVertexIndices(3);
                    elemVertexIndices[0] = elemData[3 + noOfTags];
                    elemVertexIndices[1] = elemData[3 + noOfTags + 1];
                    elemVertexIndices[2] = elemData[3 + noOfTags + 2];

                    // get the insertion indices of the low dim elements connected to this element
                    const auto lowDimElemIndices = obtainLowDimConnections_(lowDimElemVertices, elemVertexIndices);

                    // if we found some, insert the connections in the map
                    for (auto lowDimElemIdx : lowDimElemIndices)
                        lowDimCouplingMap[lowDimElemIdx].push_back(bulkElementCounter);

                    // subtract 1 from the element indices (gmsh starts at index 1)
                    std::for_each(elemVertexIndices.begin(), elemVertexIndices.end(), [] (auto& i) { i--; });

                    // insert bulk element into the factory
                    bulkFactory.insertElement(gt, elemVertexIndices);

                    // store the physical entity index this elem belongs to
                    elemToEntityMap[bulkElementCounter] = elemData[3];

                    // increase bulk element counter
                    bulkElementCounter++;
                    break;
                }

                case 3: // 4-node quadrilateral
                {
                    // the geometry type for triangles
                    static const Dune::GeometryType gt = Dune::GeometryType(Dune::GeometryType::cube, 2);

                    // we know that the low dim elements have been read already
                    if (lowDimCouplingMap.empty())
                        lowDimCouplingMap.resize(lowDimElementCounter);

                    // get the vertex indices of this bulk element
                    // we rearrange the order here to be dune compatible
                    std::vector<BulkIndexType> elemVertexIndices(4);
                    elemVertexIndices[0] = elemData[3 + noOfTags];
                    elemVertexIndices[1] = elemData[3 + noOfTags + 1];
                    elemVertexIndices[2] = elemData[3 + noOfTags + 3];
                    elemVertexIndices[3] = elemData[3 + noOfTags + 2];

                    // get the insertion indices of the low dim elements connected to this element
                    const auto lowDimElemIndices = obtainLowDimConnections_(lowDimElemVertices, elemVertexIndices);

                    // if we found some, insert the connections in the map
                    for (auto lowDimElemIdx : lowDimElemIndices)
                        lowDimCouplingMap[lowDimElemIdx].push_back(bulkElementCounter);

                    // subtract 1 from the element indices (gmsh starts at index 1)
                    std::for_each(elemVertexIndices.begin(), elemVertexIndices.end(), [] (auto& i) { i--; });

                    // insert bulk element into the factory
                    bulkFactory.insertElement(gt, elemVertexIndices);

                    // store the physical entity index this elem belongs to
                    elemToEntityMap[bulkElementCounter] = elemData[3];

                    // increase bulk element counter
                    bulkElementCounter++;
                    break;
                }

                default:
                    DUNE_THROW(Dune::NotImplemented, "GmshDualFacetGridCreator for given element type");
            }

            // get next line
            std::getline(mshFile, line);
        }

        std::cout << bulkElementCounter << " bulk elements and "
                  << lowDimElementCounter << " low dim elements have been created." << std::endl;


        bulkGridPtr_() = std::shared_ptr<BulkGrid>(bulkFactory.createGrid());
        lowDimGridPtr_() = std::shared_ptr<LowDimGrid>(lowDimFactory.createGrid());

        // set up the map from insertion to actual global indices for bulk elements
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        BulkElementMapper elementMapper(bulkGrid().leafGridView(), Dune::mcmgElementLayout());
#else
        BulkElementMapper elementMapper(bulkGrid().leafGridView());
#endif
        auto& map = bulkInsertionToGlobalIdxMap_();
        map.resize(bulkElementCounter);

        Dune::GeometryType tria; tria.makeTriangle();
        const auto quadOffset = bulkGrid().leafGridView().size(tria);
        for (const auto& element : elements(bulkGrid().leafGridView()))
        {
            const auto insertIdx = element.geometry().type() == tria ?
                                   bulkFactory.insertionIndex(element) :
                                   bulkFactory.insertionIndex(element) + quadOffset;
            map[insertIdx] = elementMapper.index(element);
        }
    }

    /*!
     * \brief Create the Grid.
     *        Specialization for bulk dim = 3 and lowDim dim = 2
     */
    template<int bd = bulkDim, int ld = lowDimDim>
    static typename std::enable_if<bd == 3 && ld == 2>::type
    makeGrid()
    {
        const std::string fileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Grid, File);

        // set up the grid factories
        auto& bulkFactory = bulkFactory_();
        auto& lowDimFactory = lowDimFactory_();

        // open the mesh file
        std::cout << "Opening " << fileName << std::endl;
        std::ifstream mshFile(fileName.c_str());
        if (mshFile.fail())
            DUNE_THROW(Dune::InvalidStateException, "Could not open the given .msh file. Make sure it exists");
        if (getFileExtension_(fileName) != "msh")
            DUNE_THROW(Dune::InvalidStateException, "Please provide a .msh file for this grid creator");

        // read file until we get to the list of nodes
        std::string line;
        std::getline(mshFile, line);
        while (line.compare(0, 6, "$Nodes") != 0)
            std::getline(mshFile, line);

        // get total number of nodes
        std::getline(mshFile, line);
        const auto numVertices = stringToNumber_(line);

        // read vertices line by line
        std::getline(mshFile, line);
        std::vector<GlobalPosition> bulkVertices(numVertices);
        unsigned int vertexCounter = 0;
        while (line.compare(0, 9, "$EndNodes") != 0)
        {
            // read the coordinates of the vertex
            std::string buf;
            std::stringstream stream(line);

            // drop first entry in line
            stream >> buf;
            std::vector<std::string> coords;
            while (stream >> buf)
                coords.push_back(buf);

            // "create" a vertex
            GlobalPosition v;
            for (int i = 0; i < dimWorld; ++i)
                v[i] = stringToNumber_(coords[i]);

            // insert vertex into container and the bulk grid factory
            bulkVertices[vertexCounter] = v;
            bulkFactory.insertVertex(v);
            std::getline(mshFile, line);
            vertexCounter++;
        }

        // we should always find as many vertices as the mesh file states
        assert(vertexCounter == numVertices);

        // read file until we get to the list of elements
        while(line.compare(0, 9, "$Elements") != 0)
            std::getline(mshFile, line);

        // get total number of elements
        std::getline(mshFile, line);
        const auto numElements = stringToNumber_(line);

        // while inserting the low dim vertices, keep track of indices that have
        // been inserted already and element vertex indices in "original" indexing
        std::vector<GlobalPosition> lowDimVertices;
        std::vector<LowDimIndexType> lowDimVertexIndices;
        std::vector<std::vector<BulkIndexType>> lowDimElemVertices;

        lowDimVertices.reserve(numVertices);
        lowDimVertexIndices.reserve(numVertices);
        lowDimElemVertices.reserve(numElements);

        // the coupling data container lowDim -> bulk
        auto& lowDimCouplingMap = lowDimCouplingElements_();

        // read in the elements line by line
        std::getline(mshFile, line);
        BulkIndexType bulkElementCounter = 0;
        LowDimIndexType lowDimElementCounter = 0;
        while (line.compare(0, 12, "$EndElements") != 0)
        {
            // read the data of this element
            std::string buf;
            std::stringstream stream(line);
            std::vector<unsigned int> elemData;
            while (stream >> buf)
                elemData.push_back(stringToNumber_(buf));

            // get number of gmsh tags for this element
            const auto noOfTags = elemData[2];

            // get element type, make container of indexvertices and insert element
            const auto elemType = elemData[1];

            switch (elemType)
            {
                case 15: // points
                {
                    // do nothing for points
                    break;
                }

                case 1: // lines
                {
                    // do nothing for lines
                    break;
                }

                case 2: // 3-node triangle
                {
                    // the geometry type for lines
                    static const Dune::GeometryType gt = Dune::GeometryType(Dune::GeometryType::simplex, 2);

                    // get the connected vertex indices
                    std::vector<LowDimIndexType> elemVertexIndices(3);
                    std::vector<BulkIndexType> originalVertexIndices(3);

                    originalVertexIndices[0] = elemData[3 + noOfTags];
                    originalVertexIndices[1] = elemData[3 + noOfTags + 1];
                    originalVertexIndices[2] = elemData[3 + noOfTags + 2];

                    elemVertexIndices[0] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[0]);
                    elemVertexIndices[1] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[1]);
                    elemVertexIndices[2] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[2]);

                    for (unsigned int i = 0; i < 3; ++i)
                    {
                        if (elemVertexIndices[i] == invalidIndex)
                        {
                            const auto vIdx = originalVertexIndices[i];
                            const auto v = bulkVertices[vIdx-1];

                            lowDimVertices.push_back(v);
                            lowDimVertexIndices.push_back(vIdx);
                            lowDimFactory.insertVertex(v);

                            elemVertexIndices[i] = lowDimVertices.size()-1;
                        }
                    }

                    // store the low dim elements' vertices
                    lowDimElemVertices.push_back(originalVertexIndices);

                    // insert element into the factory
                    lowDimFactory.insertElement(gt, elemVertexIndices);

                    // increase low dim element counter
                    lowDimElementCounter++;
                    break;
                }

                case 3: // 4-node quadrilateral
                {
                    // the geometry type for lines
                    static const Dune::GeometryType gt = Dune::GeometryType(Dune::GeometryType::cube, 2);

                    // get the connected vertex indices
                    std::vector<LowDimIndexType> elemVertexIndices(4);
                    std::vector<BulkIndexType> originalVertexIndices(4);

                    originalVertexIndices[0] = elemData[3 + noOfTags];
                    originalVertexIndices[1] = elemData[3 + noOfTags + 1];
                    originalVertexIndices[2] = elemData[3 + noOfTags + 3];
                    originalVertexIndices[3] = elemData[3 + noOfTags + 2];

                    elemVertexIndices[0] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[0]);
                    elemVertexIndices[1] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[1]);
                    elemVertexIndices[2] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[2]);
                    elemVertexIndices[3] = findIndexInVector_(lowDimVertexIndices, originalVertexIndices[3]);

                    for (unsigned int i = 0; i < 4; ++i)
                    {
                        if (elemVertexIndices[i] == invalidIndex)
                        {
                            const auto vIdx = originalVertexIndices[i];
                            const auto v = bulkVertices[vIdx-1];

                            lowDimVertices.push_back(v);
                            lowDimVertexIndices.push_back(vIdx);
                            lowDimFactory.insertVertex(v);

                            elemVertexIndices[i] = lowDimVertices.size()-1;
                        }
                    }

                    // store the low dim elements' vertices
                    lowDimElemVertices.push_back(originalVertexIndices);

                    // insert element into the factory
                    lowDimFactory.insertElement(gt, elemVertexIndices);

                    // increase low dim element counter
                    lowDimElementCounter++;
                    break;
                }


                case 4: // 4-node tetrahedron
                {
                    // the geometry type for triangles
                    static const Dune::GeometryType gt = Dune::GeometryType(Dune::GeometryType::simplex, 3);

                    // we know that the low dim elements have been read already
                    if (lowDimCouplingMap.empty())
                        lowDimCouplingMap.resize(lowDimElementCounter);

                    // get the vertex indices of this bulk element
                    std::vector<BulkIndexType> elemVertexIndices(4);
                    elemVertexIndices[0] = elemData[3 + noOfTags];
                    elemVertexIndices[1] = elemData[3 + noOfTags + 1];
                    elemVertexIndices[2] = elemData[3 + noOfTags + 2];
                    elemVertexIndices[3] = elemData[3 + noOfTags + 3];

                    // get the insertion indices of the low dim elements connected to this element
                    const auto lowDimElemIndices = obtainLowDimConnections_(lowDimElemVertices, elemVertexIndices);

                    // if we found some, insert the connections in the map
                    for (auto lowDimElemIdx : lowDimElemIndices)
                        lowDimCouplingMap[lowDimElemIdx].push_back(bulkElementCounter);

                    // subtract 1 from the element indices (gmsh starts at index 1)
                    std::for_each(elemVertexIndices.begin(), elemVertexIndices.end(), [] (auto& i) { i--; });

                    // insert bulk element into the factory
                    bulkFactory.insertElement(gt, elemVertexIndices);

                    // increase bulk element counter
                    bulkElementCounter++;
                    break;
                }

                case 5: // 8-node tetrahedron
                {
                    // the geometry type for triangles
                    static const Dune::GeometryType gt = Dune::GeometryType(Dune::GeometryType::cube, 3);

                    // we know that the low dim elements have been read already
                    if (lowDimCouplingMap.empty())
                        lowDimCouplingMap.resize(lowDimElementCounter);

                    // get the vertex indices of this bulk element
                    std::vector<BulkIndexType> elemVertexIndices(8);
                    elemVertexIndices[0] = elemData[3 + noOfTags + 1];
                    elemVertexIndices[1] = elemData[3 + noOfTags];
                    elemVertexIndices[2] = elemData[3 + noOfTags + 5];
                    elemVertexIndices[3] = elemData[3 + noOfTags + 4];
                    elemVertexIndices[4] = elemData[3 + noOfTags + 2];
                    elemVertexIndices[5] = elemData[3 + noOfTags + 3];
                    elemVertexIndices[6] = elemData[3 + noOfTags + 6];
                    elemVertexIndices[7] = elemData[3 + noOfTags + 7];

                    // get the insertion indices of the low dim elements connected to this element
                    const auto lowDimElemIndices = obtainLowDimConnections_(lowDimElemVertices, elemVertexIndices);

                    // if we found some, insert the connections in the map
                    for (auto lowDimElemIdx : lowDimElemIndices)
                        lowDimCouplingMap[lowDimElemIdx].push_back(bulkElementCounter);

                    // subtract 1 from the element indices (gmsh starts at index 1)
                    std::for_each(elemVertexIndices.begin(), elemVertexIndices.end(), [] (auto& i) { i--; });

                    // insert bulk element into the factory
                    bulkFactory.insertElement(gt, elemVertexIndices);

                    // increase bulk element counter
                    bulkElementCounter++;
                    break;
                }

                default:
                    DUNE_THROW(Dune::NotImplemented, "GmshDualFacetGridCreator for given element type");
            }

            // get next line
            std::getline(mshFile, line);
        }

        std::cout << bulkElementCounter << " bulk elements and "
                  << lowDimElementCounter << " low dim elements have been created." << std::endl;

        bulkGridPtr_() = std::shared_ptr<BulkGrid>(bulkFactory.createGrid());
        lowDimGridPtr_() = std::shared_ptr<LowDimGrid>(lowDimFactory.createGrid());

        // set up the map from insertion to actual global indices for bulk elements
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        BulkElementMapper elementMapper(bulkGrid().leafGridView(), Dune::mcmgElementLayout());
#else
        BulkElementMapper elementMapper(bulkGrid().leafGridView());
#endif
        auto& map = bulkInsertionToGlobalIdxMap_();
        map.resize(bulkElementCounter);
        for (const auto& element : elements(bulkGrid().leafGridView()))
            map[bulkFactory.insertionIndex(element)] = elementMapper.index(element);
    }

     /*!
     * \brief Returns a reference to the bulk grid.
     */
    static BulkGrid& bulkGrid()
    { return *bulkGridPtr_(); }

    /*!
     * \brief Returns a reference to the lower dimensional facet grid.
     */
    static LowDimGrid& lowDimGrid()
    { return *lowDimGridPtr_(); }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    static void loadBalance()
    {
        bulkGrid().loadBalance();
        lowDimGrid().loadBalance();
    }

    static std::vector<BulkIndexType> getCoupledBulkElementIndices(const LowDimElement& lowDimElement)
    {
        const auto& insertionIndices = lowDimCouplingElements_()[lowDimFactory_().insertionIndex(lowDimElement)];
        const auto& bulkInsertionToGlobalMap = bulkInsertionToGlobalIdxMap_();

        const auto n = insertionIndices.size();
        std::vector<BulkIndexType> globalIndices(n);
        for (unsigned int i = 0; i < n; ++i)
            globalIndices[i] = bulkInsertionToGlobalMap[insertionIndices[i]];

        return globalIndices;
    }

    static BulkIndexType getInsertionIndex(const BulkElement& element)
    { return bulkFactory_().insertionIndex(element); }

    static LowDimIndexType getInsertionIndex(const LowDimElement& element)
    { return lowDimFactory_().insertionIndex(element); }

    static unsigned int getPhysicalEntityIndex(const BulkElement& element)
    { return bulkInsertionIdxToPhysicalEntityMap_()[getInsertionIndex(element)]; }

private:

    /*!
     * \brief Returns a reference to the vector containing lists of bulk element indices
     *        that are connected to the low dim elements. Note that the access occurs
     *        with the insertion index of a low dim element.
     */
    static std::vector<std::vector<BulkIndexType>>& lowDimCouplingElements_()
    {
        static std::vector<std::vector<BulkIndexType>> couplingMap_;
        return couplingMap_;
    }

    /*!
     * \brief Returns a reference to the map mapping bulk element insertion indices
     *        to actual global indices
     */
    static std::vector<BulkIndexType>& bulkInsertionToGlobalIdxMap_()
    {
        static std::vector<BulkIndexType> insertionIdxMap_;
        return insertionIdxMap_;
    }

    /*!
     * \brief Returns a reference to the map mapping bulk element insertion indices
     *        to the physical entity the element corresponds to
     */
    static std::vector<unsigned int>& bulkInsertionIdxToPhysicalEntityMap_()
    {
        static std::vector<unsigned int> physicalEntityMap_;
        return physicalEntityMap_;
    }

    /*!
     * \brief Returns a reference to the bulk grid pointer (std::shared_ptr<BulkGrid>)
     */
    static std::shared_ptr<BulkGrid>& bulkGridPtr_()
    {
        static std::shared_ptr<BulkGrid> bulkGridPtr_;
        return bulkGridPtr_;
    }

    /*!
     * \brief Returns a reference to the lowDim grid pointer (std::shared_ptr<LowDimGrid>)
     */
    static std::shared_ptr<LowDimGrid>& lowDimGridPtr_()
    {
        static std::shared_ptr<LowDimGrid> lowDimGridPtr_;
        return lowDimGridPtr_;
    }

    static Dune::GridFactory<BulkGrid>& bulkFactory_()
    {
        static Dune::GridFactory<BulkGrid> bulkFactory_;
        return bulkFactory_;
    }

    static Dune::GridFactory<LowDimGrid>& lowDimFactory_()
    {
        static Dune::GridFactory<LowDimGrid> lowDimFactory_;
        return lowDimFactory_;
    }

    /*!
     * \brief Returns a list of low dim elements that are connected to this bulk element
     */
    static std::vector<LowDimIndexType> obtainLowDimConnections_(const std::vector<std::vector<LowDimIndexType>>& lowDimElemData,
                                                                 const std::vector<BulkIndexType>& bulkElemVertices)
    {
        std::vector<LowDimIndexType> lowDimElemIndices = std::vector<LowDimIndexType>();

        // check for each low dim element if it is contained
        unsigned int lowDimElementIdx = 0;
        for (const auto& lowDimElem : lowDimElemData)
        {
            const bool isContained = [&lowDimElem, &bulkElemVertices] ()
            {
                for (auto idx : lowDimElem)
                    if (findIndexInVector_(bulkElemVertices, idx) == invalidIndex)
                        return false;
                return true;
            } ();

            // insert in list if low dim element is in bulk element
            if (isContained)
                lowDimElemIndices.push_back(lowDimElementIdx);

            // keep track of the low dim element insertion index
            lowDimElementIdx++;
        }

        return lowDimElemIndices;
    }

    /*!
     * \brief Returns the local index of a given value in a given vector
     */
    template<typename IdxType1, typename IdxType2>
    static int findIndexInVector_(const std::vector<IdxType1>& vector, const IdxType2 globalIdx)
    {
        auto it = std::find(vector.begin(), vector.end(), globalIdx);

        // we use the -5 here to detect whether or not an index exists already
        if (it != vector.end())
            return std::distance(vector.begin(), it);
        else
            return invalidIndex;
    }

    /*!
     * \brief Returns the filename extension of a given filename
     */
    static std::string getFileExtension_(const std::string& fileName)
    {
        std::size_t i = fileName.rfind('.', fileName.length());
        if (i != std::string::npos)
            return(fileName.substr(i+1, fileName.length() - i));
        else
            DUNE_THROW(Dune::IOError, "Please provide and extension for your grid file ('"<< fileName << "')!");
    }

    /*!
     * \brief Turns a number in form of a string into a Scalar
     */
    static Scalar stringToNumber_(const std::string& numberAsString)
    {
        Scalar value;

        std::stringstream stream(numberAsString);
        stream >> value;

        if (stream.fail())
        {
            std::runtime_error e(numberAsString);
            throw e;
        }

        return value;
   }
};

} // end namespace

#endif // DUMUX_GMSHDUALFACETGRIDCREATOR_HH
