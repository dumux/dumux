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
    template<class T = TypeTag>
    static typename std::enable_if<bulkDim == 2 && lowDimDim == 1>::type
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

                case 1: // 2-node line
                {
                    // get the connected vertex indices
                    std::vector<LowDimIndexType> elemVertexIndices(2);

                    const auto vIdx1 = elemData[3 + noOfTags];
                    const auto vIdx2 = elemData[3 + noOfTags + 1];

                    elemVertexIndices[0] = findIndexInVector_(lowDimVertexIndices, vIdx1);
                    elemVertexIndices[1] = findIndexInVector_(lowDimVertexIndices, vIdx2);

                    if (elemVertexIndices[0] == invalidIndex)
                    {
                        const auto v1 = bulkVertices[vIdx1-1];

                        lowDimVertices.push_back(v1);
                        lowDimVertexIndices.push_back(vIdx1);
                        lowDimFactory.insertVertex(v1);

                        elemVertexIndices[0] = lowDimVertices.size()-1;
                    }
                    if (elemVertexIndices[1] == invalidIndex)
                    {
                        const auto v2 = bulkVertices[vIdx2-1];

                        lowDimVertices.push_back(v2);
                        lowDimVertexIndices.push_back(vIdx2);
                        lowDimFactory.insertVertex(v2);

                        elemVertexIndices[1] = lowDimVertices.size()-1;
                    }

                    // store the low dim elements' vertices
                    lowDimElemVertices.push_back(std::vector<BulkIndexType>({vIdx1, vIdx2}));

                    // insert element into the factory
                    lowDimFactory.insertElement(Dune::GeometryTypes::line, elemVertexIndices);

                    // increase low dim element counter
                    lowDimElementCounter++;
                    break;
                }

                case 2: // 3-node triangle
                {
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
                    elemVertexIndices[0] -= 1;
                    elemVertexIndices[1] -= 1;
                    elemVertexIndices[2] -= 1;

                    // insert bulk element into the factory
                    bulkFactory.insertElement(Dune::GeometryTypes::triangle, elemVertexIndices);

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
        BulkElementMapper elementMapper(bulkGrid().leafGridView(), Dune::mcmgElementLayout());
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
