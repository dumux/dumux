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
 * \ingroup InputOutput
 * \brief A vtk file reader using tinyxml2 as xml backend
 */
#ifndef DUMUX_IO_VTK_VTKREADER_HH
#define DUMUX_IO_VTK_VTKREADER_HH

#include <iostream>
#include <iterator>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/io/container.hh>
#include <dumux/io/xml/tinyxml2.h>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A vtk file reader using tinyxml2 as xml backend
 */
class VTKReader
{
public:
    /*!
     * \brief The data array types
     */
    enum class DataType { cellData, pointData };

    //! the cell / point data type for point data read from a grid file
    using Data = std::unordered_map<std::string, std::vector<double>>;

    /*!
     * \brief The contructor creates a tinyxml2::XMLDocument from file
     */
    explicit VTKReader(const std::string& fileName)
    {
        using namespace tinyxml2;
        fileName_ = Dune::MPIHelper::getCollectiveCommunication().size() > 1 ?
                        getProcessFileName_(fileName) : fileName;

        const auto eResult = doc_.LoadFile(fileName_.c_str());
        if (eResult != tinyxml2::XML_SUCCESS)
            DUNE_THROW(Dune::IOError, "Couldn't open XML file " << fileName_ << ".");

        const XMLElement* pieceNode = getPieceNode_();
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'Piece' node in " << fileName_ << ".");
    }

    /*!
     * \brief Reviews data from the vtk file to check if there is a data array with a specified name
     * \param name the name attribute of the data array to read
     * \param type the data array type
     */
    bool hasData(const std::string& name, const DataType& type) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        const XMLElement* dataNode = getDataNode_(pieceNode, type);
        if (dataNode == nullptr)
            return false;

        const XMLElement* dataArray = findDataArray_(dataNode, name);
        if (dataArray == nullptr)
            return false;

        return true;
    }

    /*!
     * \brief read data from the vtk file to a container, e.g. std::vector<double>
     * \tparam Container a container type that has begin(), end(), push_back(), e.g. std::vector<>
     * \param name the name attribute of the data array to read
     * \param type the data array type
     */
    template<class Container>
    Container readData(const std::string& name, const DataType& type) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        const XMLElement* dataNode = getDataNode_(pieceNode, type);
        if (dataNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'PointData' or 'CellData' node in " << fileName_ << ".");

        const XMLElement* dataArray = findDataArray_(dataNode, name);
        if (dataArray == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't find the data array " << name << ".");

        return parseDataArray_<Container>(dataArray);
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, ignoring cell and point data
     * \param verbose if the output should be verbose
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(bool verbose = false) const
    {
        static_assert(!Dune::Capabilities::isCartesian<Grid>::v, "Grid reader only supports unstructured grid implementations");

        if (verbose) std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;

        // make a grid factory
        Dune::GridFactory<Grid> factory;

        readGrid_(factory, verbose);

        return std::unique_ptr<Grid>(factory.createGrid());
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, ignoring cell and point data
     * \note use this signature if the factory might be needed outside to interpret the data via the factory's insertion indices
     * \param verbose if the output should be verbose
     * \param factory the (emtpy) grid factory
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(Dune::GridFactory<Grid>& factory, bool verbose = false) const
    {
        static_assert(!Dune::Capabilities::isCartesian<Grid>::v, "Grid reader only supports unstructured grid implementations");

        if (verbose) std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;

        readGrid_(factory, verbose);

        return std::unique_ptr<Grid>(factory.createGrid());
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, reading all cell and point data
     * \note the factory will be needed outside to interpret the data via the factory's insertion indices
     * \param factory the (emtpy) grid factory
     * \param cellData the cell data arrays to be filled
     * \param pointData the point data arrays to be filled
     * \param verbose if the output should be verbose
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(Dune::GridFactory<Grid>& factory, Data& cellData, Data& pointData, bool verbose = false) const
    {
        static_assert(!Dune::Capabilities::isCartesian<Grid>::v, "Grid reader only supports unstructured grid implementations");

        if (verbose) std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;

        readGrid_(factory, verbose);
        readGridData_(cellData, pointData, verbose);

        return std::unique_ptr<Grid>(factory.createGrid());
    }

private:
    /*!
     * \brief get the vtk filename for the current processor
     */
    std::string getProcessFileName_(const std::string& pvtkFileName)
    {
        using namespace tinyxml2;

        XMLDocument pDoc;
        const auto eResult = pDoc.LoadFile(pvtkFileName.c_str());
        if (eResult != XML_SUCCESS)
            DUNE_THROW(Dune::IOError, "Couldn't open XML file " << pvtkFileName << ".");

        // get the first piece node
        const XMLElement* pieceNode = getPieceNode_(pDoc, pvtkFileName);
        const auto myrank = Dune::MPIHelper::getCollectiveCommunication().rank();
        for (int rank = 0; rank < myrank; ++rank)
        {
            pieceNode = pieceNode->NextSiblingElement("Piece");
            if (pieceNode == nullptr)
                DUNE_THROW(Dune::IOError, "Couldn't find 'Piece' node for rank "
                                          << rank << " in " << pvtkFileName << ".");
        }

        const char *vtkFileName = pieceNode->Attribute("Source");
        if (vtkFileName == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'Source' attribute of 'Piece' node no. " << myrank << " in " << pvtkFileName);

        return vtkFileName;
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file
     * \param factory the (emtpy) grid factory
     * \param verbose if the output should be verbose
     */
    template<class Grid>
    void readGrid_(Dune::GridFactory<Grid>& factory, bool verbose = false) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        const XMLElement* pointsNode = pieceNode->FirstChildElement("Points")->FirstChildElement("DataArray");
        if (pointsNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get data array of points in " << fileName_ << ".");

        using Point3D = Dune::FieldVector<double, 3>;
        std::vector<Point3D> points3D;
        std::stringstream dataStream(pointsNode->GetText());
        std::istream_iterator<Point3D> it(dataStream);
        std::copy(it, std::istream_iterator<Point3D>(), std::back_inserter(points3D));

        // adapt point dimensions if grid dimension is smaller than 3
        auto points = adaptPointDimension_<Grid::dimensionworld>(std::move(points3D));

        if (verbose) std::cout << "Found " << points.size() << " vertices." << std::endl;

        // insert vertices to the grid factory
        for (auto&& point : points)
            factory.insertVertex(std::move(point));

        const XMLElement* cellsNode = pieceNode->FirstChildElement("Cells");
        const XMLElement* linesNode = pieceNode->FirstChildElement("Lines");
        if (cellsNode)
        {
            const XMLElement* connectivityNode = findDataArray_(cellsNode, "connectivity");
            const XMLElement* offsetsNode = findDataArray_(cellsNode, "offsets");
            const XMLElement* typesNode = findDataArray_(cellsNode, "types");

            const auto connectivity = parseDataArray_<std::vector<unsigned int>>(connectivityNode);
            const auto offsets = parseDataArray_<std::vector<unsigned int>>(offsetsNode);
            const auto types = parseDataArray_<std::vector<unsigned int>>(typesNode);

            if (verbose) std::cout << "Found " << offsets.size() << " element." << std::endl;

            unsigned int lastOffset = 0;
            for (unsigned int i = 0; i < offsets.size(); ++i)
            {
                const auto geomType = vtkToDuneGeomType_(types[i]);
                unsigned int offset = offsets[i];
                std::vector<unsigned int> corners; corners.resize(offset-lastOffset);
                for (unsigned int j = 0; j < offset-lastOffset; ++j)
                    corners[Dune::VTK::renumber(geomType, j)] = connectivity[lastOffset+j];
                factory.insertElement(geomType, std::move(corners));
                lastOffset = offset;
            }
        }
        // for poly data
        else if (linesNode)
        {
            // sanity check
            if (Grid::dimension != 1)
                DUNE_THROW(Dune::IOError, "Grid expects dimension " << Grid::dimension
                                           << " but " << fileName_ << " contains a 1D grid.");

            const XMLElement* connectivityNode = findDataArray_(linesNode, "connectivity");
            const XMLElement* offsetsNode = findDataArray_(linesNode, "offsets");

            const auto connectivity = parseDataArray_<std::vector<unsigned int>>(connectivityNode);
            const auto offsets = parseDataArray_<std::vector<unsigned int>>(offsetsNode);

            if (verbose) std::cout << "Found " << offsets.size() << " polylines." << std::endl;

            unsigned int lastOffset = 0;
            for (unsigned int i = 0; i < offsets.size(); ++i)
            {
                // a polyline can have many points in the VTK format
                // split the line in segments with two points
                unsigned int offset = offsets[i];
                for (unsigned int j = 0; j < offset-lastOffset-1; ++j)
                    factory.insertElement(Dune::GeometryTypes::line,
                        std::vector<unsigned int>({ connectivity[lastOffset+j], connectivity[lastOffset+j+1] }));
                lastOffset = offset;
            }
        }
        else
            DUNE_THROW(Dune::IOError, "No Cells or Lines element found in " << fileName_);
    }

    /*!
     * \brief Read a grid data from a vtk/vtu/vtp file
     * \param cellData the cell data arrays to be filled
     * \param pointData the point data arrays to be filled
     * \param verbose if the output should be verbose
     */
    void readGridData_(Data& cellData, Data& pointData, bool verbose = false) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        const XMLElement* cellDataNode = getDataNode_(pieceNode, DataType::cellData);
        if (cellDataNode != nullptr)
        {
            const XMLElement* cellsNode = pieceNode->FirstChildElement("Cells");
            const XMLElement* linesNode = pieceNode->FirstChildElement("Lines");

            if (cellsNode)
            {
                const XMLElement* dataArray = cellDataNode->FirstChildElement("DataArray");
                for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
                {
                    const char *attributeText = dataArray->Attribute("Name");

                    if (attributeText == nullptr)
                        DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a cell data array.");

                    cellData[std::string(attributeText)] = parseDataArray_<std::vector<double>>(dataArray);
                }
            }
            // for poly data
            else if (linesNode)
            {
                // first parse all the cell data (each cell in this sense can be a polyline)
                const XMLElement* dataArray = cellDataNode->FirstChildElement("DataArray");
                if (dataArray)
                {
                    Data polyLineCellData;
                    for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
                    {
                        const char *attributeText = dataArray->Attribute("Name");

                        if (attributeText == nullptr)
                            DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a cell data array.");

                        polyLineCellData[std::string(attributeText)] = parseDataArray_<std::vector<double>>(dataArray);
                    }

                    // a polyline can have many points in the VTK format
                    // we split the line in segments with two points
                    // so we also need to duplicate the cell data to fit the increased line number
                    const XMLElement* offsetsNode = findDataArray_(linesNode, "offsets");
                    const auto offsets = parseDataArray_<std::vector<unsigned int>>(offsetsNode);

                    if (offsets.size() != polyLineCellData.begin()->second.size())
                        DUNE_THROW(Dune::IOError, "Expected the same number of cell data entries (is "
                                                   << polyLineCellData.begin()->second.size()
                                                   << ") as polylines (" << offsets.size() << ")!");

                    // count the number of Dune cells to be able to resize the data vectors
                    unsigned int lastOffset = 0;
                    std::size_t numCells = 0;
                    for (unsigned int i = 0; i < offsets.size(); ++i)
                    {
                        unsigned int offset = offsets[i];
                        for (unsigned int j = 0; j < offset-lastOffset-1; ++j)
                            ++numCells;
                        lastOffset = offset;
                    }

                    // create the data arrays
                    for (const auto& dArray : polyLineCellData)
                    {
                        cellData[dArray.first] = std::vector<double>(numCells);
                        auto& cd = cellData[dArray.first];
                        const auto& pd = dArray.second;

                        lastOffset = 0;
                        std::size_t cellIdx = 0;
                        for (unsigned int i = 0; i < offsets.size(); ++i)
                        {
                            unsigned int offset = offsets[i];
                            for (unsigned int j = 0; j < offset-lastOffset-1; ++j)
                                cd[cellIdx++] = pd[i];
                            lastOffset = offset;
                        }
                    }
                }
            }
            else
                DUNE_THROW(Dune::IOError, "No Cells or Lines element found in " << fileName_);
        }

        const XMLElement* pointDataNode = getDataNode_(pieceNode, DataType::pointData);
        if (pointDataNode != nullptr)
        {
            const XMLElement* dataArray = pointDataNode->FirstChildElement("DataArray");
            for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
            {
                const char *attributeText = dataArray->Attribute("Name");

                if (attributeText == nullptr)
                    DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a point data array.");

                pointData[std::string(attributeText)] = parseDataArray_<std::vector<double>>(dataArray);
            }
        }
    }

    /*!
     * \brief Get the piece node of the XMLDocument
     * \note Returns nullptr if the piece node wasn't found
     */
    const tinyxml2::XMLElement* getPieceNode_() const
    { return getPieceNode_(doc_, fileName_); }

    /*!
     * \brief Get the piece node an xml document
     * \note Returns nullptr if the piece node wasn't found
     * \param doc an xml document
     * \param fileName a file name the doc was created from
     */
    const tinyxml2::XMLElement* getPieceNode_(const tinyxml2::XMLDocument& doc, const std::string& fileName) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = doc.FirstChildElement("VTKFile");
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'VTKFile' node in " << fileName << ".");

        pieceNode = pieceNode->FirstChildElement("UnstructuredGrid");
        if (pieceNode == nullptr)
            pieceNode = doc.FirstChildElement("VTKFile")->FirstChildElement("PolyData");
        if (pieceNode == nullptr)
            pieceNode = doc.FirstChildElement("VTKFile")->FirstChildElement("PUnstructuredGrid");
        if (pieceNode == nullptr)
            pieceNode = doc.FirstChildElement("VTKFile")->FirstChildElement("PPolyData");
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'UnstructuredGrid', 'PUnstructuredGrid', 'PolyData', or 'PPolyData' node in " << fileName << ".");

        return pieceNode->FirstChildElement("Piece");
    }

    /*!
     * \brief Get the piece node of the XMLDocument
     * \param pieceNode the pieceNode of the vtk file
     * \param type the vtk data type (cell data or point data)
     * \note Returns nullptr if the data node wasn't found
     */
    const tinyxml2::XMLElement* getDataNode_(const tinyxml2::XMLElement* pieceNode, const DataType& type) const
    {
        using namespace tinyxml2;

        const XMLElement* dataNode = nullptr;
        if (type == DataType::pointData)
            dataNode = pieceNode->FirstChildElement("PointData");
        else if (type == DataType::cellData)
            dataNode = pieceNode->FirstChildElement("CellData");
        else
            DUNE_THROW(Dune::IOError, "Only cell and point data are supported.");

        return dataNode;
    }

    /*!
     * \brief Find a data array with a specific name
     * \param dataNode a cell or point data node
     * \param name the name of the data array to be found
     * \note Returns nullptr if the data array wasn't found
     */
    const tinyxml2::XMLElement* findDataArray_(const tinyxml2::XMLElement* dataNode, const std::string& name) const
    {
        using namespace tinyxml2;

        // loop over XML node siblings to find the correct data array
        const XMLElement* dataArray = dataNode->FirstChildElement("DataArray");
        for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
        {
            const char *attributeText = dataArray->Attribute("Name");

            if (attributeText == nullptr)
                DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a data array.");

            if (attributeText == name)
                break;
        }

        return dataArray;
    }

    /*!
     * \brief Parses the text of a data array into a container
     * \tparam Container a container type that has begin(), end(), push_back(), e.g. std::vector<double>
     * \param dataArray the data array node to be parsed
     */
    template<class Container>
    Container parseDataArray_(const tinyxml2::XMLElement* dataArray) const
    {
        std::stringstream dataStream(dataArray->GetText());
        return readStreamToContainer<Container>(dataStream);
    }

    /*!
     * \brief Return the Dune::GeometryType for a given VTK geometry type
     * \param vtkCellType the vtk cell type
     */
    Dune::GeometryType vtkToDuneGeomType_(unsigned int vtkCellType) const
    {
        switch (vtkCellType)
        {
            case Dune::VTK::GeometryType::vertex: return Dune::GeometryTypes::vertex;
            case Dune::VTK::GeometryType::line: return Dune::GeometryTypes::line;
            case Dune::VTK::GeometryType::triangle: return Dune::GeometryTypes::triangle;
            case Dune::VTK::GeometryType::quadrilateral: return Dune::GeometryTypes::quadrilateral;
            case Dune::VTK::GeometryType::hexahedron: return Dune::GeometryTypes::hexahedron;
            case Dune::VTK::GeometryType::prism: return Dune::GeometryTypes::prism;
            case Dune::VTK::GeometryType::pyramid: return Dune::GeometryTypes::pyramid;
            default: DUNE_THROW(Dune::NotImplemented, "VTK cell type " << vtkCellType);
        }
    }

    template<int dim, std::enable_if_t<(dim < 3), int> = 0>
    std::vector<Dune::FieldVector<double, dim>>
    adaptPointDimension_(std::vector<Dune::FieldVector<double, 3>>&& points3D) const
    {
        std::vector<Dune::FieldVector<double, dim>> points(points3D.size());
        for (std::size_t i = 0; i < points.size(); ++i)
            for (int j = 0; j < dim; ++j)
                points[i][j] = points3D[i][j];
        return points;
    }

    template<int dim, std::enable_if_t<(dim == 3), int> = 0>
    std::vector<Dune::FieldVector<double, dim>>
    adaptPointDimension_(std::vector<Dune::FieldVector<double, 3>>&& points3D) const
    { return std::move(points3D); }

    std::string fileName_; //!< the vtk file name
    tinyxml2::XMLDocument doc_; //!< the xml document created from file with name fileName_
};

} // end namespace Dumux

#endif
