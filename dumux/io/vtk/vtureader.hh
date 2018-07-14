// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \ingroup InputOutput
 * \brief A vtu reader using tinyxml2 as xml backend
 */
#ifndef DUMUX_IO_VTK_VTUREADER_HH
#define DUMUX_IO_VTK_VTUREADER_HH

#include <iostream>
#include <iterator>
#include <algorithm>
#include <memory>
#include <unordered_map>

#include <dune/common/exceptions.hh>
#include <dumux/io/xml/tinyxml2.h>
#include <dune/grid/common/gridfactory.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A vtu reader using tinyxml2 as xml backend
 */
class VTUReader
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
    VTUReader(const std::string& fileName)
    : fileName_(fileName)
    {
        using namespace tinyxml2;
        const auto eResult = doc_.LoadFile(fileName.c_str());
        if (eResult != XML_SUCCESS)
            DUNE_THROW(Dune::IOError, "Couldn't open XML file " << fileName << ".");
    }

    /*!
     * \brief The contructor creates a tinyxml2::XMLDocument from file
     * \tparam Container a container type that has begin(), end(), push_back(), e.g. std::vector<>
     * \param name the name attribute of the data array to read
     * \param type the data array type
     */
    template<class Container>
    Container readData(const std::string& name, const DataType& type) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'Piece' node in " << fileName_ << ".");

        const XMLElement* dataNode = getDataNode_(pieceNode, type);
        if (dataNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'PointData' or 'CellData' node in " << fileName_ << ".");

        const XMLElement* dataArray = findDataArray_(dataNode, name);
        if (dataArray == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't find the data array " << name << ".");

        return parseDataArray_<Container>(dataArray);
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, without reading cell and point data
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(bool verbose = false) const
    {
        if (verbose) std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;

        // make a grid factory
        Dune::GridFactory<Grid> factory;

        readGrid_(factory, verbose);

        return std::unique_ptr<Grid>(factory.createGrid());
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, reading cell and point data
     * \note the factory will be needed outside to interpret the data via the factory's insertion indices
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(Dune::GridFactory<Grid>& factory, Data& cellData, Data& pointData, bool verbose = false) const
    {
        if (verbose) std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;

        readGrid_(factory, verbose);
        readGridData_(cellData, pointData, verbose);

        return std::unique_ptr<Grid>(factory.createGrid());
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, reading cell and point data
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(Data& cellData, Data& pointData, bool verbose = false) const
    {
        if (verbose) std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;

        // make a grid factory
        Dune::GridFactory<Grid> factory;

        readGrid_(factory, verbose);
        readGridData_(cellData, pointData, verbose);

        return std::unique_ptr<Grid>(factory.createGrid());
    }

private:
    template<class Grid>
    void readGrid_(Dune::GridFactory<Grid>& factory, bool verbose = false) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'Piece' node in " << fileName_ << ".");

        const XMLElement* pointsNode = pieceNode->FirstChildElement("Points")->FirstChildElement("DataArray");
        if (pointsNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get data array of points in " << fileName_ << ".");

        using Point = Dune::FieldVector<double, 3>;
        std::vector<Point> points;
        std::stringstream dataStream(pointsNode->GetText());
        std::istream_iterator<Point> it(dataStream);
        std::copy(it, std::istream_iterator<Point>(), std::back_inserter(points));

        if (Grid::dimensionworld < 3)
            DUNE_THROW(Dune::NotImplemented, "VTKReader for dimworld < 3");

        if (verbose) std::cout << "Found " << points.size() << " vertices." << std::endl;

        // insert vertices to the grid factory
        for (auto&& point : points)
            factory.insertVertex(std::move(point));

        const XMLElement* cellsNode = pieceNode->FirstChildElement("Cells");
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
            unsigned int offset = offsets[i];
            std::vector<unsigned int> corners; corners.reserve(offset-lastOffset);
            for (unsigned int j = lastOffset; j < offset; ++j)
                corners.emplace_back(connectivity[j]);
            insertElement_(factory, std::move(corners), types[i]);
            lastOffset = offset;
        }
    }

    void readGridData_(Data& cellData, Data& pointData, bool verbose = false) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'Piece' node in " << fileName_ << ".");

        const XMLElement* cellDataNode = getDataNode_(pieceNode, DataType::cellData);
        if (cellDataNode != nullptr)
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

    const tinyxml2::XMLElement* getPieceNode_() const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = doc_.FirstChildElement("VTKFile");
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'VTKFile' node in " << fileName_ << ".");

        pieceNode = pieceNode->FirstChildElement("UnstructuredGrid");
        if (pieceNode == nullptr)
            pieceNode = doc_.FirstChildElement("VTKFile")->FirstChildElement("PolyData");
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'UnstructuredGrid' or 'PolyData' node in " << fileName_ << ".");

        return pieceNode->FirstChildElement("Piece");
    }

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

    template<class Container>
    Container parseDataArray_(const tinyxml2::XMLElement* dataArray) const
    {
        Container data;
        std::stringstream dataStream(dataArray->GetText());
        std::istream_iterator<typename Container::value_type> it(dataStream);
        std::copy(it, std::istream_iterator<typename Container::value_type>(), std::back_inserter(data));
        return data;
    }

    template<class Grid>
    void insertElement_(Dune::GridFactory<Grid>& factory, std::vector<unsigned int>&& corners, unsigned int vtkCellType) const
    {
        switch (vtkCellType)
        {
            case 12:
                factory.insertElement(Dune::GeometryTypes::hexahedron, corners);
                break;
            default:
                DUNE_THROW(Dune::NotImplemented, "VTK cell type " << vtkCellType);
        }
    }

    const std::string fileName_; //!< the vtu file name
    tinyxml2::XMLDocument doc_; //!< the xml document created from file with name fileName_
};

} // end namespace Dumux

#endif
