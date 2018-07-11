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

#include <dune/common/exceptions.hh>
#include <dumux/io/xml/tinyxml2.h>

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
    Container readData(const std::string& name, const DataType& type)
    {
        using namespace tinyxml2;
        XMLElement *pieceNode = doc_.FirstChildElement("VTKFile")->FirstChildElement("UnstructuredGrid")->FirstChildElement("Piece");
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get Piece node in " << fileName_ << ".");

        XMLElement *dataNode = nullptr;
        if (type == DataType::pointData)
            dataNode = pieceNode->FirstChildElement("PointData");
        else if (type == DataType::cellData)
            dataNode = pieceNode->FirstChildElement("CellData");
        else
            DUNE_THROW(Dune::IOError, "Only cell and point data are supported.");

        if (dataNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get data node in " << fileName_ << ".");

        // loop over XML node siblings to find the correct data array
        XMLElement *dataArray = dataNode->FirstChildElement("DataArray");
        for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
        {
            const char *attributeText = dataArray->Attribute("Name");

            if (attributeText == nullptr)
                DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a data array.");

            if (attributeText == name)
                break;
        }
        if (dataArray == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't find the data array " << name << ".");

        Container data;
        std::stringstream dataStream(dataArray->GetText());
        std::istream_iterator<typename Container::value_type> it(dataStream);
        std::copy(it, std::istream_iterator<typename Container::value_type>(), std::back_inserter(data));
        return data;
    }

private:
    const std::string fileName_; //!< the vtu file name
    tinyxml2::XMLDocument doc_; //!< the xml document created from file with name fileName_
};

} // end namespace Dumux

#endif
