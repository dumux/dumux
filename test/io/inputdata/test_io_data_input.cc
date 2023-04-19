// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test and example of how to read custom user data in Dumux
 */
#include <config.h>

#include <algorithm>
#include <vector>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertreeparser.hh>

#include <dumux/io/container.hh>
#include <dumux/io/xml/tinyxml2.h>

////////////////////////
// the main function
////////////////////////
int main()
{
    using namespace Dumux;

    {
        std::cout << "Reading a simple list of numbers into a vector\n";
        std::cout << "-- Read numbers.txt: \n";

        // the file extension ".txt" is arbitrary, works with any other extension or none
        const auto numbers = readFileToContainer<std::vector<double>>("numbers.txt");

        std::copy(numbers.begin(), numbers.end(), std::ostream_iterator<double>(std::cout, ", "));
        std::cout << "\n" << std::endl;
    }

    {
        std::cout << "Reading a three-column list of numbers into a vector\n";
        std::cout << "-- Read coordinates.txt: \n";

        // the file extension ".txt" is arbitrary, works with any other extension or none
        const auto coordinates = readFileToContainer<std::vector<Dune::FieldVector<double, 3>>>("coordinates.txt");

        std::copy(coordinates.begin(), coordinates.end(), std::ostream_iterator<Dune::FieldVector<double, 3>>(std::cout, ", "));
        std::cout << "\n" << std::endl;
    }

    {
        std::cout << "Reading a key-value ini-style file into a Dune::ParameterTree\n";
        std::cout << "-- Read config.ini: \n";

        Dune::ParameterTree config;
        // the file extension ".ini" is arbitrary, works with any other extension or none
        Dune::ParameterTreeParser::readINITree("config.ini", config);

        config.report();
        std::cout << "\n";

        std::cout << "-- Parsing MyData.InjectionRate into a vector: \n";

        const auto injectionRates = config.get<std::vector<double>>("MyData.InjectionRate");

        std::copy(injectionRates.begin(), injectionRates.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n" << std::endl;
    }

    {
        std::cout << "Reading a XML-formatted data using tiny xml\n";

        tinyxml2::XMLDocument xmlData;
        // the file extension ".xml" is arbitrary, works with any other extension or none
        const auto returnCode = xmlData.LoadFile("mydata.xml");
        if (returnCode != tinyxml2::XML_SUCCESS)
            DUNE_THROW(Dune::IOError, "Couldn't open XML file.");
        const tinyxml2::XMLElement* inputData = xmlData.FirstChildElement("InputData");

        std::cout << "-- Read InjectionRates field: \n";

        std::stringstream injectionData(inputData->FirstChildElement("InjectionRates")->GetText());
        const auto injectionRates = readStreamToContainer<std::vector<double>>(injectionData);

        std::copy(injectionRates.begin(), injectionRates.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n";

        std::cout << "-- Read BoundaryTypes field: \n";

        std::stringstream boundaryData(inputData->FirstChildElement("BoundaryTypes")->GetText());
        const auto boundaryTypes = readStreamToContainer<std::vector<int>>(boundaryData);

        std::copy(boundaryTypes.begin(), boundaryTypes.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << "\n" << std::endl;
    }

    return 0;
}
