// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the pore-network grid manager with pores and throats given in memory
 *
 * Every vertex and element carries the parameters of the pore and throat at its position, as a grid read from a DGF
 * file of the same network does, also after grid sanitation removed a cluster not connected to the inlet.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/grid/porenetwork/networkdata.hh>

namespace {

using GridManager = Dumux::PoreNetwork::GridManager<3>;
using NetworkData = Dumux::PoreNetwork::NetworkData<3>;
using Position = Dune::FieldVector<double, 3>;

int numFailures = 0;

void checkTrue(bool condition, const std::string& name)
{
    if (!condition)
    {
        std::cout << "FAILED " << name << "\n";
        ++numFailures;
    }
}

// a branched chain from the inlet (label 1) to the outlet (label 2) and a pair of pores not connected to it
NetworkData network()
{
    NetworkData data;
    data.porePositions = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                          {1.0, 1.0, 0.0}, {5.0, 5.0, 5.0}, {6.0, 5.0, 5.0}};
    data.throatPores = {{0, 1}, {1, 2}, {2, 3}, {4, 1}, {5, 6}};
    const std::array<double, 7> labels{1.0, -1.0, -1.0, 2.0, -1.0, -1.0, -1.0};
    data.vertexParameterNames = {"PoreInscribedRadius", "PoreVolume", "PoreLabel"};
    for (std::size_t i = 0; i < data.porePositions.size(); ++i)
        data.vertexParameters.push_back({0.1 + i/7.0, 1.0/(3.0 + i), labels[i]});
    data.elementParameterNames = {"ThroatInscribedRadius", "ThroatLength", "ThroatLabel"};
    for (std::size_t t = 0; t < data.throatPores.size(); ++t)
        data.elementParameters.push_back({0.01 + t/11.0, 1.0/(1.0 + t), -1.0});
    return data;
}

void writeDgf(const std::string& fileName, const NetworkData& data)
{
    std::ofstream file(fileName);
    file << std::setprecision(17) << "DGF\n% Vertex parameters:";
    for (const auto& name : data.vertexParameterNames)
        file << " " << name;
    file << "\n% Element parameters:";
    for (const auto& name : data.elementParameterNames)
        file << " " << name;
    file << "\nVertex\nparameters " << data.vertexParameterNames.size() << "\n";
    for (std::size_t i = 0; i < data.porePositions.size(); ++i)
    {
        file << data.porePositions[i];
        for (const double value : data.vertexParameters[i])
            file << " " << value;
        file << "\n";
    }
    file << "#\nSIMPLEX\nparameters " << data.elementParameterNames.size() << "\n";
    for (std::size_t t = 0; t < data.throatPores.size(); ++t)
    {
        file << data.throatPores[t][0] << " " << data.throatPores[t][1];
        for (const double value : data.elementParameters[t])
            file << " " << value;
        file << "\n";
    }
    file << "#\n";
}

std::size_t poreAt(const NetworkData& data, const Position& position)
{
    const auto it = std::find(data.porePositions.begin(), data.porePositions.end(), position);
    if (it == data.porePositions.end())
        DUNE_THROW(Dune::InvalidStateException, "No pore at " << position);
    return std::distance(data.porePositions.begin(), it);
}

//! check that the grid holds exactly the given pores and throats, and returns the number of vertices and elements
std::array<std::size_t, 2> checkGrid(GridManager& gridManager, const NetworkData& data, const std::string& name)
{
    const auto gridView = gridManager.grid().leafGridView();
    const auto gridData = gridManager.getGridData();
    for (const auto& vertex : vertices(gridView))
    {
        const auto i = poreAt(data, vertex.geometry().center());
        for (std::size_t k = 0; k < data.vertexParameterNames.size(); ++k)
            checkTrue(gridData->getParameter(vertex, data.vertexParameterNames[k]) == data.vertexParameters[i][k],
                      name + ": " + data.vertexParameterNames[k] + " of pore " + std::to_string(i));
    }

    std::map<std::array<std::size_t, 2>, std::size_t> throatOf;
    for (std::size_t t = 0; t < data.throatPores.size(); ++t)
        throatOf[{std::min(data.throatPores[t][0], data.throatPores[t][1]), std::max(data.throatPores[t][0], data.throatPores[t][1])}] = t;
    for (const auto& element : elements(gridView))
    {
        const auto i = poreAt(data, element.geometry().corner(0));
        const auto j = poreAt(data, element.geometry().corner(1));
        const auto t = throatOf.find({std::min(i, j), std::max(i, j)});
        checkTrue(t != throatOf.end(), name + ": element between pores " + std::to_string(i) + " and " + std::to_string(j) + " is a throat");
        if (t == throatOf.end())
            continue;
        for (std::size_t k = 0; k < data.elementParameterNames.size(); ++k)
            checkTrue(gridData->getParameter(element, data.elementParameterNames[k]) == data.elementParameters[t->second][k],
                      name + ": " + data.elementParameterNames[k] + " of throat " + std::to_string(t->second));
    }
    return {static_cast<std::size_t>(gridView.size(1)), static_cast<std::size_t>(gridView.size(0))};
}

template<class Function>
void checkThrows(const Function& function, const std::string& name)
{
    bool thrown = false;
    try { function(); }
    catch (const Dune::InvalidStateException&) { thrown = true; }
    checkTrue(thrown, name);
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    writeDgf("network_data.dgf", network());
    Dumux::Parameters::init([](auto& tree) {
        tree["Dgf.Grid.File"] = "network_data.dgf";
        tree["Sanitized.Grid.Sanitize"] = "true";
        tree["Sanitized.Grid.SanitationMode"] = "UsePoreLabels";
        tree["Sanitized.Grid.PruningSeedIndices"] = "1";
    });
    const auto data = network();

    GridManager given;
    given.init(data);
    const auto [numVertices, numElements] = checkGrid(given, data, "given data");
    checkTrue(numVertices == 7 && numElements == 5, "given data: all pores and throats");
    checkThrows([&]{ given.getGridData()->assignParameters(); }, "given data: no generated parameters");

    GridManager dgf;
    dgf.init("Dgf");
    const auto [numDgfVertices, numDgfElements] = checkGrid(dgf, data, "dgf file");
    checkTrue(numDgfVertices == numVertices && numDgfElements == numElements, "dgf file: same grid size");
    for (const auto& name : data.vertexParameterNames)
        checkTrue(dgf.getGridData()->parameterIndex(name) == given.getGridData()->parameterIndex(name), "same index of " + name);
    for (const auto& name : data.elementParameterNames)
        checkTrue(dgf.getGridData()->parameterIndex(name) == given.getGridData()->parameterIndex(name), "same index of " + name);

    GridManager sanitized;
    sanitized.init(data, "Sanitized");
    const auto [numSanitizedVertices, numSanitizedElements] = checkGrid(sanitized, data, "sanitized");
    checkTrue(numSanitizedVertices == 5 && numSanitizedElements == 4, "sanitized: the pair not connected to the inlet is removed");

    auto shortRow = data;
    shortRow.elementParameters[2].pop_back();
    checkThrows([&]{ GridManager m; m.init(shortRow); }, "throat with too few parameters");
    auto missingPore = data;
    missingPore.vertexParameters.pop_back();
    checkThrows([&]{ GridManager m; m.init(missingPore); }, "pore without parameters");
    auto unknownPore = data;
    unknownPore.throatPores[0][1] = 99;
    checkThrows([&]{ GridManager m; m.init(unknownPore); }, "throat to a pore that does not exist");

    if (numFailures > 0)
        DUNE_THROW(Dune::Exception, numFailures << " checks failed");
    std::cout << "all checks passed\n";
    return 0;
}
