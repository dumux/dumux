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
 * Every vertex and element carries the parameters of the pore and throat at its position, also after writing
 * the grid to a DGF file and reading it back, and after grid sanitation removed a cluster not connected to the inlet.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <string>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/porenetwork/dgfwriter.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/grid/porenetwork/networkdata.hh>

namespace {

using GridManager = Dumux::PoreNetwork::GridManager<3>;
using NetworkData = Dumux::PoreNetwork::NetworkData<3>;

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

std::size_t poreAt(const NetworkData& data, const Dune::FieldVector<double, 3>& position)
{ return std::find(data.porePositions.begin(), data.porePositions.end(), position) - data.porePositions.begin(); }

//! the grid holds the given pores and throats with their parameters; returns its numbers of vertices and elements
std::array<int, 2> checkGrid(GridManager& gridManager, const NetworkData& data, const std::string& name)
{
    const auto gridView = gridManager.grid().leafGridView();
    const auto gridData = gridManager.getGridData();
    for (const auto& vertex : vertices(gridView))
    {
        const auto i = poreAt(data, vertex.geometry().center());
        checkTrue(i < data.porePositions.size() && gridData->parameters(vertex) == data.vertexParameters[i],
                  name + ": parameters of pore " + std::to_string(i));
    }
    for (const auto& element : elements(gridView))
    {
        const std::array<std::size_t, 2> pores{poreAt(data, element.geometry().corner(0)), poreAt(data, element.geometry().corner(1))};
        const auto t = std::find_if(data.throatPores.begin(), data.throatPores.end(), [&](const auto& p) {
            return p == pores || p == std::array<std::size_t, 2>{pores[1], pores[0]};
        }) - data.throatPores.begin();
        checkTrue(t < static_cast<long>(data.throatPores.size()) && gridData->parameters(element) == data.elementParameters[t],
                  name + ": element between pores " + std::to_string(pores[0]) + " and " + std::to_string(pores[1]));
    }
    return {gridView.size(1), gridView.size(0)};
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
    Dumux::Parameters::init([](auto& tree) {
        tree["Dgf.Grid.File"] = "network_data.dgf";
        tree["Sanitized.Grid.Sanitize"] = "true";
        tree["Sanitized.Grid.SanitationMode"] = "UsePoreLabels";
        tree["Sanitized.Grid.PruningSeedIndices"] = "1";
    });
    const auto data = network();

    GridManager given;
    given.init(data);
    checkTrue(checkGrid(given, data, "given data") == std::array<int, 2>{7, 5}, "given data: all pores and throats");
    checkThrows([&]{ given.getGridData()->assignParameters(); }, "given data: no generated parameters");

    Dumux::PoreNetwork::writeDgf("network_data.dgf", given.grid().leafGridView(), *given.getGridData());
    GridManager dgf;
    dgf.init("Dgf");
    checkTrue(checkGrid(dgf, data, "dgf file") == std::array<int, 2>{7, 5}, "dgf file: all pores and throats");

    GridManager sanitized;
    sanitized.init(data, "Sanitized");
    checkTrue(checkGrid(sanitized, data, "sanitized") == std::array<int, 2>{5, 4},
              "sanitized: the pair not connected to the inlet is removed");

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
