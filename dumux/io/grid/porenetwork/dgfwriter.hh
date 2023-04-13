// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Write pore-network grids with attached data to dgf file
 */
#ifndef DUMUX_PORE_NETWORK_DGF_WRITER_HH
#define DUMUX_PORE_NETWORK_DGF_WRITER_HH

#include <string>
#include <iostream>
#include <fstream>

namespace Dumux::PoreNetwork {

 /*!
  * \ingroup PoreNetworkModels
  * \brief Write pore-network grids with attached data to dgf file
  */
template<class GridView, class GridData>
void writeDgf(const std::string& fileName, const GridView& gridView, const GridData& gridData)
{
    const auto someElement = *(elements(gridView).begin());
    const auto someVertex = *(vertices(gridView).begin());
    const auto numVertexParams = gridData.parameters(someVertex).size();
    const auto numElementParams = gridData.parameters(someElement).size();

    std::ofstream dgfFile(fileName);
    dgfFile << "DGF\nVertex % Coordinates, volumes and boundary flags of the pore bodies\nparameters " << numVertexParams << "\n";
    dgfFile << "% Vertex parameters: ";
    for (const auto& p : gridData.vertexParameterNames())
        dgfFile << p << " ";
    dgfFile << "\n% Element parameters: ";
    for (const auto& p : gridData.elementParameterNames())
        dgfFile << p << " ";
    dgfFile << std::endl;

    for (const auto& vertex : vertices(gridView))
    {
        dgfFile << vertex.geometry().center() << " ";

        const auto& params = gridData.parameters(vertex);
        for (int i = 0; i < params.size(); ++i)
        {
            dgfFile << params[i];

            if (i < params.size() - 1)
                dgfFile << " ";
        }

        dgfFile << std::endl;
    }

    dgfFile << "#\nSIMPLEX % Connections of the pore bodies (pore throats)\nparameters " << numElementParams << "\n";

    for (const auto& element : elements(gridView))
    {
        dgfFile << gridView.indexSet().subIndex(element, 0, 1) << " ";
        dgfFile << gridView.indexSet().subIndex(element, 1, 1) << " ";

        const auto& params = gridData.parameters(element);
        for (int i = 0; i < params.size(); ++i)
        {
            dgfFile << params[i];

            if (i < params.size() - 1)
                dgfFile << " ";
        }

        dgfFile << std::endl;
    }

    dgfFile << "#";
}

} // end namespace Dumux::PoreNetwork

#endif
