// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModelGrids
 * \brief Pores and throats of a pore network with their parameters, to create a grid without a grid file
 */
#ifndef DUMUX_IO_PORENETWORKGRID_NETWORK_DATA_HH
#define DUMUX_IO_PORENETWORKGRID_NETWORK_DATA_HH

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModelGrids
 * \brief Pores and throats of a pore network with their parameters, to create a grid without a grid file
 *
 * Throat t connects the pores throatPores[t]. The parameters of pore i are vertexParameters[i] in the order
 * of vertexParameterNames, those of throat t are elementParameters[t] in the order of elementParameterNames.
 * The names are those of a DGF file, e.g. PoreInscribedRadius, PoreVolume, PoreLabel, ThroatInscribedRadius,
 * ThroatLength and ThroatLabel.
 */
template<int dimWorld, class Scalar = double>
struct NetworkData
{
    std::vector<Dune::FieldVector<Scalar, dimWorld>> porePositions;
    std::vector<std::array<std::size_t, 2>> throatPores;
    std::vector<std::string> vertexParameterNames;
    std::vector<std::vector<Scalar>> vertexParameters;
    std::vector<std::string> elementParameterNames;
    std::vector<std::vector<Scalar>> elementParameters;
};

} // end namespace Dumux::PoreNetwork

#endif
