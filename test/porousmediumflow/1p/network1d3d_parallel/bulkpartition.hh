// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Partition of a network given by a distributed YaspGrid around it, as a bulk grid would give
 */
#ifndef DUMUX_TEST_1P_NETWORK_BULK_PARTITION_HH
#define DUMUX_TEST_1P_NETWORK_BULK_PARTITION_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <limits>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/multidomain/embedded/lowdimpartition.hh>

namespace Dumux {

template<class NetworkGrid>
std::vector<int> partitionByBoxAround(const NetworkGrid& network, int cellsPerDirection)
{
    constexpr int dimworld = NetworkGrid::dimensionworld;
    using Point = Dune::FieldVector<double, dimworld>;
    Point lo(std::numeric_limits<double>::max()), hi(std::numeric_limits<double>::lowest());
    for (const auto& v : vertices(network.leafGridView()))
    {
        const auto x = v.geometry().center();
        for (int d = 0; d < dimworld; ++d)
        {
            lo[d] = std::min(lo[d], x[d]);
            hi[d] = std::max(hi[d], x[d]);
        }
    }
    for (int d = 0; d < dimworld; ++d)
    {
        const double margin = 0.05*(hi[d] - lo[d]) + 1e-12;
        lo[d] -= margin;
        hi[d] += margin;
    }

    using Bulk = Dune::YaspGrid<dimworld, Dune::EquidistantOffsetCoordinates<double, dimworld>>;
    std::array<int, dimworld> cells;
    cells.fill(cellsPerDirection);
    const Bulk bulk(lo, hi, cells, std::bitset<dimworld>(), 1, Dune::MPIHelper::getCommunicator());
    return lowDimPartition(bulk.leafGridView(), network.leafGridView());
}

} // end namespace Dumux

#endif
