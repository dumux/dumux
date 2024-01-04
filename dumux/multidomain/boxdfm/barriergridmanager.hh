// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup BoxDFMModel
 * \copydoc Dumux::BoxDfmBarrierGridManager.
 */

#ifndef DUMUX_MULTIDOMAIN_BOXDFM_BARRIER_GRID_MANAGER_HH
#define DUMUX_MULTIDOMAIN_BOXDFM_BARRIER_GRID_MANAGER_HH

#include <config.h>

#include <vector>
#include <algorithm>

#if HAVE_DUNE_SUBGRID

#include <dumux/io/grid/gridmanager_sub.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup BoxDFMModel
 * \brief Grid manager for constructing a grid consisting of only the barrier elements.
 */
template<class HostGrid>
class BoxDfmBarrierGridManager
{
public:
    using Grid = Dune::SubGrid<HostGrid::dimension, HostGrid>;

    template<class FractureIntersections>
    BoxDfmBarrierGridManager(HostGrid& host, const FractureIntersections& isections)
    {
        std::vector<bool> isBarrier(host.leafGridView().size(0), false);
        isections.visitBarriers([&] (const auto& e) {
            isBarrier[host.leafGridView().indexSet().index(e)] = true;
        });
        gridManager_.init(host, [&] (const auto& e) {
            return isBarrier[host.leafGridView().indexSet().index(e)];
        });
    }

    Grid& grid() { return gridManager_.grid(); }
    const Grid& grid() const { return gridManager_.grid(); }

private:
    GridManager<Grid> gridManager_;
};

} // end namespace Dumux

#else // HAVE_DUNE_SUBGRID

template<class... Args>
class BoxDfmBarrierGridManager
{ static_assert(false, "dune-subgrid required for this grid manager"); };

#endif // HAVE_DUNE_SUBGRID
#endif
