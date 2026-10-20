// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCDiscretization
 * \brief Connectivity map for cell-centered schemes whose flux over a face depends on
 *        the two adjacent elements and on all face neighbors of both
 */
#ifndef DUMUX_CC_TWO_RING_CONNECTIVITY_MAP_HH
#define DUMUX_CC_TWO_RING_CONNECTIVITY_MAP_HH

#include <vector>
#include <algorithm>

#include <dune/common/exceptions.hh>

#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Connectivity map for cell-centered schemes whose flux over a face depends on
 *        the two adjacent elements and on all face neighbors of both, as for fluxes
 *        between states reconstructed with cell gradients from the face neighbors.
 *        The flux over a boundary face depends on the inside element and its face neighbors.
 *
 * For each element i, the map lists the elements j (other than i) with at least one scvf
 * whose flux depends on the degrees of freedom of i, together with these scvfs of j.
 *
 * \tparam GridGeometry the finite volume grid geometry
 */
template<class GridGeometry>
class CCTwoRingConnectivityMap
{
    using GridView = typename GridGeometry::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    static constexpr int maxElemStencilSize = GridGeometry::maxElementStencilSize;

    struct DataJ
    {
        GridIndexType globalJ;
        std::vector<GridIndexType> scvfsJ;
    };

public:

    /*!
     * \brief Initialize the ConnectivityMap object.
     *
     * \param gridGeometry The grid's finite volume geometry.
     */
    void update(const GridGeometry& gridGeometry)
    {
        const auto numElements = gridGeometry.gridView().size(0);
        std::vector<std::vector<GridIndexType>> neighbors(numElements);

        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
                if (!scvf.boundary())
                    for (std::size_t k = 0; k < scvf.numOutsideScvs(); ++k)
                        neighbors[eIdx].push_back(scvf.outsideScvIdx(k));
        }

        map_.clear();
        map_.resize(numElements);

        std::vector<GridIndexType> stencil;
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto globalJ = gridGeometry.elementMapper().index(element);
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                stencil = neighbors[globalJ];
                if (!scvf.boundary())
                    for (std::size_t k = 0; k < scvf.numOutsideScvs(); ++k)
                    {
                        const auto& neighborsK = neighbors[scvf.outsideScvIdx(k)];
                        stencil.insert(stencil.end(), neighborsK.begin(), neighborsK.end());
                    }

                std::sort(stencil.begin(), stencil.end());
                stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());

                // the scvfs of element j are visited consecutively, so an existing entry for j is the last one
                for (const auto globalI : stencil)
                {
                    if (globalI == globalJ)
                        continue;

                    auto& dataJForI = map_[globalI];
                    if (dataJForI.empty() || dataJForI.back().globalJ != globalJ)
                        dataJForI.push_back(DataJ{globalJ, {scvf.index()}});
                    else
                        dataJForI.back().scvfsJ.push_back(scvf.index());
                }
            }
        }

        for (const auto& dataJForI : map_)
            if (dataJForI.size() > static_cast<std::size_t>(maxElemStencilSize - 1))
                DUNE_THROW(Dune::InvalidStateException, "Maximum admissible stencil size (" << maxElemStencilSize-1
                                                         << ") is surpassed (" << dataJForI.size() << "). "
                                                         << "Please adjust the GridGeometry traits accordingly!");
    }

    const std::vector<DataJ>& operator[] (const GridIndexType globalI) const
    { return map_[globalI]; }

private:
    std::vector<std::vector<DataJ>> map_;
};

} // end namespace Dumux

#endif
