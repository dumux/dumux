// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredConnectivityMap
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_CONNECTIVITY_MAP_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_CONNECTIVITY_MAP_HH

#include <algorithm>
#include <vector>

#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Stores the dof indices corresponding to the neighboring scvs
 *        that contribute to the derivative calculation.
 */
template<class GridGeometry>
class FaceCenteredStaggeredConnectivityMap
{
    using GridView = typename GridGeometry::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using Stencil = std::vector<GridIndexType>;
    using Map = std::vector<Stencil>;
    static constexpr bool useHigherOrder = GridGeometry::useHigherOrder;

public:

    //! Update the map and prepare the stencils
    void update(const GridGeometry& gridGeometry)
    {
        map_.clear();
        map_.resize(gridGeometry.numScv());

        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            if (element.partitionType() == Dune::InteriorEntity)
                continue;

            // restrict the FvGeometry locally and bind to the element
            fvGeometry.bind(element);
            // loop over sub control faces
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.isFrontal() && !scvf.boundary())
                {
                    const auto& ownScv = fvGeometry.scv(scvf.insideScvIdx());
                    const auto& facet = element.template subEntity <1> (ownScv.indexInElement());
                    if (facet.partitionType() == 1)
                    {
                        const auto& oppositeScv =  fvGeometry.scv(scvf.outsideScvIdx());
                        map_[ownScv.index()].push_back(oppositeScv.index());
                    }
                }
            }
        }

        for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
        {
            fvGeometry.bind(element);

            // loop over sub control faces
            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto& ownScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto ownDofIndex = ownScv.dofIndex();
                const auto ownScvIndex = ownScv.index();

                if (scvf.isFrontal())
                {
                    if (!scvf.boundary()) // opposite dof
                    {
                        map_[ownScvIndex].push_back(fvGeometry.scv(scvf.outsideScvIdx()).index());
                        if constexpr (useHigherOrder)
                        {
                            // add the forward SCV index, if possible
                            if (fvGeometry.hasForwardNeighbor(scvf))
                                map_[ownScvIndex].push_back(fvGeometry.forwardScvIdx(scvf));
                            // add the backward SCV index, if possible
                            if (fvGeometry.hasBackwardNeighbor(scvf))
                                map_[ownScvIndex].push_back(fvGeometry.backwardScvIdx(scvf));
                        }
                    }
                    else
                    {
                        // treat frontal faces on boundaries
                        for (const auto& scv : scvs(fvGeometry))
                        {
                            const auto otherDofIndex = scv.dofIndex();
                            if (ownDofIndex != otherDofIndex)
                                map_[scv.index()].push_back(ownScvIndex);
                        }
                    }

                }
                else // lateral face
                {
                    // the parallel DOF scv
                    // outsideScv   insideScv
                    // ###y######|s|#####x####
                    //           |c|
                    //           |v|
                    //           |f|
                    if (fvGeometry.hasParallelNeighbor(scvf))
                    {
                        map_[ownScvIndex].push_back(fvGeometry.parallelScvIdx(scvf));
                        if constexpr (useHigherOrder)
                        {
                            // add the second parallel SCV index, if possible
                            if (fvGeometry.hasSecondParallelNeighbor(scvf))
                                map_[ownScvIndex].push_back(fvGeometry.secondParallelScvIdx(scvf));
                        }
                    }

                    // the normal DOF scv
                    //
                    //     outsideScv
                    //
                    // |s|orthogonalScvf
                    // |c|
                    // |v| insideScv
                    // |f|
                    const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                    assert(orthogonalScvf.isLateral());
                    map_[ownScvIndex].push_back(fvGeometry.scv(orthogonalScvf.insideScvIdx()).index());
                    if (!orthogonalScvf.boundary())
                        map_[ownScvIndex].push_back(fvGeometry.scv(orthogonalScvf.outsideScvIdx()).index());
                }
            }
        }

        // make stencils unique
        for (auto& stencil : map_)
        {
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }
    }

    const Stencil& operator[] (const GridIndexType globalI) const
    { return map_[globalI]; }

private:
    Map map_;
};

} // end namespace Dumux

#endif
