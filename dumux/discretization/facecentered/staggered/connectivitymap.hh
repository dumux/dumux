// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredConnectivityMap
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_CONNECTIVITY_MAP_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_CONNECTIVITY_MAP_HH

#include <algorithm>
#include <vector>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/gridenums.hh>
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

            assert(element.partitionType() == Dune::OverlapEntity);

            // restrict the FvGeometry locally and bind to the element
            fvGeometry.bind(element);
            // loop over sub control faces
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.isFrontal() && !scvf.boundary() && !scvf.processorBoundary())
                {
                    const auto& ownScv = fvGeometry.scv(scvf.insideScvIdx());
                    const auto& facet = element.template subEntity <1> (ownScv.indexInElement());
                    if (facet.partitionType() == Dune::BorderEntity)
                        map_[ownScv.index()].push_back(scvf.outsideScvIdx());
                }
            }
        }

        for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
        {
            fvGeometry.bind(element);

            // loop over sub control faces
            for (const auto& scvf : scvfs(fvGeometry))
            {
                assert(!scvf.processorBoundary());
                const auto& ownScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto ownDofIndex = ownScv.dofIndex();
                const auto ownScvIndex = ownScv.index();

                if (scvf.isFrontal())
                {
                    if (!scvf.boundary()) // opposite dof
                        map_[ownScvIndex].push_back(scvf.outsideScvIdx());
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
                    if (!scvf.boundary())
                        map_[ownScvIndex].push_back(scvf.outsideScvIdx());


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
                    map_[ownScvIndex].push_back(orthogonalScvf.insideScvIdx());
                    if (!orthogonalScvf.boundary())
                        map_[ownScvIndex].push_back(orthogonalScvf.outsideScvIdx());
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
