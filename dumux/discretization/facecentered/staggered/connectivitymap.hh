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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredConnectivityMap
 */
#ifndef DUMUX_FACECENTERED_STAGGERED_CONNECTIVITY_MAP_HH
#define DUMUX_FACECENTERED_STAGGERED_CONNECTIVITY_MAP_HH


#include <algorithm>
#include <vector>
#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Stores the dof indices corresponding to the neighboring scvs
 *        that contribute to the derivative calculation.
 */
template<class GridGeometry, int upwindSchemeOrder = 1>
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

        for (const auto& element: elements(gridGeometry.gridView()))
        {
            // restrict the FvGeometry locally and bind to the element
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto& ownScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto ownDofIndex = ownScv.dofIndex();
                const auto ownScvIndex = ownScv.index();

                if (scvf.isFrontal())
                {
                    if (!scvf.boundary()) // opposite dof
                         map_[ownScvIndex].push_back(fvGeometry.scv(scvf.outsideScvIdx()).index());
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
                        map_[ownScvIndex].push_back(fvGeometry.scv(scvf.outsideScvIdx()).index());


                    // the normal DOF scv
                    //
                    //     outsideScv
                    //
                    // |s|orthogonalScvf
                    // |c|
                    // |v| insideScv
                    // |f|
                    const auto& orthogonalScvf = fvGeometry.scvfWithCommonEntity(scvf);
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

#endif // DUMUX_FACECENTERED_STAGGERED_CONNECTIVITY_MAP_HH
