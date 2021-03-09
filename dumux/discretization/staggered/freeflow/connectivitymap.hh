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
 * \copydoc Dumux::StaggeredFreeFlowConnectivityMap
 */
#ifndef DUMUX_STAGGERED_FREEFLOW_CONNECTIVITY_MAP_HH
#define DUMUX_STAGGERED_FREEFLOW_CONNECTIVITY_MAP_HH

#include <vector>
#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Stores the dof indices corresponding to the neighboring cell centers and faces
 *        that contribute to the derivative calculation. Specialization for the staggered free flow model.
 */
template<class GridGeometry>
class StaggeredFreeFlowConnectivityMap
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    using CellCenterIdxType = typename GridGeometry::DofTypeIndices::CellCenterIdx;
    using FaceIdxType = typename GridGeometry::DofTypeIndices::FaceIdx;

    using SmallLocalIndex = typename IndexTraits<GridView>::SmallLocalIndex;

    using Stencil = std::vector<GridIndexType>;
    using Map = std::vector<Stencil>;

    static constexpr SmallLocalIndex upwindSchemeOrder = GridGeometry::upwindSchemeOrder;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

public:

    //! Update the map and prepare the stencils
    void update(const GridGeometry& gridGeometry)
    {
        const auto numDofsCC = gridGeometry.gridView().size(0);
        const auto numDofsFace = gridGeometry.gridView().size(1);
        const auto numBoundaryFacets = gridGeometry.numBoundaryScvf();
        cellCenterToCellCenterMap_.resize(numDofsCC);
        cellCenterToFaceMap_.resize(numDofsCC);
        faceToCellCenterMap_.resize(2*numDofsFace - numBoundaryFacets);
        faceToFaceMap_.resize(2*numDofsFace - numBoundaryFacets);

        for(auto&& element: elements(gridGeometry.gridView()))
        {
            // restrict the FvGeometry locally and bind to the element
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // handle the cell center dof stencils first
                const auto dofIdxCellCenter = gridGeometry.elementMapper().index(element);

                // the stencil for cell center dofs w.r.t. to other cell center dofs,
                // includes all neighboring element indices
                if (!scvf.boundary())
                    cellCenterToCellCenterMap_[dofIdxCellCenter].push_back(scvf.outsideScvIdx());

                // the stencil for cell center dofs w.r.t. face dofs, includes the face dof indices of the current element
                cellCenterToFaceMap_[dofIdxCellCenter].push_back(scvf.dofIndex());

                // handle the face dof stencils
                const auto scvfIdx = scvf.index();
                computeFaceToCellCenterStencil_(faceToCellCenterMap_[scvfIdx], fvGeometry, scvf);
                computeFaceToFaceStencil_(faceToFaceMap_[scvfIdx], fvGeometry, scvf);
            }
        }
    }

    //! Returns the stencil of a cell center dof w.r.t. other cell center dofs
    const Stencil& operator() (CellCenterIdxType, CellCenterIdxType, const GridIndexType globalI) const
    {
        return cellCenterToCellCenterMap_[globalI];
    }

    //! Returns the stencil of a cell center dof w.r.t. face dofs
    const Stencil& operator() (CellCenterIdxType, FaceIdxType, const GridIndexType globalI) const
    {
        return cellCenterToFaceMap_[globalI];
    }

    //! Returns the stencil of a face dof w.r.t. cell center dofs
    const Stencil& operator() (FaceIdxType, CellCenterIdxType, const GridIndexType globalI) const
    {
        return faceToCellCenterMap_[globalI];
    }

    //! Returns the stencil of a face dof w.r.t. other face dofs
    const Stencil& operator() (FaceIdxType, FaceIdxType, const GridIndexType globalI) const
    {
        return faceToFaceMap_[globalI];
    }

private:

    /*
     * \brief Computes the stencil for face dofs w.r.t to cell center dofs.
     *        Basically, these are the dof indices of the elements adjacent to the face and those of
     *        the elements adjacent to the faces parallel to the own face.
     */
    void computeFaceToCellCenterStencil_(Stencil& stencil,
                                         const FVElementGeometry& fvGeometry,
                                         const SubControlVolumeFace& scvf)
    {
        const auto eIdx = scvf.insideScvIdx();
        stencil.push_back(eIdx);

        for (const auto& data : scvf.pairData())
        {
            auto& lateralFace = fvGeometry.scvf(eIdx, data.localLateralFaceIdx);
            if (!lateralFace.boundary())
            {
                const auto firstParallelElementDofIdx = lateralFace.outsideScvIdx();
                stencil.push_back(firstParallelElementDofIdx);
            }
        }
    }

    /*
     * \brief Computes the stencil for face dofs w.r.t to face dofs.
     *        For a full description of the stencil, please see the document under dumux/doc/docextra/staggered
     */
    void computeFaceToFaceStencil_(Stencil& stencil,
                                   const FVElementGeometry& fvGeometry,
                                   const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.axisData().oppositeDof);
        addHigherOrderInAxisDofs_(scvf, stencil, std::integral_constant<bool, useHigherOrder>{});

        for (const auto& data : scvf.pairData())
        {
            // add normal dofs
            stencil.push_back(data.lateralPair.first);
            if (!scvf.boundary())
                stencil.push_back(data.lateralPair.second);

            // add parallel dofs
            for (SmallLocalIndex i = 0; i < upwindSchemeOrder; i++)
            {
                if (data.hasParallelNeighbor[i])
                    stencil.push_back(data.parallelDofs[i]);
            }
        }
    }

    void addHigherOrderInAxisDofs_(const SubControlVolumeFace& scvf, Stencil& stencil, std::false_type) {}

    void addHigherOrderInAxisDofs_(const SubControlVolumeFace& scvf, Stencil& stencil, std::true_type)
    {
        for (SmallLocalIndex i = 0; i < upwindSchemeOrder - 1; i++)
        {
            if (scvf.hasBackwardNeighbor(i))
                stencil.push_back(scvf.axisData().inAxisBackwardDofs[i]);

            if (scvf.hasForwardNeighbor(i))
                stencil.push_back(scvf.axisData().inAxisForwardDofs[i]);
        }
    }

    Map cellCenterToCellCenterMap_;
    Map cellCenterToFaceMap_;
    Map faceToCellCenterMap_;
    Map faceToFaceMap_;
};

} // end namespace Dumux

#endif
