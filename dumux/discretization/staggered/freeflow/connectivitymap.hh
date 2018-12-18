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
template<class FVGridGeometry>
class StaggeredFreeFlowConnectivityMap
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    using CellCenterIdxType = typename FVGridGeometry::DofTypeIndices::CellCenterIdx;
    using FaceIdxType = typename FVGridGeometry::DofTypeIndices::FaceIdx;

    using CellCenterToCellCenterMap = std::vector<std::vector<GridIndexType>>;
    using CellCenterToFaceMap = std::vector<std::vector<GridIndexType>>;
    using FaceToCellCenterMap = std::vector<std::vector<GridIndexType>>;
    using FaceToFaceMap = std::vector<std::vector<GridIndexType>>;

    using Stencil = std::vector<GridIndexType>;

public:

    //! Update the map and prepare the stencils
    void update(const FVGridGeometry& fvGridGeometry)
    {
        const auto numDofsCC = fvGridGeometry.gridView().size(0);
        const auto numDofsFace = fvGridGeometry.gridView().size(1);
        const auto numBoundaryFacets = fvGridGeometry.numBoundaryScvf();
        cellCenterToCellCenterMap_.resize(numDofsCC);
        cellCenterToFaceMap_.resize(numDofsCC);
        faceToCellCenterMap_.resize(2*numDofsFace - numBoundaryFacets);
        faceToFaceMap_.resize(2*numDofsFace - numBoundaryFacets);

        std::vector<Stencil> fullFaceToCellCenterStencils;
        fullFaceToCellCenterStencils.resize(numDofsFace);
        std::vector<Stencil> fullfaceToFaceStencils;
        fullfaceToFaceStencils.resize(numDofsFace);

        for(auto&& element: elements(fvGridGeometry.gridView()))
        {
            // restrict the FvGeometry locally and bind to the element
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto dofIdxCellCenter = fvGridGeometry.elementMapper().index(element);
                computeCellCenterToCellCenterStencil_(cellCenterToCellCenterMap_[dofIdxCellCenter], element, fvGeometry, scvf);
                computeCellCenterToFaceStencil_(cellCenterToFaceMap_[dofIdxCellCenter], element, fvGeometry, scvf);

                const auto scvfIdx = scvf.index();
                computeFaceToCellCenterStencil_(faceToCellCenterMap_[scvfIdx], fvGeometry, scvf);
                computeFaceToFaceStencil_(faceToFaceMap_[scvfIdx], fvGeometry, scvf);
            }
        }
    }

    //! Returns the stencil of a cell center dof w.r.t. other cell center dofs
    const std::vector<GridIndexType>& operator() (CellCenterIdxType, CellCenterIdxType, const GridIndexType globalI) const
    {
        return cellCenterToCellCenterMap_[globalI];
    }

    //! Returns the stencil of a cell center dof w.r.t. face dofs
    const std::vector<GridIndexType>& operator() (CellCenterIdxType, FaceIdxType, const GridIndexType globalI) const
    {
        return cellCenterToFaceMap_[globalI];
    }

    //! Returns the stencil of a face dof w.r.t. cell center dofs
    const std::vector<GridIndexType>& operator() (FaceIdxType, CellCenterIdxType, const GridIndexType globalI) const
    {
        return faceToCellCenterMap_[globalI];
    }

    //! Returns the stencil of a face dof w.r.t. other face dofs
    const std::vector<GridIndexType>& operator() (FaceIdxType, FaceIdxType, const GridIndexType globalI) const
    {
        return faceToFaceMap_[globalI];
    }

private:

    /*
     * \brief Computes the stencil for cell center dofs w.r.t to other cell center dofs.
     *        Basically, these are the dof indices of the neighboring elements plus the dof index of the element itself.
     */
    void computeCellCenterToCellCenterStencil_(Stencil& stencil,
                                               const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const SubControlVolumeFace& scvf)
    {
        // the first entry is always the cc dofIdx itself
        if(stencil.empty())
            stencil.push_back(scvf.insideScvIdx());
        if(!scvf.boundary())
            stencil.push_back(scvf.outsideScvIdx());
    }

    /*
     * \brief Computes the stencil for cell center dofs w.r.t to face dofs.
     *        Basically, these are the dof indices of the element's faces.
     */
    void computeCellCenterToFaceStencil_(Stencil& stencil,
                                         const Element& element,
                                         const FVElementGeometry& fvGeometry,
                                         const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.axisData().selfDof);
    }

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
        stencil.push_back(scvf.insideScvIdx());

        for(const auto& data : scvf.pairData())
        {
            auto& normalFace = fvGeometry.scvf(eIdx, data.localNormalFaceIdx);
            if(!normalFace.boundary())
            {
                const auto firstParallelElementDofIdx = normalFace.outsideScvIdx();
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
        if(stencil.empty())
        {
            for(int i = 0; i < scvf.axisData().inAxisBackwardDofs.size(); i++)
            {
                if(scvf.hasBackwardNeighbor(i))
                {
                    stencil.push_back(scvf.axisData().inAxisBackwardDofs[i]);
                }
            }

            stencil.push_back(scvf.axisData().selfDof);
            stencil.push_back(scvf.axisData().oppositeDof);

            for(int i = 0; i < scvf.axisData().inAxisForwardDofs.size(); i++)
            {
                if(scvf.hasForwardNeighbor(i))
                {
                    stencil.push_back(scvf.axisData().inAxisForwardDofs[i]);
                }
            }
        }

        for(const auto& data : scvf.pairData())
        {
            // add normal dofs
            stencil.push_back(data.normalPair.first);
            if(!scvf.boundary())
                stencil.push_back(data.normalPair.second);

            // add parallel dofs
            for (int i = 0; i < data.parallelDofs.size(); i++)
            {
                if(!(data.parallelDofs[i] < 0))
                {
                    stencil.push_back(data.parallelDofs[i]);
                }
            }
        }
    }

    CellCenterToCellCenterMap cellCenterToCellCenterMap_;
    CellCenterToFaceMap cellCenterToFaceMap_;
    FaceToCellCenterMap faceToCellCenterMap_;
    FaceToFaceMap faceToFaceMap_;
};

} // end namespace Dumux

#endif // DUMUX_STAGGERED_FREEFLOW_CONNECTIVITY_MAP_HH
