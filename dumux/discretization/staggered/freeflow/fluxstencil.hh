// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \copydoc Dumux::StaggeredNavierStokesFluxStencils
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_FLUX_STENCIL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_FLUX_STENCIL_HH

#include <vector>
#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \brief Computes the stencil for the Navier-Stokes specific staggered grid discretization.
*         For a full description of the stencils, please see the document under dumux/doc/docextra/staggered
 */
template<class TypeTag>
class StaggeredNavierStokesFluxStencils
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using Stencil = std::vector<IndexType>;

public:

    /*
     * \brief Computes the stencil for cell center dofs w.r.t to other cell center dofs.
     *        Basically, these are the dof indices of the neighboring elements plus the dof index of the element itself.
     */
    static void computeCellCenterToCellCenterStencil(Stencil& stencil,
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
    static void computeCellCenterToFaceStencil(Stencil& stencil,
                                               const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.dofIndex());
    }

    /*
     * \brief Computes the stencil for face dofs w.r.t to cell center dofs.
     *        Basically, these are the dof indices of the elements adjacent to the face and those of
     *        the elements adjacent to the faces parallel to the own face.
     */
    static void computeFaceToCellCenterStencil(Stencil& stencil,
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
                const auto outerParallelElementDofIdx = normalFace.outsideScvIdx();
                stencil.push_back(outerParallelElementDofIdx);
            }
        }
    }

    /*
     * \brief Computes the stencil for face dofs w.r.t to face dofs.
     *        For a full description of the stencil, please see the document under dumux/doc/docextra/staggered
     */
    static void computeFaceToFaceStencil(Stencil& stencil,
                                         const FVElementGeometry& fvGeometry,
                                         const SubControlVolumeFace& scvf)
    {
        // the first entries are always the face dofIdx itself and the one of the opposing face
        if(stencil.empty())
        {
            stencil.push_back(scvf.dofIndex());
            stencil.push_back(scvf.dofIndexOpposingFace());
        }

        for(const auto& data : scvf.pairData())
        {
            stencil.push_back(data.normalPair.first);
            const auto outerParallelFaceDofIdx = data.outerParallelFaceDofIdx;
            if(outerParallelFaceDofIdx >= 0)
                stencil.push_back(outerParallelFaceDofIdx);
            if(!scvf.boundary())
                stencil.push_back(data.normalPair.second);
        }
    }
};

} // end namespace Dumux

#endif //DUMUX_STAGGERED_NAVIERSTOKES_FLUX_STENCIL_HH
