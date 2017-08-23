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
 * \brief Stores the face indices corresponding to the neighbors of an element
 *        that contribute to the derivative calculation. This is used for
 *        finite-volume schemes with symmetric sparsity pattern in the global matrix.
 */
#ifndef DUMUX_STAGGERED_ASSEMBLY_MAP_HH
#define DUMUX_STAGGERED_ASSEMBLY_MAP_HH

#include <dune/istl/bcrsmatrix.hh>

#include <dumux/implicit/properties.hh>

namespace Dumux
{

template<class TypeTag>
class StaggeredAssemblyMap
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using IndexType = typename GridView::IndexSet::IndexType;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using CellCenterIdxType = typename DofTypeIndices::CellCenterIdx;
    using FaceIdxType = typename DofTypeIndices::FaceIdx;

    using CellCenterToCellCenterMap = std::vector<std::vector<IndexType>>;
    using CellCenterToFaceMap = std::vector<std::vector<IndexType>>;
    using FaceToCellCenterMap = std::vector<std::vector<IndexType>>;
    using FaceToFaceMap = std::vector<std::vector<IndexType>>;

    using Stencil = std::vector<IndexType>;

public:

    /*!
     * \brief Initialize the AssemblyMap object.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(const Problem& problem)
    {
        const auto numDofsCC = problem.model().numCellCenterDofs();
        const auto numDofsFace = problem.model().numFaceDofs();
        const auto numBoundaryFacets = problem.model().fvGridGeometry().numBoundaryScvf();
        cellCenterToCellCenterMap_.resize(numDofsCC);
        cellCenterToFaceMap_.resize(numDofsCC);
        faceToCellCenterMap_.resize(2*numDofsFace - numBoundaryFacets);
        faceToFaceMap_.resize(2*numDofsFace - numBoundaryFacets);

        std::vector<Stencil> fullFaceToCellCenterStencils;
        fullFaceToCellCenterStencils.resize(numDofsFace);
        std::vector<Stencil> fullfaceToFaceStencils;
        fullfaceToFaceStencils.resize(numDofsFace);

        FluxVariables fluxVars;
        for(auto&& element: elements(problem.gridView()))
        {
            // restrict the FvGeometry locally and bind to the element
            auto fvGeometry = localView(problem.model().fvGridGeometry());
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto dofIdxCellCenter = problem.elementMapper().index(element);
                fluxVars.computeCellCenterToCellCenterStencil(cellCenterToCellCenterMap_[dofIdxCellCenter], problem, element, fvGeometry, scvf);
                fluxVars.computeCellCenterToFaceStencil(cellCenterToFaceMap_[dofIdxCellCenter], problem, element, fvGeometry, scvf);

                const auto scvfIdx = scvf.index();
                fluxVars.computeFaceToCellCenterStencil(faceToCellCenterMap_[scvfIdx],problem, fvGeometry, scvf);
                fluxVars.computeFaceToFaceStencil(faceToFaceMap_[scvfIdx],problem, fvGeometry, scvf);
            }
        }
    }

    const std::vector<IndexType>& operator() (CellCenterIdxType, CellCenterIdxType, const IndexType globalI) const
    {
        return cellCenterToCellCenterMap_[globalI];
    }

    const std::vector<IndexType>& operator() (CellCenterIdxType, FaceIdxType, const IndexType globalI) const
    {
        return cellCenterToFaceMap_[globalI];
    }

    const std::vector<IndexType>& operator() (FaceIdxType, CellCenterIdxType, const IndexType globalI) const
    {
        return faceToCellCenterMap_[globalI];
    }

    const std::vector<IndexType>& operator() (FaceIdxType, FaceIdxType, const IndexType globalI) const
    {
        return faceToFaceMap_[globalI];
    }

private:

    CellCenterToCellCenterMap cellCenterToCellCenterMap_;
    CellCenterToFaceMap cellCenterToFaceMap_;
    FaceToCellCenterMap faceToCellCenterMap_;
    FaceToFaceMap faceToFaceMap_;
};

} // end namespace Dumux

#endif
