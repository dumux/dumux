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
 * \brief Implements the notion of stencils for cell-centered models
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_STENCILS_HH
#define DUMUX_DISCRETIZATION_STAGGERED_STENCILS_HH

#include <set>
#include <dumux/implicit/staggered/properties.hh>

namespace Dumux
{

/*!
 * \brief Element-related stencils
 */
template<class TypeTag>
class StaggeredCellCenterStencils
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    // TODO a separate stencil class all stencils can derive from?
    using Stencil = std::vector<IndexType>;
public:
    //! Update the stencil. We expect a bound fvGeometry
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry)
    {

        const auto globalI = problem.elementMapper().index(element);

        cellCenterToCellCenterStencil_.clear();
        cellCenterToFaceStencil_.clear();

        // the cell center dof indices
        cellCenterToCellCenterStencil_.push_back(globalI);

        // loop over sub control faces
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary())
                cellCenterToCellCenterStencil_.push_back(scvf.outsideScvIdx());

            cellCenterToFaceStencil_.push_back(scvf.dofIndexSelf());
        }

        // make values unique
        std::sort(cellCenterToCellCenterStencil_.begin(), cellCenterToCellCenterStencil_.end());
        cellCenterToCellCenterStencil_.erase(std::unique(cellCenterToCellCenterStencil_.begin(), cellCenterToCellCenterStencil_.end()), cellCenterToCellCenterStencil_.end());

        std::sort(cellCenterToFaceStencil_.begin(), cellCenterToFaceStencil_.end());
        cellCenterToFaceStencil_.erase(std::unique(cellCenterToFaceStencil_.begin(), cellCenterToFaceStencil_.end()), cellCenterToFaceStencil_.end());
    }


    //! The full element stencil (all element this element is interacting with)
    const Stencil& cellCenterToCellCenterStencil() const
    {
        return cellCenterToCellCenterStencil_;
    }

    //! The full element stencil (all element this element is interacting with)
    const Stencil& cellCenterToFaceStencil() const
    {
        return cellCenterToFaceStencil_;
    }



private:
    Stencil cellCenterToCellCenterStencil_;
    Stencil cellCenterToFaceStencil_;
};

template<class TypeTag>
class StaggeredFaceStencils
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    // TODO a separate stencil class all stencils can derive from?
    using Stencil = std::vector<IndexType>;

public:
    void update(const Problem& problem,
                const SubControlVolumeFace& scvf,
                const FVElementGeometry& fvGeometry)
    {
        faceToCellCenterStencil_.clear();
        faceToFaceStencil_.clear();
        const int eIdx = scvf.insideScvIdx();

        faceToCellCenterStencil_.push_back(eIdx);

        faceToFaceStencil_.push_back(scvf.dofIndexSelf());
        faceToFaceStencil_.push_back(scvf.dofIndexOpposite());

        for(const auto& data : scvf.pairData())
        {
            auto& normalFace = fvGeometry.scvf(eIdx, data.localNormalFaceIdx);
            const auto outerParallelElementDofIdx = normalFace.outsideScvIdx();
            if(!normalFace.boundary())
                faceToCellCenterStencil_.push_back(outerParallelElementDofIdx);

            const auto& outerParallelFaceDofIdx = data.outerParallelFaceDofIdx;
            if(outerParallelFaceDofIdx >= 0)
                faceToFaceStencil_.push_back(outerParallelFaceDofIdx);

            faceToFaceStencil_.push_back(data.normalPair.first);
            if(!scvf.boundary())
                faceToFaceStencil_.push_back(data.normalPair.second);
        }
    }

     //! The full face stencil (all dofs this face is interacting with)
    const Stencil& faceToCellCenterStencil() const
    {
        return faceToCellCenterStencil_;
    }

     //! The full face stencil (all dofs this face is interacting with)
    const Stencil& faceToFaceStencil() const
    {
        return faceToFaceStencil_;
    }

private:
    Stencil faceToCellCenterStencil_;
    Stencil faceToFaceStencil_;
};


/*!
 * \ingroup CCModel
 * \brief The global stencil container class
 */
template<class TypeTag>
class StaggeredStencilsVector
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using StencilType = StaggeredCellCenterStencils<TypeTag>;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;


public:
    void update(Problem& problem)
    {
        problemPtr_ = &problem;

        const IndexType numElements = problem.gridView().size(0);
        const IndexType numFacets = problem.gridView().size(1);
        const IndexType numBoundaryFacets = problem.gridView().grid().numBoundarySegments();
        cellCenterStencils_.resize(numElements);
        faceStencils_.resize(2*numFacets - numBoundaryFacets);

        std::vector<Stencil> fullFaceToCellCenterStencils;
        fullFaceToCellCenterStencils.resize(numFacets);
        std::vector<Stencil> fullfaceToFaceStencils;
        fullfaceToFaceStencils.resize(numFacets);

        for (const auto& element : elements(problem.gridView()))
        {
            auto eIdx = problem.elementMapper().index(element);
            // restrict the FvGeometry locally and bind to the element
            auto fvGeometry = localView(problem.model().globalFvGeometry());
            fvGeometry.bindElement(element); // for TPFA bind element is enough
            cellCenterStencils_[eIdx].update(problem, element, fvGeometry);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                faceStencils_[scvf.index()].update(problem, scvf, fvGeometry);

                const IndexType idx = scvf.dofIndexSelf();

                const auto& faceToCellCenterStencil = faceStencils_[scvf.index()].faceToCellCenterStencil();
                fullFaceToCellCenterStencils[idx].insert(fullFaceToCellCenterStencils[idx].end(), faceToCellCenterStencil.begin(), faceToCellCenterStencil.end());
                std::sort(fullFaceToCellCenterStencils[idx].begin(), fullFaceToCellCenterStencils[idx].end());
                fullFaceToCellCenterStencils[idx].erase(std::unique(fullFaceToCellCenterStencils[idx].begin(), fullFaceToCellCenterStencils[idx].end()), fullFaceToCellCenterStencils[idx].end());

                const auto& faceToFaceStencil = faceStencils_[scvf.index()].faceToFaceStencil();
                fullfaceToFaceStencils[idx].insert(fullfaceToFaceStencils[idx].end(), faceToFaceStencil.begin(), faceToFaceStencil.end());
                std::sort(fullfaceToFaceStencils[idx].begin(), fullfaceToFaceStencils[idx].end());
                fullfaceToFaceStencils[idx].erase(std::unique(fullfaceToFaceStencils[idx].begin(), fullfaceToFaceStencils[idx].end()), fullfaceToFaceStencils[idx].end());
            }
        }
        // TODO: is this a good idea?
        fullFaceToCellCenterStencils_ = std::make_unique<decltype(fullFaceToCellCenterStencils)>(fullFaceToCellCenterStencils);
        fullfaceToFaceStencils_ = std::make_unique<decltype(fullfaceToFaceStencils)>(fullfaceToFaceStencils);
    }


    //! overload for elements
    auto& get(const Element& entity) const
    {
        return cellCenterStencils_[problemPtr_->elementMapper().index(entity)];
    }

    //! overload for faces
    auto& get(const SubControlVolumeFace& scvFace) const
    {
        const IndexType numElements = problemPtr_->gridView().size(0);
        return faceStencils_[scvFace.dofIndexSelf()];
    }

    /*!
    * \brief Returns the size of a complete face dof stencil
    */
    size_t fullFaceToCellCenterStencilSize(const int idx) const
    {
        assert(fullFaceToCellCenterStencils_ && "fullFaceToCellCenterStencils_ has already been called and deleted!");
        return fullFaceToCellCenterStencils_.get()[0][idx].size();
    }

    /*!
    * \brief Returns the size of a complete face dof stencil
    */
    size_t fullfaceToFaceStencilSize(const int idx) const
    {
        assert(fullfaceToFaceStencils_ && "fullfaceToFaceStencils_ has already been called and deleted!");
        return fullfaceToFaceStencils_.get()[0][idx].size();
    }

    /*!
    * \brief Returns a unique pointer to the complete face dof stencils which is used once for setting up the global matrix and deleted afterwards
    */
    auto getFullFaceToCellCenterStencilsPtr()
    {
        assert(fullFaceToCellCenterStencils_ && "fullFaceToCellCenterStencils_ has already been used and deleted!");
        return std::move(fullFaceToCellCenterStencils_);
    }

    /*!
    * \brief Returns a unique pointer to the complete face dof stencils which is used once for setting up the global matrix and deleted afterwards
    */
    auto getFullfaceToFaceStencilsPtr()
    {
        assert(fullfaceToFaceStencils_ && "fullfaceToFaceStencils_ has already been used and deleted!");
        return std::move(fullfaceToFaceStencils_);
    }


private:
    std::vector<StaggeredCellCenterStencils<TypeTag>> cellCenterStencils_;
    std::vector<StaggeredFaceStencils<TypeTag>> faceStencils_;
    std::unique_ptr<std::vector<Stencil>> fullFaceToCellCenterStencils_;
    std::unique_ptr<std::vector<Stencil>> fullfaceToFaceStencils_;


    const Problem* problemPtr_;
};

} // end namespace

#endif
