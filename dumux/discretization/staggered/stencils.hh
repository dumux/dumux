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
#include <dumux/implicit/cellcentered/properties.hh>

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
        elementStencil_.clear();

        // loop over sub control faces
        for (auto&& scvf : scvfs(fvGeometry))
        {
            FluxVariables fluxVars;
            const auto& centerStencil = fluxVars.computeCellCenterStencil(problem, element, fvGeometry, scvf);
            elementStencil_.insert(elementStencil_.end(), centerStencil.begin(), centerStencil.end());
        }
        // make values in elementstencil unique
        std::sort(elementStencil_.begin(), elementStencil_.end());
        elementStencil_.erase(std::unique(elementStencil_.begin(), elementStencil_.end()), elementStencil_.end());

        auto globalI = problem.elementMapper().index(element);
        neighborStencil_ = elementStencil_;

        // remove the element itself and possible ghost neighbors from the neighbor stencil
        neighborStencil_.erase(std::remove_if(neighborStencil_.begin(), neighborStencil_.end(),
                                             [globalI](int i){ return (i == globalI); }),
                               neighborStencil_.end());
    }

    //! The full element stencil (all element this element is interacting with)
    const Stencil& elementStencil() const
    {
        return elementStencil_;
    }

    //! The full element stencil without this element
    const Stencil& neighborStencil() const
    {
        return neighborStencil_;
    }

private:
    Stencil elementStencil_;
    Stencil neighborStencil_;
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
                const SubControlVolumeFace& scvf)
    {
        faceStencil_.clear();
        FluxVariables fluxVars;

        const auto& faceStencil = fluxVars.computeFaceStencil(problem, scvf);
        faceStencil_.insert(faceStencil_.end(), faceStencil.begin(), faceStencil.end());
        std::sort(faceStencil_.begin(), faceStencil_.end());
        faceStencil_.erase(std::unique(faceStencil_.begin(), faceStencil_.end()), faceStencil_.end());
    }

     //! The full face stencil (all dofs this face is interacting with)
    const Stencil& faceStencil() const
    {
        return faceStencil_;
    }

private:
    Stencil faceStencil_;
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
        elementStencils_.resize(numElements);
        faceStencils_.resize(2*numFacets - numBoundaryFacets);

        std::vector<Stencil> fullFaceStencils;
        fullFaceStencils.resize(numFacets);

        for (const auto& element : elements(problem.gridView()))
        {
            auto eIdx = problem.elementMapper().index(element);
            // restrict the FvGeometry locally and bind to the element
            auto fvGeometry = localView(problem.model().globalFvGeometry());
            fvGeometry.bindElement(element); // for TPFA bind element is enough
            elementStencils_[eIdx].update(problem, element, fvGeometry);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                faceStencils_[scvf.index()].update(problem, scvf);

                const IndexType idx = scvf.dofIndexSelf() - numElements;
                const auto& faceStencil = faceStencils_[scvf.index()].faceStencil();
                fullFaceStencils[idx].insert(fullFaceStencils[idx].end(), faceStencil.begin(), faceStencil.end());
                std::sort(fullFaceStencils[idx].begin(), fullFaceStencils[idx].end());
                fullFaceStencils[idx].erase(std::unique(fullFaceStencils[idx].begin(), fullFaceStencils[idx].end()), fullFaceStencils[idx].end());
            }
        }
        // TODO: is this a good idea?
        fullFaceDofStencils_ = std::make_unique<decltype(fullFaceStencils)>(fullFaceStencils);
    }


    //! overload for elements
    auto& get(const Element& entity) const
    {
        return elementStencils_[problemPtr_->elementMapper().index(entity)];
    }

    //! overload for faces
    auto& get(const SubControlVolumeFace& scvFace) const
    {
        const IndexType numElements = problemPtr_->gridView().size(0);
        return faceStencils_[scvFace.dofIndexSelf() - numElements];
    }

    /*!
    * \brief Returns the size of a complete face dof stencil
    */
    size_t completeFaceDofStencilSize(const int idx) const
    {
//         const IndexType numElements = problemPtr_->gridView().size(0);
        assert(fullFaceDofStencils_ && "fullFaceDofStencils_ has already been called and deleted!");
        return fullFaceDofStencils_.get()[0][idx/*-numElements*/].size();
        // TODO: why does this not work?
    }

    /*!
    * \brief Returns a unique pointer to the complete face dof stencils which is used once for setting up the global matrix and deleted afterwards
    */
    auto getFullFaceDofStencilsPtr()
    {
        assert(fullFaceDofStencils_ && "fullFaceDofStencils_ has already been used and deleted!");
        return std::move(fullFaceDofStencils_);
    }


private:
    std::vector<StaggeredCellCenterStencils<TypeTag>> elementStencils_;
    std::vector<StaggeredFaceStencils<TypeTag>> faceStencils_;
    std::unique_ptr<std::vector<Stencil>> fullFaceDofStencils_;


    const Problem* problemPtr_;
};

} // end namespace

#endif
