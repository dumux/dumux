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
 * \brief Base class for the volume variables vector
 */
#ifndef DUMUX_IMPLICIT_VOLVARSVECTOR_HH
#define DUMUX_IMPLICIT_VOLVARSVECTOR_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool useOldSol, bool enableVolVarsCache>
class VolumeVariablesVector
{};

// specialization in case of storing the volume variables
template<class TypeTag, bool useOldSol>
class VolumeVariablesVector<TypeTag, useOldSol,/*enableVolVarCaching*/true> : public std::vector<typename GET_PROP_TYPE(TypeTag, VolumeVariables)>
{
    friend VolumeVariablesVector<TypeTag, !useOldSol, true>;
    friend ImplicitModel<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;

    enum{ isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    VolumeVariablesVector<TypeTag, useOldSol, true>& operator= (const VolumeVariablesVector<TypeTag, useOldSol, true>& other) = default;

    VolumeVariablesVector<TypeTag, useOldSol, true>& operator= (const VolumeVariablesVector<TypeTag, !useOldSol, true>& other)
    {
        // do the copy
        numScvs_ = other.numScvs_;
        numBoundaryScvs_ = other.numBoundaryScvs_;
        volumeVariables_ = other.volumeVariables_;
        boundaryVolumeVariables_ = other.boundaryVolumeVariables_;

        // return the existing object
        return *this;
    };

public:
    void update(Problem& problem, const SolutionVector& sol)
    {
        numScvs_ = problem.model().fvGeometries().numScv();
        if (!isBox)
            numBoundaryScvs_ = problem.model().fvGeometries().numBoundaryScvf();
        else
            numBoundaryScvs_ = 0;

        volumeVariables_.resize(numScvs_);
        boundaryVolumeVariables_.resize(numBoundaryScvs_);
        for (const auto& element : elements(problem.gridView()))
        {
            problem.model().fvGeometries_().bind(element);
            const auto& fvGeometry = problem.model().fvGeometries(element);
            for (const auto& scv : fvGeometry.scvs())
            {
                (*this)[scv.index()].update(sol[scv.dofIndex()],
                                            problem,
                                            element,
                                            scv);
            }

            if (!isBox)
            {
                for (const auto& scvFace : fvGeometry.scvfs())
                {
                    // if we are not on a boundary, skip the rest
                    if (!scvFace.boundary())
                        continue;

                    // When complex boundary handling is inactive, we only use BC vol vars on pure Dirichlet boundaries
                    const auto bcTypes = problem.boundaryTypes(element, scvFace);
                    if (/*TODO !GET_PROP_VALUE(TypeTag, BoundaryReconstruction) && */!(bcTypes.hasDirichlet() && !bcTypes.hasNeumann()))
                        continue;

                    const auto insideScvIdx = scvFace.insideScvIdx();
                    const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
                    const auto dirichletPriVars = problem.dirichlet(element, scvFace);

                    (*this)[scvFace.outsideScvIdx()].update(dirichletPriVars, problem, element, insideScv);
                }
            }
        }
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    {
        assert(scvIdx < numScvs_ + numBoundaryScvs_);

        if (scvIdx < numScvs_)
            return volumeVariables_[scvIdx];
        else
            return boundaryVolumeVariables_[scvIdx-numScvs_];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        assert(scvIdx < numScvs_ + numBoundaryScvs_);

        if (scvIdx < numScvs_)
            return volumeVariables_[scvIdx];
        else
            return boundaryVolumeVariables_[scvIdx-numScvs_];
    }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element) const {}
    // function to prepare the vol vars within the element
    void bindElement(const Element& element) const {}

private:
    IndexType numScvs_;
    IndexType numBoundaryScvs_;
    std::vector<VolumeVariables> volumeVariables_;
    std::vector<VolumeVariables> boundaryVolumeVariables_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class VolumeVariablesVector<TypeTag, /*isOldSol*/false, /*enableVolVarCaching*/false>
{
    // prev vol vars have to be a friend class in order for the assignment operator to work
    friend VolumeVariablesVector<TypeTag, true, false>;
    friend ImplicitModel<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;

    VolumeVariablesVector& operator= (const VolumeVariablesVector<TypeTag, /*isOldSol*/false, /*enableVolVarCaching*/false>& other) = default;

    void operator= (const VolumeVariablesVector<TypeTag, /*isOldSol*/true, /*enableVolVarCaching*/false>& other)
    {
        eIdxBound_ = -1;
    }

    VolumeVariablesVector() : problemPtr_(nullptr), eIdxBound_(-1) {}


public:

    void update(Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;
    }


    // Binding of an element, prepares the volume variables within the element stencil
    // called by the local jacobian to prepare element assembly
    // specialization for cc models
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(TypeTag, ImplicitIsBox)>::type
    bind(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);
        if (eIdx == eIdxBound_ && volVarIndices_.size() > 1)
            return;

        eIdxBound_ = eIdx;

        // make sure the FVElementGeometry is bound to the element
        problem_().model().fvGeometries_().bind(element);
        const auto& fvGeometry = problem_().model().fvGeometries(eIdx);

        // get stencil information
        const auto& elementStencil = problem_().model().stencils(element).elementStencil();
        const auto& neighborStencil = problem_().model().stencils(element).neighborStencil();
        const auto numDofs = elementStencil.size();

        // resize volume variables to the required size
        volumeVariables_.resize(numDofs);
        volVarIndices_.resize(numDofs);
        int localIdx = 0;

        // update the volume variables of the element at hand
        const auto& scvI = problem_().model().fvGeometries().subControlVolume(eIdx);
        const auto& solI = problem_().model().curSol()[eIdx];
        volumeVariables_[localIdx].update(solI, problem_(), element, scvI);
        volVarIndices_[localIdx] = scvI.index();
        localIdx++;

        // Update the volume variables of the neighboring elements and the boundary
        for (auto globalJ : neighborStencil)
        {
            const auto& elementJ = problem_().model().fvGeometries().element(globalJ);
            const auto& scvJ = problem_().model().fvGeometries().subControlVolume(globalJ);
            const auto& solJ = problem_().model().curSol()[globalJ];
            volumeVariables_[localIdx].update(solJ, problem_(), elementJ, scvJ);
            volVarIndices_[localIdx] = scvJ.index();
            localIdx++;
        }

        // Check for boundary volume variables
        for (const auto& scvFace : fvGeometry.scvfs())
        {
            // if we are not on a boundary, skip the rest
            if (!scvFace.boundary())
                continue;

            // When complex boundary handling is inactive, we only use BC vol vars on pure Dirichlet boundaries
            const auto bcTypes = problem_().boundaryTypes(element, scvFace);
            if (/*TODO !GET_PROP_VALUE(TypeTag, BoundaryReconstruction) && */!(bcTypes.hasDirichlet() && !bcTypes.hasNeumann()))
                continue;

            volumeVariables_.resize(localIdx+1);
            volVarIndices_.resize(localIdx+1);
            const auto dirichletPriVars = problem_().dirichlet(element, scvFace);
            volumeVariables_[localIdx].update(dirichletPriVars, problem_(), element, scvI);
            volVarIndices_[localIdx] = scvFace.outsideScvIdx();
            localIdx++;
        }
    }

    // Binding of an element, prepares only the volume variables of the element
    // specialization for cc models
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(TypeTag, ImplicitIsBox)>::type
    bindElement(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);
        if (eIdx == eIdxBound_ || std::find(volVarIndices_.begin(), volVarIndices_.end(), eIdx) != volVarIndices_.end())
            return;

        volumeVariables_.resize(1);
        volVarIndices_.resize(1);
        eIdxBound_ = eIdx;

        // make sure the FVElementGeometry is bound to the element
        problem_().model().fvGeometries_().bindElement(element);

        // update the volume variables of the element at hand
        const auto& scv = problem_().model().fvGeometries().subControlVolume(eIdx);
        const auto& sol = problem_().model().curSol()[eIdx];
        volumeVariables_[0].update(sol, problem_(), element, scv);
        volVarIndices_[0] = scv.index();
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    {
        return volumeVariables_[getLocalIdx_(scvIdx)];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        return volumeVariables_[getLocalIdx_(scvIdx)];
    }

private:

    void release_()
    {
        volumeVariables_.clear();
        volVarIndices_.clear();
    }

    const int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(volVarIndices_.begin(), volVarIndices_.end(), volVarIdx);

        if (it != volVarIndices_.end())
            return std::distance(volVarIndices_.begin(), it);
        else
            DUNE_THROW(Dune::InvalidStateException, "Could not find the current volume variables for volVarIdx = " << volVarIdx <<
                                                    ", make sure to properly bind the volume variables to the element before using them");
    }

    Problem& problem_() const
    { return *problemPtr_;}

    Problem* problemPtr_;
    IndexType eIdxBound_;
    std::vector<IndexType> volVarIndices_;
    std::vector<VolumeVariables> volumeVariables_;
};

// Specialization when the previous volume variables are not stored
template<class TypeTag>
class VolumeVariablesVector<TypeTag, /*isOldSol*/true, /*enableVolVarCaching*/false>
{
    // current vol vars have to be a friend class in order for the assignment operator to work
    friend VolumeVariablesVector<TypeTag, false, false>;
    friend ImplicitModel<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;

    VolumeVariablesVector& operator= (const VolumeVariablesVector<TypeTag, /*isOldSol*/false, /*enableVolVarCaching*/false>& other)
    {
        release_();
        problemPtr_ = other.problemPtr_;
        eIdxBound_ = -1;
        return *this;
    }

    VolumeVariablesVector() : problemPtr_(nullptr), eIdxBound_(-1) {}


public:

    void update(const Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;
    }

    // Binding of an element, prepares the volume variables of only the element
    // specialization for cc models
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(TypeTag, ImplicitIsBox)>::type
    bindElement(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);
        if (eIdx == eIdxBound_)
            return;

        volumeVariables_.resize(1);
        volVarIndices_.resize(1);
        eIdxBound_ = eIdx;

        // make sure FVElementGeometry is bound to the element
        problem_().model().fvGeometries_().bindElement(element);

        // update the volume variables of the element at hand
        const auto& scv = problem_().model().fvGeometries().subControlVolume(eIdx);
        const auto& sol = problem_().model().prevSol()[eIdx];
        volumeVariables_[0].update(sol, problem_(), element, scv);
        volVarIndices_[0] = scv.index();
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    {
        return volumeVariables_[getLocalIdx_(scvIdx)];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        return volumeVariables_[getLocalIdx_(scvIdx)];
    }

private:

    void release_()
    {
        volumeVariables_.clear();
        volVarIndices_.clear();
    }

    const int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(volVarIndices_.begin(), volVarIndices_.end(), volVarIdx);

        if (it != volVarIndices_.end())
            return std::distance(volVarIndices_.begin(), it);
        else
            DUNE_THROW(Dune::InvalidStateException, "Could not find the previous volume variables for volVarIdx = " << volVarIdx <<
                                                    ", make sure to properly bind the volume variables to the element before using them");
    }

    Problem& problem_() const
    { return *problemPtr_;}

    Problem* problemPtr_;
    IndexType eIdxBound_;
    std::vector<IndexType> volVarIndices_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace

#endif
