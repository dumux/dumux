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
#ifndef DUMUX_DISCRETIZATION_BOX_VOLVARSVECTOR_HH
#define DUMUX_DISCRETIZATION_BOX_VOLVARSVECTOR_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool useOldSol, bool enableVolVarsCache>
class BoxVolumeVariablesVector
{};

// specialization in case of storing the volume variables
template<class TypeTag, bool useOldSol>
class BoxVolumeVariablesVector<TypeTag, useOldSol,/*enableVolVarCaching*/true> : public std::vector<typename GET_PROP_TYPE(TypeTag, VolumeVariables)>
{
    friend BoxVolumeVariablesVector<TypeTag, !useOldSol, true>;
    friend ImplicitModel<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;

    enum{ isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    BoxVolumeVariablesVector<TypeTag, useOldSol, true>& operator= (const BoxVolumeVariablesVector<TypeTag, useOldSol, true>& other) = default;

    BoxVolumeVariablesVector<TypeTag, useOldSol, true>& operator= (const BoxVolumeVariablesVector<TypeTag, !useOldSol, true>& other)
    {
        // do the copy
        numScvs_ = other.numScvs_;
        volumeVariables_ = other.volumeVariables_;

        // return the existing object
        return *this;
    };

public:
    void update(Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;

        numScvs_ = problem.model().fvGeometries().numScv();
        volumeVariables_.resize(numScvs_);
        for (const auto& element : elements(problem.gridView()))
        {
            problem.model().fvGeometries_().bindElement(element);
            const auto& fvGeometry = problem.model().fvGeometries(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                (*this)[scv].update(sol[scv.dofIndex()],
                                        problem,
                                        element,
                                        scv);
            }
        }
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    {
        const auto& scv = problem_().model().fvGeometries().subControlVolume(scvIdx);
        return (*this)[scv];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        const auto& scv = problem_().model().fvGeometries().subControlVolume(scvIdx);
        return (*this)[scv];
    }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        return volumeVariables_[scv.index()];
    }

    VolumeVariables& operator [](const SubControlVolume& scv)
    {
        return volumeVariables_[scv.index()];
    }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element) {}
    // In the box method, the vol vars within the elements adjacent to a vertex need to be bound
    void bind(const Vertex& vertex) {}
    // function to prepare the vol vars within the element
    void bindElement(const Element& element) {}

private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    IndexType eIdxBound_;
    IndexType numScvs_;
    std::vector<VolumeVariables> volumeVariables_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class BoxVolumeVariablesVector<TypeTag, /*isOldSol*/false, /*enableVolVarCaching*/false>
{
    // prev vol vars have to be a friend class in order for the assignment operator to work
    friend BoxVolumeVariablesVector<TypeTag, true, false>;
    friend ImplicitModel<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;

    BoxVolumeVariablesVector& operator= (const BoxVolumeVariablesVector<TypeTag, /*isOldSol*/false, /*enableVolVarCaching*/false>& other) = default;

    // operator curVolVars = prevVolVars
    void operator= (const BoxVolumeVariablesVector<TypeTag, /*isOldSol*/true, /*enableVolVarCaching*/false>& other)
    {
        eIdxBound_ = -1;
    }

    BoxVolumeVariablesVector() : problemPtr_(nullptr), eIdxBound_(-1) {}


public:

    void update(Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;
    }


    // Binding of an element, prepares the volume variables within the element stencil
    // called by the local jacobian to prepare element assembly
    // specialization for cc models
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox)>::type
    bind(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);
        if (eIdx == eIdxBound_ && volVarIndices_.size() > 1)
            return;

        eIdxBound_ = eIdx;

        // make sure the FVElementGeometry is bound to the element
        problem_().model().fvGeometries_().bind(element);
        const auto& fvGeometry = problem_().model().fvGeometries(eIdx);

        // stencil information
        const auto& elementStencil = problem_().model().stencils(element).elementStencil();
        const auto& neighborStencil = problem_().model().stencils(element).neighborStencil();

        // maximum number of possible vol vars to be created
        const auto numDofs = elementStencil.size();;// + fvGeometry.numScvf();

        // resize local containers to the required size
        volumeVariables_.resize(numDofs);
        volVarIndices_.resize(numDofs);
        int localIdx = 0;

        // update the volume variables of the element at hand
        const auto& scvI = problem_().model().fvGeometries().subControlVolume(eIdx);
        const auto& solI = problem_().model().curSol()[eIdx];
        // VolumeVariables tmp;
        // tmp.update(solI, problem_(), element, scvI);
        volumeVariables_[localIdx].update(solI, problem_(), element, scvI);
        volVarIndices_[localIdx] = scvI.index();
        // volumeVariables_.push_back(tmp);
        // volVarIndices_.push_back(scvI.index());
        localIdx++;

        // Update the volume variables of the neighboring elements
        for (auto globalJ : neighborStencil)
        {
            const auto& elementJ = problem_().model().fvGeometries().element(globalJ);
            const auto& scvJ = problem_().model().fvGeometries().subControlVolume(globalJ);
            const auto& solJ = problem_().model().curSol()[globalJ];
            // tmp.update(solJ, problem_(), elementJ, scvJ);
            // volumeVariables_.push_back(tmp);
            // volVarIndices_.push_back(scvJ.index());
            volumeVariables_[localIdx].update(solJ, problem_(), elementJ, scvJ);
            volVarIndices_[localIdx] = scvJ.index();
            localIdx++;
        }

        // Update boundary volume variables
        for (const auto& scvFace : scvfs(fvGeometry))
        {
            // if we are not on a boundary, skip the rest
            if (!scvFace.boundary())
                continue;

            // When complex boundary handling is inactive, we only use BC vol vars on pure Dirichlet boundaries
            const auto bcTypes = problem_().boundaryTypes(element, scvFace);
            if (/*TODO !GET_PROP_VALUE(TypeTag, BoundaryReconstruction) && */bcTypes.hasNeumann() || bcTypes.hasOutflow())
                continue;

            const auto dirichletPriVars = problem_().dirichlet(element, scvFace);
            // tmp.update(dirichletPriVars, problem_(), element, scvI);
            // volumeVariables_.push_back(tmp);
            // volVarIndices_.push_back(scvI.index());
            volumeVariables_.resize(localIdx+1);
            volVarIndices_.resize(localIdx+1);
            volumeVariables_[localIdx].update(dirichletPriVars, problem_(), element, scvI);
            volVarIndices_[localIdx] = scvFace.outsideScvIdx();
            localIdx++;
        }
    }

    // specialization for box models, simply forwards to the bindElement method
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox)>::type
    bind(const Element& element)
    {
        bindElement(element);
    }

    // Binding of an element, prepares only the volume variables of the element
    // specialization for cc models
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox)>::type
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

    // In the box method, the vol vars within the elements adjacent to a vertex need to be bound (TODO: IMPLEMENT)
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox)>::type
    bind(const Vertex& vertex) const {}

    // specialization for box models
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox)>::type
    bindElement(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);
        if (eIdx == eIdxBound_)
            return;

        eIdxBound_ = eIdx;

        // make sure the FVElementGeometry is bound to the element
        problem_().model().fvGeometries_().bind(element);
        const auto& fvGeometry = problem_().model().fvGeometries(element);

        // get stencil information
        const auto numDofs = element.subEntities(dim);

        // resize volume variables to the required size
        volumeVariables_.resize(numDofs);
        // volVarIndices_.resize(numDofs);

        int localIdx = 0;
        for (const auto& scv : scvs(fvGeometry))
        {
            // std::cout << "scv index: " << scv.index() << ", dofIdx: " << scv.dofIndex() << ", localIdx: " << scv.indexInElement() << std::endl;
            const auto& sol = problem_().model().curSol()[scv.dofIndex()];
            // let the interface solver update the volvars
            // volVarIndices_[localIdx] = scv.index();
            // TODO: INTERFACE SOLVER
            // problem_().model().boxInterfaceConditionSolver().updateScvVolVars(element, scv, sol);
            volumeVariables_[scv.indexInElement()].update(sol, problem_(), element, scv);
            localIdx++;
        }
        // std::cout << "finished\n";
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    {
        return volumeVariables_[scvIdx];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        return volumeVariables_[scvIdx];
    }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        return volumeVariables_[scv.indexInElement()];
    }

    VolumeVariables& operator [](const SubControlVolume& scv)
    {
        return volumeVariables_[scv.indexInElement()];
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
class BoxVolumeVariablesVector<TypeTag, /*isOldSol*/true, /*enableVolVarCaching*/false>
{
    // current vol vars have to be a friend class in order for the assignment operator to work
    friend BoxVolumeVariablesVector<TypeTag, false, false>;
    friend ImplicitModel<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;

    BoxVolumeVariablesVector& operator= (const BoxVolumeVariablesVector<TypeTag, /*isOldSol*/false, /*enableVolVarCaching*/false>& other)
    {
        release_();
        problemPtr_ = other.problemPtr_;
        eIdxBound_ = -1;
        return *this;
    }

    BoxVolumeVariablesVector() : problemPtr_(nullptr), eIdxBound_(-1) {}


public:

    void update(const Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;
    }

    // Binding of an element, prepares the volume variables of only the element
    // specialization for cc models
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox)>::type
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

    // specialization for box models
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox)>::type
    bindElement(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);
        if (eIdx == eIdxBound_)
            return;

        eIdxBound_ = eIdx;

        // make sure the FVElementGeometry is bound to the element
        problem_().model().fvGeometries_().bind(element);
        const auto& fvGeometry = problem_().model().fvGeometries(element);

        // get stencil information
        const auto numDofs = element.subEntities(dim);

        // resize volume variables to the required size
        volumeVariables_.resize(numDofs);
        // volVarIndices_.resize(numDofs);

        int localIdx = 0;
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& sol = problem_().model().prevSol()[scv.dofIndex()];
            // let the interface solver update the volvars
            // volVarIndices_[localIdx] = scv.index();
            //TODO: INTERFACE SOLVER?
            // problem_().model().boxInterfaceConditionSolver().updateScvVolVars(element, scv, sol);
            volumeVariables_[localIdx].update(sol, problem_(), element, scv);
            localIdx++;
        }
    }

    // In the box method, the vol vars within the elements adjacent to a vertex need to be bound (TODO: IMPLEMENT)
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox)>::type
    bind(const Vertex& vertex) const {}

    // const VolumeVariables& operator [](IndexType scvIdx) const
    // {
    //     return volumeVariables_[getLocalIdx_(scvIdx)];
    // }

    // VolumeVariables& operator [](IndexType scvIdx)
    // {
    //     return volumeVariables_[getLocalIdx_(scvIdx)];
    // }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        return volumeVariables_[scv.indexInElement()];
    }

    VolumeVariables& operator [](const SubControlVolume& scv)
    {
        return volumeVariables_[scv.indexInElement()];
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
