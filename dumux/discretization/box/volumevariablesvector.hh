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
class BoxVolumeVariablesVector<TypeTag, useOldSol,/*enableVolVarCaching*/true>
{
    friend BoxVolumeVariablesVector<TypeTag, !useOldSol, true>;
    friend ImplicitModel<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    enum{ isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    BoxVolumeVariablesVector<TypeTag, useOldSol, true>& operator= (const BoxVolumeVariablesVector<TypeTag, useOldSol, true>& other) = default;

    BoxVolumeVariablesVector<TypeTag, useOldSol, true>& operator= (const BoxVolumeVariablesVector<TypeTag, !useOldSol, true>& other)
    {
        volVars_ = other.volVars_;
        eIdx_ = other.eIdx_;
        return *this;
    };

public:
    void update(Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;

        volVars_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            auto fvGeometry = localView(problem.model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
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
        return volVars_[eIdx_][scvIdx];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        return volVars_[eIdx_][scvIdx];
    }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        return volVars_[eIdx_][scv.index()];
    }

    VolumeVariables& operator [](const SubControlVolume& scv)
    {
        return volVars_[eIdx_][scv.index()];
    }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element)
    { bindElement(element); }
    // function to prepare the vol vars within the element
    void bindElement(const Element& element)
    { eIdx_ = problem_().elementMapper().index(element); }

private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    IndexType eIdx_;
    std::vector<std::vector<VolumeVariables>> volVars_;
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
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    BoxVolumeVariablesVector& operator= (const BoxVolumeVariablesVector<TypeTag, false, false>& other) = default;

    // operator curVolVars = prevVolVars
    void operator= (const BoxVolumeVariablesVector<TypeTag, true, false>& other)
    {}

public:

    void update(Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;
    }

    // specialization for box models, simply forwards to the bindElement method
    void bind(const Element& element, const FVElementGeometry& fvGeometry)
    {
        bindElement(element, fvGeometry);
    }

    // specialization for box models
    void bindElement(const Element& element, const FVElementGeometry& fvGeometry)
    {
        // resize volume variables to the required size
        volVars_.resize(fvGeometry.numScv());

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& sol = problem_().model().curSol()[scv.dofIndex()];
            // TODO: INTERFACE SOLVER
            // problem_().model().boxInterfaceConditionSolver().updateScvVolVars(element, scv, sol);
            volVars_[scv.index()].update(sol, problem_(), element, scv);
        }
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    {
        return volVars_[scvIdx];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        return volVars_[scvIdx];
    }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        return volVars_[scv.index()];
    }

    VolumeVariables& operator [](const SubControlVolume& scv)
    {
        return volVars_[scv.index()];
    }

private:

    Problem& problem_() const
    { return *problemPtr_;}

    Problem* problemPtr_;
    std::vector<VolumeVariables> volVars_;
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
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    BoxVolumeVariablesVector& operator= (const BoxVolumeVariablesVector<TypeTag, /*isOldSol*/false, /*enableVolVarCaching*/false>& other)
    {
        volVars_.clear();
        problemPtr_ = other.problemPtr_;
        return *this;
    }

public:

    void update(const Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;
    }

    // bind an element
    void bindElement(const Element& element, const FVElementGeometry& fvGeometry)
    {
        // resize volume variables to the required size
        volVars_.resize(fvGeometry.numScv());

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& sol = problem_().model().prevSol()[scv.dofIndex()];
            //TODO: INTERFACE SOLVER?
            // problem_().model().boxInterfaceConditionSolver().updateScvVolVars(element, scv, sol);
            volVars_[scv.index()].update(sol, problem_(), element, scv);
        }
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    {
        return volVars_[scvIdx];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        return volVars_[scvIdx];
    }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        return volVars_[scv.index()];
    }

    VolumeVariables& operator [](const SubControlVolume& scv)
    {
        return volVars_[scv.index()];
    }

private:

    Problem& problem_() const
    { return *problemPtr_;}

    Problem* problemPtr_;

    std::vector<VolumeVariables> volVars_;
};

} // end namespace

#endif
