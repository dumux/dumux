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
 * \brief The local (stencil) volume variables class for cell centered models
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CCTPFA_ELEMENT_VOLUMEVARIABLES_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool enableGlobalVolVarsCache>
class CCTpfaElementVolumeVariables
{};

// specialization in case of storing the volume variables globally
template<class TypeTag>
class CCTpfaElementVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Constructor
    CCTpfaElementVolumeVariables(const GlobalVolumeVariables& globalVolVars)
    : globalVolVarsPtr_(&globalVolVars) {}

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return globalVolVars().volVars(scv.dofIndex()); }

    // operator for the access with an index
    // needed for cc methods for the access to the boundary volume variables
    const VolumeVariables& operator [](const IndexType scvIdx) const
    { return globalVolVars().volVars(scvIdx); }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {}

    // function to prepare the vol vars within the element
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}

    //! The global volume variables object we are a restriction of
    const GlobalVolumeVariables& globalVolVars() const
    { return *globalVolVarsPtr_; }

private:
    const GlobalVolumeVariables* globalVolVarsPtr_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class CCTpfaElementVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    //! Constructor
    CCTpfaElementVolumeVariables(const GlobalVolumeVariables& globalVolVars)
    : globalVolVarsPtr_(&globalVolVars) {}

    // Binding of an element, prepares the volume variables within the element stencil
    // called by the local jacobian to prepare element assembly
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        clear();

        const auto& problem = globalVolVars().problem_();
        const auto globalI = problem.elementMapper().index(element);
        const auto& assemblyMapI = problem.model().localJacobian().assemblyMap()[globalI];
        const auto numDofs = assemblyMapI.size() + 1;

        // resize local containers to the required size (for internal elements)
        volumeVariables_.resize(numDofs);
        volVarIndices_.resize(numDofs);
        int localIdx = 0;

        // update the volume variables of the element at hand
        auto&& scvI = fvGeometry.scv(globalI);
        volumeVariables_[localIdx].update(problem.model().elementSolution(element, sol),
                                          problem,
                                          element,
                                          scvI);
        volVarIndices_[localIdx] = scvI.dofIndex();
        ++localIdx;

        // Update the volume variables of the neighboring elements
        for (const auto& dataJ : assemblyMapI)
        {
            const auto& elementJ = fvGeometry.globalFvGeometry().element(dataJ.globalJ);
            auto&& scvJ = fvGeometry.scv(dataJ.globalJ);
            volumeVariables_[localIdx].update(problem.model().elementSolution(elementJ, sol),
                                              problem,
                                              elementJ,
                                              scvJ);
            volVarIndices_[localIdx] = scvJ.dofIndex();
            ++localIdx;
        }

        // Update boundary volume variables
        for (auto&& scvf : scvfs(fvGeometry))
        {
            // if we are not on a boundary, skip to the next scvf
            if (!scvf.boundary())
                continue;

            // check if boundary is a pure dirichlet boundary
            const auto bcTypes = problem.boundaryTypes(element, scvf);
            if (bcTypes.hasOnlyDirichlet())
            {
                const ElementSolution dirichletPriVars({problem.dirichlet(element, scvf)});

                volumeVariables_.resize(localIdx+1);
                volVarIndices_.resize(localIdx+1);
                volumeVariables_[localIdx].update(dirichletPriVars,
                                                  problem,
                                                  element,
                                                  scvI);
                volVarIndices_[localIdx] = scvf.outsideScvIdx();
                ++localIdx;
            }
        }

        //! Check if user added additional DOF dependencies, i.e. the residual of DOF globalI depends
        //! on additional DOFs not included in the discretization schemes' occupation pattern
        const auto& additionalDofDependencies = problem.getAdditionalDofDependencies(globalI);
        if (!additionalDofDependencies.empty())
        {
            volumeVariables_.resize(volumeVariables_.size() + additionalDofDependencies.size());
            volVarIndices_.resize(volVarIndices_.size() + additionalDofDependencies.size());
            for (auto globalJ : additionalDofDependencies)
            {
                const auto& elementJ = fvGeometry.globalFvGeometry().element(globalJ);
                auto&& scvJ = fvGeometry.scv(globalJ);

                volumeVariables_[localIdx].update(problem.model().elementSolution(elementJ, sol),
                                                  problem,
                                                  elementJ,
                                                  scvJ);
                volVarIndices_[localIdx] = scvJ.dofIndex();
                ++localIdx;
            }
        }
    }

    // Binding of an element, prepares only the volume variables of the element
    // specialization for cc models
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        clear();

        auto eIdx = globalVolVars().problem_().elementMapper().index(element);
        volumeVariables_.resize(1);
        volVarIndices_.resize(1);

        // update the volume variables of the element
        auto&& scv = fvGeometry.scv(eIdx);
        volumeVariables_[0].update(globalVolVars().problem_().model().elementSolution(element, sol),
                                   globalVolVars().problem_(),
                                   element,
                                   scv);
        volVarIndices_[0] = scv.dofIndex();
    }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    const VolumeVariables& operator [](IndexType scvIdx) const
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    VolumeVariables& operator [](IndexType scvIdx)
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! The global volume variables object we are a restriction of
    const GlobalVolumeVariables& globalVolVars() const
    { return *globalVolVarsPtr_; }

    //! Clear all local storage
    void clear()
    {
        volVarIndices_.clear();
        volumeVariables_.clear();
    }
private:
    const GlobalVolumeVariables* globalVolVarsPtr_;

    int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(volVarIndices_.begin(), volVarIndices_.end(), volVarIdx);
        assert(it != volVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(volVarIndices_.begin(), it);
    }

    std::vector<IndexType> volVarIndices_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace

#endif
