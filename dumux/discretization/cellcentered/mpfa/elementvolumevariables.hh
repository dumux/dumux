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
 * \brief The local (stencil) volume variables class for cell centered mpfa models
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_VOLUMEVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include "facetypes.hh"

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the local volume variables vector
 */
template<class TypeTag, bool enableGlobalVolVarsCache>
class CCMpfaElementVolumeVariables
{};

// specialization in case of storing the volume variables globally
template<class TypeTag>
class CCMpfaElementVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Constructor
    CCMpfaElementVolumeVariables(const GlobalVolumeVariables& globalVolVars)
    : globalVolVarsPtr_(&globalVolVars) {}

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return globalVolVars().volVars(scv.index()); }

    // operator for the access with an index
    // needed for cc methods for the access to the boundary volume variables
    const VolumeVariables& operator [](const IndexType scvIdx) const
    { return globalVolVars().volVars(scvIdx); }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) {}

    // function to prepare the vol vars within the element
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) {}

    //! The global volume variables object we are a restriction of
    const GlobalVolumeVariables& globalVolVars() const
    { return *globalVolVarsPtr_; }

private:
    const GlobalVolumeVariables* globalVolVarsPtr_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class CCMpfaElementVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    //! Constructor
    CCMpfaElementVolumeVariables(const GlobalVolumeVariables& globalVolVars)
    : globalVolVarsPtr_(&globalVolVars) {}

    // Binding of an element, prepares the volume variables within the element stencil
    // called by the local jacobian to prepare element assembly
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        const auto& problem = globalVolVars().problem_();
        const auto& globalFvGeometry = problem.model().globalFvGeometry();

        // stencil information
        const auto globalI = problem.elementMapper().index(element);
        const auto& assemblyMapI = problem.model().localJacobian().assemblyMap()[globalI];
        const auto numDofs = assemblyMapI.size() + 1;

        // resize local containers to the required size (for internal elements)
        volumeVariables_.resize(numDofs);
        volVarIndices_.resize(numDofs);
        int localIdx = 0;

        // update the volume variables of the element at hand
        auto eIdx = problem.elementMapper().index(element);
        auto&& scvI = fvGeometry.scv(eIdx);
        volumeVariables_[localIdx].update(problem.model().elementSolution(element, sol),
                                          problem,
                                          element,
                                          scvI);
        volVarIndices_[localIdx] = scvI.index();
        ++localIdx;

        // Update the volume variables of the neighboring elements
        for (auto&& dataJ : assemblyMapI)
        {
            const auto& elementJ = globalFvGeometry.element(dataJ.globalJ);
            auto&& scvJ = fvGeometry.scv(dataJ.globalJ);
            volumeVariables_[localIdx].update(problem.model().elementSolution(elementJ, sol),
                                              problem,
                                              elementJ,
                                              scvJ);
            volVarIndices_[localIdx] = scvJ.index();
            ++localIdx;
        }

        // eventually prepare boundary volume variables
        auto estimate = boundaryVolVarsEstimate_(element, fvGeometry);
        if (estimate > 0)
        {
            volumeVariables_.reserve(numDofs+estimate);
            volVarIndices_.reserve(numDofs+estimate);

            // treat the BCs inside the element
            if (element.hasBoundaryIntersections())
            {
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    // if we are not on a boundary, skip to the next scvf
                    if (!scvf.boundary())
                        continue;

                    // on dirichlet boundaries use dirichlet values
                    if (MpfaHelper::getMpfaFaceType(problem, element, scvf) == MpfaFaceTypes::dirichlet)
                    {
                        // boundary volume variables
                        VolumeVariables dirichletVolVars;
                        dirichletVolVars.update(ElementSolution({problem.dirichlet(element, scvf)}),
                                                problem,
                                                element,
                                                scvI);

                        volumeVariables_.emplace_back(std::move(dirichletVolVars));
                        volVarIndices_.push_back(scvf.outsideScvIdx());
                    }
                    // use the inside volume variables for neumann boundaries
                    else if (!useTpfaBoundary)
                    {
                        volumeVariables_.emplace_back(volumeVariables_[0]);
                        volVarIndices_.push_back(scvf.outsideScvIdx());
                    }
                }
            }

            // Update boundary volume variables in the neighbors
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // skip the rest if the scvf does not touch a boundary
                if (!globalFvGeometry.scvfTouchesBoundary(scvf))
                    continue;

                // loop over all the scvfs in the interaction region
                const auto& ivSeed = globalFvGeometry.boundaryInteractionVolumeSeed(scvf);
                for (auto scvfIdx : ivSeed.globalScvfIndices())
                {
                    auto&& ivScvf = fvGeometry.scvf(scvfIdx);
                    // only proceed for scvfs on the boundary and not in the inside element
                    if (!ivScvf.boundary() || ivScvf.insideScvIdx() == eIdx)
                        continue;

                    auto insideScvIdx = ivScvf.insideScvIdx();
                    auto insideElement = globalFvGeometry.element(insideScvIdx);

                    // on dirichlet boundaries use dirichlet values
                    if (MpfaHelper::getMpfaFaceType(problem, insideElement, ivScvf) == MpfaFaceTypes::dirichlet)
                    {
                        // boundary volume variables
                        VolumeVariables dirichletVolVars;
                        auto&& ivScv = fvGeometry.scv(insideScvIdx);
                        dirichletVolVars.update(ElementSolution({problem.dirichlet(insideElement, ivScvf)}),
                                                problem,
                                                insideElement,
                                                ivScv);

                        volumeVariables_.emplace_back(std::move(dirichletVolVars));
                        volVarIndices_.push_back(ivScvf.outsideScvIdx());
                    }
                    // use the inside volume variables for neumann boundaries
                    else if (!useTpfaBoundary)
                    {
                        volumeVariables_.emplace_back((*this)[insideScvIdx]);
                        volVarIndices_.push_back(ivScvf.outsideScvIdx());
                    }
                }
            }
        }
    }

    // Binding of an element, prepares only the volume variables of the element
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        auto eIdx = globalVolVars().problem_().elementMapper().index(element);
        volumeVariables_.resize(1);
        volVarIndices_.resize(1);

        // update the volume variables of the element
        auto&& scv = fvGeometry.scv(eIdx);
        volumeVariables_[0].update(globalVolVars().problem_().model().elementSolution(element, sol),
                                   globalVolVars().problem_(),
                                   element,
                                   scv);
        volVarIndices_[0] = scv.index();
    }

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[getLocalIdx_(scv.index())]; }

    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[getLocalIdx_(scv.index())]; }

    const VolumeVariables& operator [](IndexType scvIdx) const
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    VolumeVariables& operator [](IndexType scvIdx)
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! The global volume variables object we are a restriction of
    const GlobalVolumeVariables& globalVolVars() const
    { return *globalVolVarsPtr_; }

private:
    const GlobalVolumeVariables* globalVolVarsPtr_;

    //! checks whether an scvf touches the boundary and returns an estimate of how many
    //! boundary vol vars will be necessary. In 2d, this is the sum of faces touching
    //! the boundary, which should be correct. In 3d, we count each face double - probably
    //! too much for hexahedrons but might be even too little for simplices.
    int boundaryVolVarsEstimate_(const Element& element,
                                 const FVElementGeometry& fvGeometry)
    {
        int bVolVarEstimate = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            bool boundary = scvf.boundary();
            if (boundary || (!boundary && fvGeometry.globalFvGeometry().scvfTouchesBoundary(scvf)))
                bVolVarEstimate += dim-1;
        }

        return bVolVarEstimate;
    }

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
