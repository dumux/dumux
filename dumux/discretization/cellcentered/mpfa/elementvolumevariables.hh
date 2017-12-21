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

#include <utility>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the local volume variables vector
 */
template<class TypeTag, bool enableGridVolVarsCache>
class CCMpfaElementVolumeVariables
{};

// specialization in case of storing the volume variables globally
template<class TypeTag>
class CCMpfaElementVolumeVariables<TypeTag, /*enableGridVolVarsCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Constructor
    CCMpfaElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return gridVolVars().volVars(scv.dofIndex()); }

    // operator for the access with an index
    // needed for cc methods for the access to the boundary volume variables
    const VolumeVariables& operator [](const IndexType scvIdx) const
    { return gridVolVars().volVars(scvIdx); }

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
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class CCMpfaElementVolumeVariables<TypeTag, /*enableGridVolVarsCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    //! Constructor
    CCMpfaElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    // Binding of an element, prepares the volume variables within the element stencil
    // called by the local jacobian to prepare element assembly
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        const auto& problem = gridVolVars().problem();
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        const auto& gridIvIndexSets = fvGridGeometry.gridInteractionVolumeIndexSets();

        // stencil information
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& assemblyMapI = fvGridGeometry.connectivityMap()[globalI];
        const auto numVolVars = assemblyMapI.size() + 1;

        // resize local containers to the required size (for internal elements)
        volumeVariables_.resize(numVolVars);
        volVarIndices_.resize(numVolVars);
        unsigned int localIdx = 0;

        // update the volume variables of the element at hand
        const auto& scvI = fvGeometry.scv(globalI);
        volumeVariables_[localIdx].update(ElementSolution(sol[globalI]),
                                          problem,
                                          element,
                                          scvI);
        volVarIndices_[localIdx] = scvI.dofIndex();
        ++localIdx;

        // Update the volume variables of the neighboring elements
        for (auto&& dataJ : assemblyMapI)
        {
            const auto& elementJ = fvGridGeometry.element(dataJ.globalJ);
            const auto& scvJ = fvGeometry.scv(dataJ.globalJ);
            volumeVariables_[localIdx].update(ElementSolution(sol[dataJ.globalJ]),
                                              problem,
                                              elementJ,
                                              scvJ);
            volVarIndices_[localIdx] = scvJ.dofIndex();
            ++localIdx;
        }

        // eventually prepare boundary volume variables
        const auto maxNumBoundaryVolVars = maxNumBoundaryVolVars_(fvGeometry);
        if (maxNumBoundaryVolVars > 0)
        {
            volumeVariables_.reserve(numVolVars+maxNumBoundaryVolVars);
            volVarIndices_.reserve(numVolVars+maxNumBoundaryVolVars);

            // treat the BCs inside the element
            for (const auto& scvf : scvfs(fvGeometry))
            {
                // if we are not on a boundary, skip to the next scvf
                if (!scvf.boundary())
                    continue;

                const auto bcTypes = problem.boundaryTypes(element, scvf);

                // on dirichlet boundaries use dirichlet values
                if (bcTypes.hasOnlyDirichlet())
                {
                    // boundary volume variables
                    VolumeVariables dirichletVolVars;
                    dirichletVolVars.update(ElementSolution(problem.dirichlet(element, scvf)),
                                            problem,
                                            element,
                                            scvI);

                    volumeVariables_.emplace_back(std::move(dirichletVolVars));
                    volVarIndices_.push_back(scvf.outsideScvIdx());
                }
                // use the inside volume variables for neumann boundaries
                else
                {
                    volumeVariables_.emplace_back(volumeVariables_[0]);
                    volVarIndices_.push_back(scvf.outsideScvIdx());
                }
            }

            // Update boundary volume variables in the neighbors
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (fvGridGeometry.vertexUsesPrimaryInteractionVolume(scvf.vertexIndex()))
                {
                    const auto& nodalIndexSet = gridIvIndexSets.primaryIndexSet(scvf).nodalIndexSet();

                    // if present, insert boundary vol vars
                    if (nodalIndexSet.numBoundaryScvfs() > 0)
                        addBoundaryVolVars_(problem, fvGeometry, nodalIndexSet);

                }
                else
                {
                    const auto& nodalIndexSet = gridIvIndexSets.secondaryIndexSet(scvf).nodalIndexSet();

                    // if present, insert boundary vol vars
                    if (nodalIndexSet.numBoundaryScvfs() > 0)
                        addBoundaryVolVars_(problem, fvGeometry, nodalIndexSet);
                }
            }
        }

        // //! TODO Check if user added additional DOF dependencies, i.e. the residual of DOF globalI depends
        // //! on additional DOFs not included in the discretization schemes' occupation pattern
        // const auto& additionalDofDependencies = problem.getAdditionalDofDependencies(globalI);
        // if (!additionalDofDependencies.empty())
        // {
        //     volumeVariables_.reserve(volumeVariables_.size() + additionalDofDependencies.size());
        //     volVarIndices_.reserve(volVarIndices_.size() + additionalDofDependencies.size());
        //     for (auto globalJ : additionalDofDependencies)
        //     {
        //         const auto& elementJ = fvGridGeometry.element(globalJ);
        //         const auto& scvJ = fvGeometry.scv(globalJ);

        //         VolumeVariables additionalVolVars;
        //         additionalVolVars.update(ElementSolution(elementJ, sol, fvGridGeometry),
        //                                  problem,
        //                                  elementJ,
        //                                  scvJ);

        //         volumeVariables_.emplace_back(std::move(additionalVolVars));
        //         volVarIndices_.push_back(globalJ);
        //     }
        // }
    }

    // Binding of an element, prepares only the volume variables of the element
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        auto eIdx = fvGeometry.fvGridGeometry().elementMapper().index(element);
        volumeVariables_.resize(1);
        volVarIndices_.resize(1);

        // update the volume variables of the element
        const auto& scv = fvGeometry.scv(eIdx);
        volumeVariables_[0].update(ElementSolution(sol[scv.dofIndex()]),
                                   gridVolVars().problem(),
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
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;

    //! Computes how many boundary vol vars come into play for flux calculations
    //! on this element. This number is probably almost always higher than the actually
    //! needed number of volume variables. However, memory is not an issue for the global
    //! caching being deactivated and we want to make sure we reserve enough memory here.
    // TODO: What about non-symmetric schemes? Is there a better way for estimating this?
    std::size_t maxNumBoundaryVolVars_(const FVElementGeometry& fvGeometry)
    {
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        const auto& gridIvIndexSets = fvGridGeometry.gridInteractionVolumeIndexSets();

        std::size_t numBoundaryVolVars = 0;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (fvGridGeometry.vertexUsesPrimaryInteractionVolume(scvf.vertexIndex()))
                numBoundaryVolVars += gridIvIndexSets.primaryIndexSet(scvf).nodalIndexSet().numBoundaryScvfs();
            else
                numBoundaryVolVars += gridIvIndexSets.secondaryIndexSet(scvf).nodalIndexSet().numBoundaryScvfs();
        }

        return numBoundaryVolVars;
    }

    //! adds the volume variables from a given nodal index set to the local containers
    template<class NodalIndexSet>
    void addBoundaryVolVars_(const Problem& problem, const FVElementGeometry& fvGeometry, const NodalIndexSet& nodalIndexSet)
    {
        // check each scvf in the index set for boundary presence
        for (auto scvfIdx : nodalIndexSet.globalScvfIndices())
        {
            const auto& ivScvf = fvGeometry.scvf(scvfIdx);

            // only proceed for scvfs on the boundary and not in the inside element
            if (!ivScvf.boundary() || ivScvf.insideScvIdx() == volVarIndices_[0])
                continue;

            const auto insideScvIdx = ivScvf.insideScvIdx();
            const auto insideElement = fvGeometry.fvGridGeometry().element(insideScvIdx);
            const auto bcTypes = problem.boundaryTypes(insideElement, ivScvf);

            // on dirichlet boundaries use dirichlet values
            if (bcTypes.hasOnlyDirichlet())
            {
                // boundary volume variables
                VolumeVariables dirichletVolVars;
                const auto& ivScv = fvGeometry.scv(insideScvIdx);
                dirichletVolVars.update(ElementSolution(problem.dirichlet(insideElement, ivScvf)),
                                        problem,
                                        insideElement,
                                        ivScv);

                volumeVariables_.emplace_back(std::move(dirichletVolVars));
                volVarIndices_.push_back(ivScvf.outsideScvIdx());
            }
            // use the inside volume variables for neumann boundaries
            else
            {
                volumeVariables_.emplace_back((*this)[insideScvIdx]);
                volVarIndices_.push_back(ivScvf.outsideScvIdx());
            }
        }
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
