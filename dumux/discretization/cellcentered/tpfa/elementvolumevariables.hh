// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup CCTpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered tpfa models
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CCTPFA_ELEMENT_VOLUMEVARIABLES_HH

#include <algorithm>
#include <type_traits>
#include <vector>

#include <dumux/discretization/cellcentered/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered tpfa models
 * \note The class is specilized for versions with and without caching
 * \tparam GVV the grid volume variables type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GVV, bool cachingEnabled>
class CCTpfaElementVolumeVariables
{};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered tpfa models with caching
 * \note the volume variables are stored for the whole grid view in the corresponding GridVolumeVariables class
 */
template<class GVV>
class CCTpfaElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    CCTpfaElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars)
    , numScv_(gridVolVars.problem().gridGeometry().numScv())
    {}

    //! operator for the access with an scv
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        if (scv.dofIndex() < numScv_)
            return gridVolVars().volVars(scv.dofIndex());
        else
            return boundaryVolumeVariables_[getLocalIdx_(scv.dofIndex())];
    }

    //! operator for the access with an index
    const VolumeVariables& operator [](const std::size_t scvIdx) const
    {
        if (scvIdx < numScv_)
            return gridVolVars().volVars(scvIdx);
        else
            return boundaryVolumeVariables_[getLocalIdx_(scvIdx)];
    }

    //! precompute all boundary volume variables in a stencil of an element, the remaining ones are cached
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        if (!fvGeometry.hasBoundaryScvf())
            return;

        clear_();
        boundaryVolVarIndices_.reserve(fvGeometry.numScvf());
        boundaryVolumeVariables_.reserve(fvGeometry.numScvf());

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary())
                continue;

            // check if boundary is a pure dirichlet boundary
            const auto& problem = gridVolVars().problem();
            const auto bcTypes = problem.boundaryTypes(element, scvf);
            if (bcTypes.hasOnlyDirichlet())
            {
                const auto dirichletPriVars = elementSolution<FVElementGeometry>(problem.dirichlet(element, scvf));
                auto&& scvI = fvGeometry.scv(scvf.insideScvIdx());

                VolumeVariables volVars;
                volVars.update(dirichletPriVars,
                               problem,
                               element,
                               scvI);

                boundaryVolumeVariables_.emplace_back(std::move(volVars));
                boundaryVolVarIndices_.push_back(scvf.outsideScvIdx());
            }
        }
    }

    //! precompute the volume variables of an element - do nothing: volVars are cached
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    //! Clear all local storage
    void clear_()
    {
        boundaryVolVarIndices_.clear();
        boundaryVolumeVariables_.clear();
    }

    const GridVolumeVariables* gridVolVarsPtr_;

    //! map a global scv index to the local storage index
    int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(boundaryVolVarIndices_.begin(), boundaryVolVarIndices_.end(), volVarIdx);
        assert(it != boundaryVolVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(boundaryVolVarIndices_.begin(), it);
    }

    std::vector<std::size_t> boundaryVolVarIndices_;
    std::vector<VolumeVariables> boundaryVolumeVariables_;
    const std::size_t numScv_;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered tpfa models with caching
 */
template<class GVV>
class CCTpfaElementVolumeVariables<GVV, /*cachingEnabled*/false>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    CCTpfaElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    //! Prepares the volume variables within the element stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        clear_();

        const auto& problem = gridVolVars().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& connectivityMapI = gridGeometry.connectivityMap()[globalI];
        const auto numDofs = connectivityMapI.size() + 1;

        // resize local containers to the required size (for internal elements)
        volumeVariables_.resize(numDofs);
        volVarIndices_.resize(numDofs);
        int localIdx = 0;

        // update the volume variables of the element at hand
        auto&& scvI = fvGeometry.scv(globalI);
        volumeVariables_[localIdx].update(elementSolution(element, sol, gridGeometry),
                                          problem,
                                          element,
                                          scvI);
        volVarIndices_[localIdx] = scvI.dofIndex();
        ++localIdx;

        // Update the volume variables of the neighboring elements
        for (const auto& dataJ : connectivityMapI)
        {
            const auto& elementJ = gridGeometry.element(dataJ.globalJ);
            auto&& scvJ = fvGeometry.scv(dataJ.globalJ);
            volumeVariables_[localIdx].update(elementSolution(elementJ, sol, gridGeometry),
                                              problem,
                                              elementJ,
                                              scvJ);
            volVarIndices_[localIdx] = scvJ.dofIndex();
            ++localIdx;
        }

        if (fvGeometry.hasBoundaryScvf())
        {
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
                    const auto dirichletPriVars = elementSolution<FVElementGeometry>(problem.dirichlet(element, scvf));

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
        }

        //! Check if user added additional DOF dependencies, i.e. the residual of DOF globalI depends
        //! on additional DOFs not included in the discretization schemes' occupation pattern
        // const auto& additionalDofDependencies = problem.getAdditionalDofDependencies(globalI);
        // if (!additionalDofDependencies.empty())
        // {
        //     volumeVariables_.resize(volumeVariables_.size() + additionalDofDependencies.size());
        //     volVarIndices_.resize(volVarIndices_.size() + additionalDofDependencies.size());
        //     for (auto globalJ : additionalDofDependencies)
        //     {
        //         const auto& elementJ = gridGeometry.element(globalJ);
        //         auto&& scvJ = fvGeometry.scv(globalJ);

        //         volumeVariables_[localIdx].update(elementSolution(elementJ, sol, gridGeometry),
        //                                           problem,
        //                                           elementJ,
        //                                           scvJ);
        //         volVarIndices_[localIdx] = scvJ.dofIndex();
        //         ++localIdx;
        //     }
        // }
    }

    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        clear_();

        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        volumeVariables_.resize(1);
        volVarIndices_.resize(1);

        // update the volume variables of the element
        auto&& scv = fvGeometry.scv(eIdx);
        volumeVariables_[0].update(elementSolution(element, sol, fvGeometry.gridGeometry()),
                                   gridVolVars().problem(),
                                   element,
                                   scv);
        volVarIndices_[0] = scv.dofIndex();
    }

    //! access operator with scv
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    //! access operator with scv
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    //! access operator with scv index
    const VolumeVariables& operator [](std::size_t scvIdx) const
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! access operator with scv index
    VolumeVariables& operator [](std::size_t scvIdx)
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    //! Clear all local storage
    void clear_()
    {
        volVarIndices_.clear();
        volumeVariables_.clear();
    }

    const GridVolumeVariables* gridVolVarsPtr_;

    //! map a global scv index to the local storage index
    int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(volVarIndices_.begin(), volVarIndices_.end(), volVarIdx);
        assert(it != volVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(volVarIndices_.begin(), it);
    }

    std::vector<std::size_t> volVarIndices_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace Dumux

#endif
