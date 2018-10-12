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
 * \ingroup GodunovDiscretization
 * \brief The local (stencil) volume variables class for cell centered godunov models
 */
#ifndef DUMUX_DISCRETIZATION_Godunov_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_Godunov_ELEMENT_VOLUMEVARIABLES_HH

#include <algorithm>
#include <type_traits>
#include <vector>

#include <dumux/discretization/cellcentered/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup GodunovDiscretization
 * \brief The local (stencil) volume variables class for cell centered godunov models
 * \note The class is specilized for versions with and without caching
 * \tparam GVV the grid volume variables type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GVV, bool cachingEnabled>
class GodunovElementVolumeVariables
{};

/*!
 * \ingroup GodunovDiscretization
 * \brief The local (stencil) volume variables class for cell centered godunov models with caching
 * \note the volume variables are stored for the whole grid view in the corresponding GridVolumeVariables class
 */
template<class GVV>
class GodunovElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    GodunovElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    //! operator for the access with an scv
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return gridVolVars().volVars(scv.dofIndex()); }

    //! operator for the access with an index
    const VolumeVariables& operator [](const std::size_t scvIdx) const
    { return gridVolVars().volVars(scvIdx); }

    //! precompute all volume variables in a stencil of an element - do nothing volVars: are cached
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {}

    //! precompute the volume variables of an element - do nothing: volVars are cached
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
};

/*!
 * \ingroup GodunovDiscretization
 * \brief The local (stencil) volume variables class for cell centered godunov models with caching
 */
template<class GVV>
class GodunovElementVolumeVariables<GVV, /*cachingEnabled*/false>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    GodunovElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    //! Prepares the volume variables within the element stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        clear();

        const auto& problem = gridVolVars().problem();
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        const auto globalI = fvGridGeometry.elementMapper().index(element);
        const auto& connectivityMapI = fvGridGeometry.connectivityMap()[globalI];
        const auto numDofs = connectivityMapI.size() + 1;

        // resize local containers to the required size (for internal elements)
        volumeVariables_.resize(numDofs);
        volVarIndices_.resize(numDofs);
        int localIdx = 0;

        // update the volume variables of the element at hand
        auto&& scvI = fvGeometry.scv(globalI);
        volumeVariables_[localIdx].update(elementSolution(element, sol, fvGridGeometry),
                                          problem,
                                          element,
                                          scvI);
        volVarIndices_[localIdx] = scvI.dofIndex();
        ++localIdx;

        // Update the volume variables of the neighboring elements
        for (const auto& dataJ : connectivityMapI)
        {
            const auto& elementJ = fvGridGeometry.element(dataJ.globalJ);
            auto&& scvJ = fvGeometry.scv(dataJ.globalJ);
            volumeVariables_[localIdx].update(elementSolution(elementJ, sol, fvGridGeometry),
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

        //! Check if user added additional DOF dependencies, i.e. the residual of DOF globalI depends
        //! on additional DOFs not included in the discretization schemes' occupation pattern
        // const auto& additionalDofDependencies = problem.getAdditionalDofDependencies(globalI);
        // if (!additionalDofDependencies.empty())
        // {
        //     volumeVariables_.resize(volumeVariables_.size() + additionalDofDependencies.size());
        //     volVarIndices_.resize(volVarIndices_.size() + additionalDofDependencies.size());
        //     for (auto globalJ : additionalDofDependencies)
        //     {
        //         const auto& elementJ = fvGridGeometry.element(globalJ);
        //         auto&& scvJ = fvGeometry.scv(globalJ);

        //         volumeVariables_[localIdx].update(elementSolution(elementJ, sol, fvGridGeometry),
        //                                           problem,
        //                                           elementJ,
        //                                           scvJ);
        //         volVarIndices_[localIdx] = scvJ.dofIndex();
        //         ++localIdx;
        //     }
        // }
    }

    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        clear();

        const auto eIdx = fvGeometry.fvGridGeometry().elementMapper().index(element);
        volumeVariables_.resize(1);
        volVarIndices_.resize(1);

        // update the volume variables of the element
        auto&& scv = fvGeometry.scv(eIdx);
        volumeVariables_[0].update(elementSolution(element, sol, fvGeometry.fvGridGeometry()),
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

    //! Clear all local storage
    void clear()
    {
        volVarIndices_.clear();
        volumeVariables_.clear();
    }
private:
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
