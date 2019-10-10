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
 * \ingroup CCMpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered mpfa models
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_VOLUMEVARIABLES_HH

#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>

#include <dumux/discretization/cellcentered/elementsolution.hh>

namespace Dumux {
namespace CCMpfa {

    /*!
     * \ingroup CCMpfaDiscretization
     * \brief Computes how many boundary vol vars come into play for flux calculations
     *        on an element (for a given element finite volume geometry). This number here
     *        is probably always higher than the actually needed number of volume variables.
     *        However, we want to make sure it is high enough so that enough memory is reserved
     *        in the element volume variables below.
     * \todo TODO What about non-symmetric schemes? Is there a better way for estimating this?
     *
     * \param fvGeometry the element finite volume geometry
     */
    template<class FVElementGeometry>
    std::size_t maxNumBoundaryVolVars(const FVElementGeometry& fvGeometry)
    {
        const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto& gridIvIndexSets = gridGeometry.gridInteractionVolumeIndexSets();

        std::size_t numBoundaryVolVars = 0;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!gridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
                numBoundaryVolVars += gridIvIndexSets.primaryIndexSet(scvf).nodalIndexSet().numBoundaryScvfs();
            else
                numBoundaryVolVars += gridIvIndexSets.secondaryIndexSet(scvf).nodalIndexSet().numBoundaryScvfs();
        }

        return numBoundaryVolVars;
    }

    /*!
     * \ingroup CCMpfaDiscretization
     * \brief Adds the boundary volume variables found within a given nodal index set
     *        into the provided containers and stores the indices associated with them.
     * \note  It only adds those boundary vol vars that do not live on scvfs that are
     *        inside the bound element. These have to be added separately
     *
     * \param volVars       The container where the volume variables are stored
     * \param volVarIndices The container where the volume variable indices are stored
     * \param problem       The problem containing the Dirichlet boundary conditions
     * \param element       The element to which the finite volume geometry is bound
     * \param fvGeometry    The element finite volume geometry
     * \param nodalIndexSet The dual grid index set around a node
     */
    template<class VolumeVariables, class IndexType, class Problem, class FVElemGeom, class NodalIndexSet>
    void addBoundaryVolVarsAtNode(std::vector<VolumeVariables>& volVars,
                                  std::vector<IndexType>& volVarIndices,
                                  const Problem& problem,
                                  const typename FVElemGeom::GridGeometry::GridView::template Codim<0>::Entity& element,
                                  const FVElemGeom& fvGeometry,
                                  const NodalIndexSet& nodalIndexSet)
    {
        if (nodalIndexSet.numBoundaryScvfs() == 0)
            return;

        // index of the element the fvGeometry was bound to
        const auto boundElemIdx = fvGeometry.gridGeometry().elementMapper().index(element);

        // check each scvf in the index set for boundary presence
        for (auto scvfIdx : nodalIndexSet.gridScvfIndices())
        {
            const auto& ivScvf = fvGeometry.scvf(scvfIdx);

            // only proceed for scvfs on the boundary and not in the bound element
            if (!ivScvf.boundary() || ivScvf.insideScvIdx() == boundElemIdx)
                continue;

            const auto insideScvIdx = ivScvf.insideScvIdx();
            const auto insideElement = fvGeometry.gridGeometry().element(insideScvIdx);
            const auto bcTypes = problem.boundaryTypes(insideElement, ivScvf);

            // Only proceed on dirichlet boundaries. On Neumann
            // boundaries the "outside" vol vars cannot be properly defined.
            if (bcTypes.hasOnlyDirichlet())
            {
                VolumeVariables dirichletVolVars;
                dirichletVolVars.update(elementSolution<FVElemGeom>(problem.dirichlet(insideElement, ivScvf)),
                                        problem,
                                        insideElement,
                                        fvGeometry.scv(insideScvIdx));

                volVars.emplace_back(std::move(dirichletVolVars));
                volVarIndices.push_back(ivScvf.outsideScvIdx());
            }
        }
    }

    /*!
     * \ingroup CCMpfaDiscretization
     * \brief Adds the boundary volume variables found within the stencil to the
     *        provided containers and stores the indices associated with them.
     *
     * \param volVars       The container where the volume variables are stored
     * \param volVarIndices The container where the volume variable indices are stored
     * \param problem       The problem containing the Dirichlet boundary conditions
     * \param element        The element to which the finite volume geometry was bound
     * \param fvGeometry    The element finite volume geometry
     */
    template<class VolumeVariables, class IndexType, class Problem, class FVElemGeom>
    void addBoundaryVolVars(std::vector<VolumeVariables>& volVars,
                            std::vector<IndexType>& volVarIndices,
                            const Problem& problem,
                            const typename FVElemGeom::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const FVElemGeom& fvGeometry)
    {
        const auto& gridGeometry = fvGeometry.gridGeometry();

        // treat the BCs inside the element
        if (fvGeometry.hasBoundaryScvf())
        {
            const auto boundElemIdx = gridGeometry.elementMapper().index(element);
            const auto& scvI = fvGeometry.scv(boundElemIdx);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (!scvf.boundary())
                    continue;

                // Only proceed on dirichlet boundaries. On Neumann
                // boundaries the "outside" vol vars cannot be properly defined.
                if (problem.boundaryTypes(element, scvf).hasOnlyDirichlet())
                {
                    VolumeVariables dirichletVolVars;
                    dirichletVolVars.update(elementSolution<FVElemGeom>(problem.dirichlet(element, scvf)),
                                            problem,
                                            element,
                                            scvI);

                    volVars.emplace_back(std::move(dirichletVolVars));
                    volVarIndices.push_back(scvf.outsideScvIdx());
                }
            }
        }

        // Update boundary volume variables in the neighbors
        const auto& gridIvIndexSets = gridGeometry.gridInteractionVolumeIndexSets();
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!gridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
                addBoundaryVolVarsAtNode( volVars, volVarIndices, problem, element, fvGeometry,
                                          gridIvIndexSets.primaryIndexSet(scvf).nodalIndexSet() );
            else
                addBoundaryVolVarsAtNode( volVars, volVarIndices, problem, element, fvGeometry,
                                          gridIvIndexSets.secondaryIndexSet(scvf).nodalIndexSet() );
        }
    }
} // end namespace CCMpfa

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered mpfa models
 * \note The class is specilized for versions with and without caching
 * \tparam GVV the grid volume variables type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GVV, bool cachingEnabled>
class CCMpfaElementVolumeVariables;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered mpfa models with caching
 * \note the volume variables are stored for the whole grid view in the corresponding GridVolumeVariables class
 */
template<class GVV>
class CCMpfaElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    CCMpfaElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars)
    , numScv_(gridVolVars.problem().gridGeometry().numScv())
    {}

    //! operator for the access with an scv
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        return scv.dofIndex() < numScv_ ? gridVolVars().volVars(scv.dofIndex())
                                        : boundaryVolVars_[getLocalIdx_(scv.dofIndex())];
    }

    //! operator for the access with an index
    const VolumeVariables& operator [](const std::size_t scvIdx) const
    {
        return scvIdx < numScv_ ? gridVolVars().volVars(scvIdx)
                                : boundaryVolVars_[getLocalIdx_(scvIdx)];
    }

    //! precompute all volume variables in a stencil of an element - bind Dirichlet vol vars in the stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        clear();

        // maybe prepare boundary volume variables
        const auto maxNumBoundaryVolVars = CCMpfa::maxNumBoundaryVolVars(fvGeometry);
        if (maxNumBoundaryVolVars > 0)
        {
            boundaryVolVars_.reserve(maxNumBoundaryVolVars);
            boundaryVolVarIndices_.reserve(maxNumBoundaryVolVars);
            CCMpfa::addBoundaryVolVars(boundaryVolVars_, boundaryVolVarIndices_, gridVolVars().problem(), element, fvGeometry);
        }
    }

    //! precompute the volume variables of an element - do nothing: volVars are cached
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}

    //! Clear all local storage
    void clear()
    {
        boundaryVolVarIndices_.clear();
        boundaryVolVars_.clear();
    }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    //! map a global scv index to the local storage index
    int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(boundaryVolVarIndices_.begin(), boundaryVolVarIndices_.end(), volVarIdx);
        assert(it != boundaryVolVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(boundaryVolVarIndices_.begin(), it);
    }

    const GridVolumeVariables* gridVolVarsPtr_;

    std::size_t numScv_;
    std::vector<std::size_t> boundaryVolVarIndices_;
    std::vector<VolumeVariables> boundaryVolVars_;
};


/*!
 * \ingroup CCMpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered tpfa models with caching
 */
template<class GVV>
class CCMpfaElementVolumeVariables<GVV, /*cachingEnabled*/false>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    CCMpfaElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    //! Prepares the volume variables within the element stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        clear();

        const auto& problem = gridVolVars().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();

        // stencil information
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& assemblyMapI = gridGeometry.connectivityMap()[globalI];
        const auto numVolVars = assemblyMapI.size() + 1;

        // resize local containers to the required size (for internal elements)
        const auto maxNumBoundaryVolVars = CCMpfa::maxNumBoundaryVolVars(fvGeometry);
        volumeVariables_.reserve(numVolVars+maxNumBoundaryVolVars);
        volVarIndices_.reserve(numVolVars+maxNumBoundaryVolVars);

        VolumeVariables volVars;
        const auto& scvI = fvGeometry.scv(globalI);
        volVars.update(elementSolution(element, sol, gridGeometry),
                       problem,
                       element,
                       scvI);

        volVarIndices_.push_back(scvI.dofIndex());
        volumeVariables_.emplace_back(std::move(volVars));

        // Update the volume variables of the neighboring elements
        for (auto&& dataJ : assemblyMapI)
        {
            const auto& elementJ = gridGeometry.element(dataJ.globalJ);
            const auto& scvJ = fvGeometry.scv(dataJ.globalJ);
            VolumeVariables volVarsJ;
            volVarsJ.update(elementSolution(elementJ, sol, gridGeometry),
                            problem,
                            elementJ,
                            scvJ);

            volVarIndices_.push_back(scvJ.dofIndex());
            volumeVariables_.emplace_back(std::move(volVarsJ));
        }

        // maybe prepare boundary volume variables
        if (maxNumBoundaryVolVars > 0)
            CCMpfa::addBoundaryVolVars(volumeVariables_, volVarIndices_, problem, element, fvGeometry);

        // //! TODO Check if user added additional DOF dependencies, i.e. the residual of DOF globalI depends
        // //! on additional DOFs not included in the discretization schemes' occupation pattern
        // const auto& additionalDofDependencies = problem.getAdditionalDofDependencies(globalI);
        // if (!additionalDofDependencies.empty())
        // {
        //     volumeVariables_.reserve(volumeVariables_.size() + additionalDofDependencies.size());
        //     volVarIndices_.reserve(volVarIndices_.size() + additionalDofDependencies.size());
        //     for (auto globalJ : additionalDofDependencies)
        //     {
        //         const auto& elementJ = gridGeometry.element(globalJ);
        //         const auto& scvJ = fvGeometry.scv(globalJ);

        //         VolumeVariables additionalVolVars;
        //         additionalVolVars.update(elementSolution(elementJ, sol, gridGeometry),
        //                                  problem,
        //                                  elementJ,
        //                                  scvJ);

        //         volumeVariables_.emplace_back(std::move(additionalVolVars));
        //         volVarIndices_.push_back(globalJ);
        //     }
        // }
    }

    //! Prepares the volume variables of an element
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        clear();

        const auto& gridGeometry = fvGeometry.gridGeometry();
        auto eIdx = gridGeometry.elementMapper().index(element);
        volumeVariables_.resize(1);
        volVarIndices_.resize(1);

        // update the volume variables of the element
        const auto& scv = fvGeometry.scv(eIdx);
        volumeVariables_[0].update(elementSolution(element, sol, gridGeometry),
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
