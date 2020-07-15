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
 * \ingroup CCWMpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered tpfa models
 */
#ifndef DUMUX_DISCRETIZATION_CCWMPFA_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CCWMPFA_ELEMENT_VOLUMEVARIABLES_HH

#include <algorithm>
#include <type_traits>
#include <vector>

#include <dumux/discretization/cellcentered/elementsolution.hh>

namespace Dumux {

namespace CCWMpfa {

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

        std::size_t numBoundaryVolVars = fvGeometry.numScvf();
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if(!scvf.boundary())
            {
                const auto outsideScvIdx = scvf.outsideScvIdx();
                const auto outsideElement = gridGeometry.element(outsideScvIdx);
                FVElementGeometry fvGeometryJ;
                fvGeometryJ.bindElement(outsideElement);
                if(fvGeometryJ.hasBoundaryScvf())
                    numBoundaryVolVars += fvGeometryJ.numScvf();
            }
        }

        return numBoundaryVolVars;
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
    std::pair<std::vector<VolumeVariables>, std::vector<IndexType>>
                         boundaryVolVars(const Problem& problem,
                                         const typename FVElemGeom::GridGeometry::GridView::template Codim<0>::Entity& element,
                                         const FVElemGeom& fvGeometry)
    {
        std::vector<VolumeVariables> volVars; volVars.reserve(fvGeometry.numScvf());
        std::vector<IndexType> volVarIndices; volVarIndices.reserve(fvGeometry.numScvf());
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
                else
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
        return std::make_pair(volVars, volVarIndices);
    }
} // end namespace CCMpfa

/*!
 * \ingroup CCWMpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered tpfa models
 * \note The class is specilized for versions with and without caching
 * \tparam GVV the grid volume variables type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GVV, bool cachingEnabled>
class CCWMpfaElementVolumeVariables
{};

/*!
 * \ingroup CCWMpfaDiscretization
 * \brief The local (stencil) volume variables class for cell centered tpfa models with caching
 * \note the volume variables are stored for the whole grid view in the corresponding GridVolumeVariables class
 */
template<class GVV>
class CCWMpfaElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    CCWMpfaElementVolumeVariables(const GridVolumeVariables& gridVolVars)
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
        clear_();

        if (fvGeometry.hasBoundaryScvf())
        {
            auto&& [volVars, indices] = CCWMpfa::boundaryVolVars<VolumeVariables, std::size_t>(gridVolVars().problem(), element, fvGeometry);
            std::move(volVars.begin(), volVars.end(), std::back_inserter(boundaryVolumeVariables_));
            std::move(indices.begin(), indices.end(), std::back_inserter(boundaryVolVarIndices_));
        }

        const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& assemblyMapI = gridGeometry.connectivityMap()[globalI];

        // add boundary volVars of neighbors
        for (const auto& dataJ : assemblyMapI)
        {
            const auto& elementJ = gridGeometry.element(dataJ.globalJ);
            auto fvGeometryJ = localView(gridGeometry);
            fvGeometryJ.bind(elementJ);

            if (!fvGeometryJ.hasBoundaryScvf())
                continue;

            auto&& [volVars, indices]  = CCWMpfa::boundaryVolVars<VolumeVariables, std::size_t>(gridVolVars().problem(), elementJ, fvGeometryJ);
            std::move(volVars.begin(), volVars.end(), std::back_inserter(boundaryVolumeVariables_));
            std::move(indices.begin(), indices.end(), std::back_inserter(boundaryVolVarIndices_));
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

} // end namespace Dumux

#endif
