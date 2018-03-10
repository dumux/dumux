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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridVolumeVariables
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GRID_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GRID_VOLUMEVARIABLES_HH

#include <dune/common/exceptions.hh>

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Grid volume variables class for staggered models
 */
template<class FVGridGeometry, class Traits, bool enableGridVolVarsCache>
class StaggeredGridVolumeVariables;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Grid volume variables class for staggered models.
          Specialization in case of storing the volume variables
 */
template<class FVGridGeometry, class Traits>
class StaggeredGridVolumeVariables<FVGridGeometry, Traits, /*enableGridVolVarsCache*/true>
{
    using ThisType = StaggeredGridVolumeVariables<FVGridGeometry, Traits, true>;
    using Problem = typename Traits::Problem;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Indices = typename Traits::Indices;
    using IndexType = typename FVGridGeometry::GridView::IndexSet::IndexType;

public:
    //! export the type of the VolumeVariables
    using VolumeVariables = typename Traits::VolumeVariables;
    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<FVGridGeometry, ThisType, true>;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    StaggeredGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! Update all volume variables
    template<class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol)
    {
        using CellCenterPrimaryVariables = typename SolutionVector::value_type;
        auto numScv = fvGridGeometry.numScv();
        auto numBoundaryScvf = fvGridGeometry.numBoundaryScvf();

        volumeVariables_.resize(numScv + numBoundaryScvf);
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                CellCenterPrimaryVariables priVars = sol[scv.dofIndex()];
                auto elemSol = elementSolution<FVElementGeometry>(std::move(priVars));
                volumeVariables_[scv.dofIndex()].update(elemSol, problem(), element, scv);
            }

            // handle the boundary volume variables
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // if we are not on a boundary, skip the rest
                if (!scvf.boundary())
                    continue;

                const auto bcTypes = problem().boundaryTypes(element, scvf);
                const auto insideScvIdx = scvf.insideScvIdx();
                const auto& insideScv = fvGeometry.scv(insideScvIdx);

                CellCenterPrimaryVariables boundaryPriVars(0.0);

                for(int eqIdx = 0; eqIdx < CellCenterPrimaryVariables::dimension; ++eqIdx)
                {
                    if(bcTypes.isDirichlet(eqIdx) || bcTypes.isDirichletCell(eqIdx))
                        boundaryPriVars[eqIdx] = problem().dirichlet(element, scvf)[eqIdx];
                    else if(bcTypes.isNeumann(eqIdx) || bcTypes.isOutflow(eqIdx) || bcTypes.isSymmetry())
                        boundaryPriVars[eqIdx] = sol[scvf.insideScvIdx()][eqIdx];
                    //TODO: this assumes a zero-gradient for e.g. the pressure on the boundary
                    // could be made more general by allowing a non-zero-gradient, provided in problem file
                    else
                        if(eqIdx == Indices::pressureIdx)
                            DUNE_THROW(Dune::InvalidStateException, "Face at: " << scvf.center() << " has neither Dirichlet nor Neumann BC.");
                }
                auto elemSol = elementSolution<FVElementGeometry>(std::move(boundaryPriVars));
                volumeVariables_[scvf.outsideScvIdx()].update(elemSol, problem(), element, insideScv);
            }
        }
    }

    const VolumeVariables& volVars(const IndexType scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& volVars(const IndexType scvIdx)
    { return volumeVariables_[scvIdx]; }

    const VolumeVariables& volVars(const SubControlVolume& scv) const
    { return volumeVariables_[scv.dofIndex()]; }

    VolumeVariables& volVars(const SubControlVolume& scv)
    { return volumeVariables_[scv.dofIndex()]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
    std::vector<VolumeVariables> volumeVariables_;
};


/*!
 * \ingroup StaggeredDiscretization
 * \brief Grid volume variables class for staggered models.
          Specialization in case of not storing the volume variables
 */
template<class FVGridGeometry, class Traits>
class StaggeredGridVolumeVariables<FVGridGeometry, Traits, /*enableGridVolVarsCache*/false>
{
    using ThisType = StaggeredGridVolumeVariables<FVGridGeometry, Traits, false>;
    using Problem = typename Traits::Problem;

public:
    //! export the type of the VolumeVariables
    using VolumeVariables = typename Traits::VolumeVariables;
    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<FVGridGeometry, ThisType, true>;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    StaggeredGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_;}

private:

    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
