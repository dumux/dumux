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
#include <dumux/common/properties.hh>

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Grid volume variables class for staggered models
 */
template<class TypeTag, bool enableGridVolVarsCache>
class StaggeredGridVolumeVariables;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Grid volume variables class for staggered models.
          Specialization in case of storing the volume variables
 */
template<class TypeTag>
class StaggeredGridVolumeVariables<TypeTag, /*enableGridVolVarsCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using CellCenterSolutionVector = typename GET_PROP_TYPE(TypeTag, CellCenterSolutionVector);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using IndexType = typename GridView::IndexSet::IndexType;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    StaggeredGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! update all volume variables using the complete solution vector
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol)
    {
        update(fvGridGeometry, sol[cellCenterIdx]);
    }

    //! update all volume variables using the sub vector for cell center dofs
    void update(const FVGridGeometry& fvGridGeometry, const CellCenterSolutionVector& sol)
    {
        auto numScv = fvGridGeometry.numScv();
        auto numBoundaryScvf = fvGridGeometry.numBoundaryScvf();

        volumeVariables_.resize(numScv + numBoundaryScvf);
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                CellCenterPrimaryVariables priVars(0.0);
                priVars = sol[scv.dofIndex()];
                ElementSolutionVector elemSol{std::move(priVars)};
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
                ElementSolutionVector elemSol{std::move(boundaryPriVars)};
                volumeVariables_[scvf.outsideScvIdx()].update(elemSol, problem(), element, insideScv);
            }
        }
    }

    const VolumeVariables& volVars(const IndexType scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& volVars(const IndexType scvIdx)
    { return volumeVariables_[scvIdx]; }

    const VolumeVariables& volVars(const SubControlVolume scv) const
    { return volumeVariables_[scv.dofIndex()]; }

    VolumeVariables& volVars(const SubControlVolume scv)
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
template<class TypeTag>
class StaggeredGridVolumeVariables<TypeTag, /*enableGridVolVarsCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

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
