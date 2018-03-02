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
 * \ingroup CCDiscretization
 * \brief The grid volume variables class for cell centered models
 */
#ifndef DUMUX_DISCRETIZATION_CC_GRID_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CC_GRID_VOLUMEVARIABLES_HH

#include <dumux/common/properties.hh>

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Base class for the grid volume variables
 * \note This class has a cached version and a non-cached version
 * \tparam TypeTag the TypeTag
 * \tparam enableGridVolVarsCache if the cache is enabled
 */
template<class TypeTag, bool enableGridVolVarsCache>
class CCGridVolumeVariables
{};

//! specialization in case of storing the volume variables
template<class TypeTag>
class CCGridVolumeVariables<TypeTag, /*enableGridVolVarsCache*/true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    CCGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol)
    {
        const auto numScv = fvGridGeometry.numScv();
        const auto numBoundaryScvf = fvGridGeometry.numBoundaryScvf();

        volumeVariables_.resize(numScv + numBoundaryScvf);
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                const ElementSolution elemSol({sol[scv.dofIndex()]});
                volumeVariables_[scv.dofIndex()].update(elemSol, problem(), element, scv);
            }

            // handle the boundary volume variables
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // if we are not on a boundary, skip the rest
                if (!scvf.boundary())
                    continue;

                // check if boundary is a pure dirichlet boundary
                const auto bcTypes = problem().boundaryTypes(element, scvf);
                if (bcTypes.hasOnlyDirichlet())
                {
                    const auto insideScvIdx = scvf.insideScvIdx();
                    const auto& insideScv = fvGeometry.scv(insideScvIdx);
                    const ElementSolution dirichletPriVars({problem().dirichlet(element, scvf)});

                    volumeVariables_[scvf.outsideScvIdx()].update(dirichletPriVars,
                                                                  problem(),
                                                                  element,
                                                                  insideScv);
                }
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

    // required for compatibility with the box method
    const VolumeVariables& volVars(const IndexType scvIdx, const IndexType localIdx) const
    { return volumeVariables_[scvIdx]; }

    // required for compatibility with the box method
    VolumeVariables& volVars(const IndexType scvIdx, const IndexType localIdx)
    { return volumeVariables_[scvIdx]; }

    //! The problem we are solving
    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
    std::vector<VolumeVariables> volumeVariables_;
};


//! Specialization when the current volume variables are not stored globally
template<class TypeTag>
class CCGridVolumeVariables<TypeTag, /*enableGridVolVarsCache*/false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    CCGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol) {}

    //! The problem we are solving
    const Problem& problem() const
    { return *problemPtr_;}

private:
    const Problem* problemPtr_;
};

} // end namespace

#endif
