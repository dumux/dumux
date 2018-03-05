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
 * \copydoc Dumux::StaggeredGridFluxVariablesCache
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GRID_FLUXVARSCACHE_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/staggered/elementfluxvariablescache.hh>

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models
 */
template<class TypeTag, bool EnableGridFluxVariablesCache>
class StaggeredGridFluxVariablesCache;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models.
          Specialization in case of storing the flux cache.
 */
template<class TypeTag>
class StaggeredGridFluxVariablesCache<TypeTag, true>
{
    // the local class needs access to the problem
    friend StaggeredElementFluxVariablesCache<TypeTag, true>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    StaggeredGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // fluxVarsCache_.resize(fvGridGeometry.numScvf());
        // for (const auto& element : elements(fvGridGeometry.gridView()))
        // {
        //     // Prepare the geometries within the elements of the stencil
        //     auto fvGeometry = localView(fvGridGeometry);
        //     fvGeometry.bind(element);
        //
        //     auto elemVolVars = localView(gridVolVars);
        //     elemVolVars.bind(element, fvGeometry, sol);
        //
        //     for (auto&& scvf : scvfs(fvGeometry))
        //     {
        //         fluxVarsCache_[scvf.index()].update(problem, element, fvGeometry, elemVolVars, scvf);
        //     }
        // }
    }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    // access operators in the case of caching
    const FluxVariablesCache& operator [](IndexType scvfIdx) const
    { return fluxVarsCache_[scvfIdx]; }

    FluxVariablesCache& operator [](IndexType scvfIdx)
    { return fluxVarsCache_[scvfIdx]; }

    const Problem* problemPtr_;

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models.
          Specialization in case of not storing the flux cache.
 */
template<class TypeTag>
class StaggeredGridFluxVariablesCache<TypeTag, false>
{
    // the local class needs access to the problem
    friend StaggeredElementFluxVariablesCache<TypeTag, false>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    // When global flux variables caching is disabled, we don't need to update the cache
    void update(Problem& problem)
    { problemPtr_ = &problem; }

private:

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
