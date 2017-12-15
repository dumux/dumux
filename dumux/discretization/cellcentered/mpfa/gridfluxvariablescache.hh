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
 * \brief The global object of flux var caches
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCMPFA_GRID_FLUXVARSCACHE_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/fluxvariablescachefiller.hh>

namespace Dumux
{

/*!
 * \ingroup Mpfa
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag, bool EnableGridFluxVariablesCache>
class CCMpfaGridFluxVariablesCache;


/*!
 * \ingroup Mpfa
 * \brief Spezialization when caching globally
 */
template<class TypeTag>
class CCMpfaGridFluxVariablesCache<TypeTag, true>
{
    // the flux variables cache filler needs to be friend to fill
    // the interaction volumes and data handles
    friend CCMpfaFluxVariablesCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
    using DataHandle = typename PrimaryInteractionVolume::Traits::DataHandle;
    using FluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFiller<TypeTag>;

public:
    CCMpfaGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    // When global caching is enabled, precompute transmissibilities for all scv faces
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // only do the update if fluxes are solution dependent or if update is forced
        if (FluxVariablesCacheFiller::isSolDependent || forceUpdate)
        {
            // clear data if forced update is desired
            if (forceUpdate)
            {
                clear_();

                const auto& gridIvIndexSets = fvGridGeometry.gridInteractionVolumeIndexSets();
                const auto numPrimaryIvs = gridIvIndexSets.numPrimaryInteractionVolumes();
                const auto numSecondaryIVs = gridIvIndexSets.numSecondaryInteractionVolumes();
                primaryInteractionVolumes_.reserve(numPrimaryIvs);
                secondaryInteractionVolumes_.reserve(numSecondaryIVs);
                primaryIvDataHandles_.reserve(numPrimaryIvs);
                secondaryIvDataHandles_.reserve(numSecondaryIVs);
            }

            // reserve memory estimate for caches, interaction volumes and corresponding data
            fluxVarsCache_.resize(fvGridGeometry.numScvf());

            // instantiate helper class to fill the caches
            FluxVariablesCacheFiller filler(problem());

            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
                // Prepare the geometries within the elements of the stencil
                auto fvGeometry = localView(fvGridGeometry);
                fvGeometry.bind(element);

                auto elemVolVars = localView(gridVolVars);
                elemVolVars.bind(element, fvGeometry, sol);

                // prepare all the caches of the scvfs inside the corresponding interaction volume
                for (const auto& scvf : scvfs(fvGeometry))
                    if (!fluxVarsCache_[scvf.index()].isUpdated())
                        filler.fill(*this, fluxVarsCache_[scvf.index()], element, fvGeometry, elemVolVars, scvf, forceUpdate);
            }
        }
    }

    void updateElement(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars)
    {
        // update only if transmissibilities are solution-dependent
        if (FluxVariablesCacheFiller::isSolDependent)
        {
            const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
            const auto& assemblyMapI = fvGridGeometry.connectivityMap()[fvGridGeometry.elementMapper().index(element)];

            // helper class to fill flux variables caches
            FluxVariablesCacheFiller filler(problem());

            // first, set all the caches to "outdated"
            for (const auto& scvf : scvfs(fvGeometry))
                fluxVarsCache_[scvf.index()].setUpdateStatus(false);
            for (const auto& dataJ : assemblyMapI)
                for (const auto scvfIdx : dataJ.scvfsJ)
                    fluxVarsCache_[scvfIdx].setUpdateStatus(false);

            // go through the caches maybe update them
            for (const auto& scvf : scvfs(fvGeometry))
            {
                auto& scvfCache = fluxVarsCache_[scvf.index()];
                if (!scvfCache.isUpdated())
                    filler.fill(*this, scvfCache, element, fvGeometry, elemVolVars, scvf);
            }

            for (const auto& dataJ : assemblyMapI)
            {
                const auto elementJ = fvGridGeometry.element(dataJ.globalJ);
                for (const auto scvfIdx : dataJ.scvfsJ)
                {
                    auto& scvfCache = fluxVarsCache_[scvfIdx];
                    if (!scvfCache.isUpdated())
                    {
                        // update cache
                        const auto& scvf = fvGeometry.scvf(scvfIdx);
                        filler.fill(*this, scvfCache, elementJ, fvGeometry, elemVolVars, scvf);
                    }
                }
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementFluxVariablesCache localView(const CCMpfaGridFluxVariablesCache& global)
    { return ElementFluxVariablesCache(global); }

    // access operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:

    void clear_()
    {
        fluxVarsCache_.clear();
        primaryInteractionVolumes_.clear();
        secondaryInteractionVolumes_.clear();
        primaryIvDataHandles_.clear();
        secondaryIvDataHandles_.clear();
    }

    const Problem* problemPtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;

    // store the interaction volumes and handles
    std::vector<PrimaryInteractionVolume> primaryInteractionVolumes_;
    std::vector<SecondaryInteractionVolume> secondaryInteractionVolumes_;
    std::vector<DataHandle> primaryIvDataHandles_;
    std::vector<DataHandle> secondaryIvDataHandles_;
};

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when not using global caching
 */
template<class TypeTag>
class CCMpfaGridFluxVariablesCache<TypeTag, false>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:
    CCMpfaGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    // When global flux variables caching is disabled, we don't need to update the cache
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false) {}

    void updateElement(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars) {}

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementFluxVariablesCache localView(const CCMpfaGridFluxVariablesCache& global)
    { return ElementFluxVariablesCache(global); }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace

#endif
