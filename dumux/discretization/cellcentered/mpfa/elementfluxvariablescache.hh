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
 * \ingroup CCMpfaDiscretization
 * \brief The element-local object of flux variables caches
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_FLUXVARSCACHE_HH

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The flux variables caches for an element
 * \note The class is specialized for a version with and without caching.
 *       If grid caching is enabled the flux caches are stored for the whole
 *       gridview in the corresponding GridFluxVariablesCache. This is memory
 *       intensive but faster. For caching disabled, the flux caches are locally
 *       computed for each element whenever needed.
 */
template<class GFVC, bool cachingEnabled>
class CCMpfaElementFluxVariablesCache;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The flux variables caches for an element with caching enabled
 */
template<class GFVC>
class CCMpfaElementFluxVariablesCache<GFVC, true>
{
public:
    //! export the interaction volume types
    using PrimaryInteractionVolume = typename GFVC::PrimaryInteractionVolume;
    using SecondaryInteractionVolume = typename GFVC::SecondaryInteractionVolume;

    //! export the data handle types used
    using PrimaryIvDataHandle = typename GFVC::PrimaryIvDataHandle;
    using SecondaryIvDataHandle = typename GFVC::SecondaryIvDataHandle;

    //! export the flux variable cache type
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    //! export the flux variable cache filler type
    using FluxVariablesCacheFiller = typename GFVC::FluxVariablesCacheFiller;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the grid flux variables
    using GridFluxVariablesCache = GFVC;

    //! The constructor
    CCMpfaElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) {}

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) {}

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf) {}

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void update(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars) {}

    //! access operators in the case of caching
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return (*gridFluxVarsCachePtr_)[scvf]; }

    //! access to the stored interaction volumes
    const std::vector<PrimaryInteractionVolume>& primaryInteractionVolumes() const
    { return gridFluxVarsCachePtr_->primaryInteractionVolumes(); }

    //! access to the stored data handles
    const std::vector<PrimaryIvDataHandle>& primaryDataHandles() const
    { return gridFluxVarsCachePtr_->primaryDataHandles(); }

    //! access to the stored interaction volumes
    const std::vector<SecondaryInteractionVolume>& secondaryInteractionVolumes() const
    { return gridFluxVarsCachePtr_->secondaryInteractionVolumes(); }

    //! access to the stored data handles
    const std::vector<SecondaryIvDataHandle>& secondaryDataHandles() const
    { return gridFluxVarsCachePtr_->secondaryDataHandles(); }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The flux variables caches for an element with caching disabled
 */
template<class GFVC>
class CCMpfaElementFluxVariablesCache<GFVC, false>
{
public:
    //! export the interaction volume types
    using PrimaryInteractionVolume = typename GFVC::PrimaryInteractionVolume;
    using SecondaryInteractionVolume = typename GFVC::SecondaryInteractionVolume;

    //! export the data handle types used
    using PrimaryIvDataHandle = typename GFVC::PrimaryIvDataHandle;
    using SecondaryIvDataHandle = typename GFVC::SecondaryIvDataHandle;

    //! export the flux variable cache type
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;

    //! export the flux variable cache filler type
    using FluxVariablesCacheFiller = typename GFVC::FluxVariablesCacheFiller;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the grid flux variables
    using GridFluxVariablesCache = GFVC;

    CCMpfaElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    /*!
     * \brief Prepares the transmissibilities of the scv faces in an element
     * \note the fvGeometry is assumed to be bound to the same element
     * \note this function has to be called prior to flux calculations on the element.
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars)
    {
        // For mpfa schemes we will have to prepare the caches of all scvfs that are
        // embedded in the interaction volumes in which the element-local scvfs are embedded
        DUNE_THROW(Dune::NotImplemented, "Local element binding of the flux variables cache in mpfa schemes");
    }

    /*!
     * \brief Prepares the transmissibilities of the scv faces in the stencil of an element
     * \note the fvGeometry is assumed to be bound to the same element
     * \note this function has to be called prior to flux calculations on the element.
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars)
    {
        clear_();

        // some references for convenience
        const auto& problem = gridFluxVarsCache().problem();
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();

        // the assembly map of the given element
        const auto& assemblyMapI = fvGridGeometry.connectivityMap()[fvGridGeometry.elementMapper().index(element)];

        // reserve memory for scvf index container
        unsigned int numNeighborScvfs = 0;
        for (const auto& dataJ : assemblyMapI)
            numNeighborScvfs += dataJ.scvfsJ.size();
        globalScvfIndices_.resize(fvGeometry.numScvf() + numNeighborScvfs);

        // set the scvf indices in scvf index container
        unsigned int i = 0;
        for (const auto& scvf : scvfs(fvGeometry))
            globalScvfIndices_[i++] = scvf.index();
        for (const auto& dataJ : assemblyMapI)
            for (auto scvfIdx : dataJ.scvfsJ)
                globalScvfIndices_[i++] = scvfIdx;

        // Reserve memory (over-) estimate for interaction volumes and corresponding data.
        // The overestimate doesn't hurt as we are not in a memory-limited configuration.
        // We need to avoid reallocation because in the caches we store pointers to the data handles.
        // Default -> each facet has two neighbors (local adaption) and all scvfs belongs to different ivs.
        // If you want to use higher local differences change the parameter below.
        constexpr auto numIvEstimate = FVElementGeometry::maxNumElementScvfs
                                       * GridFluxVariablesCache::Traits::maxLocalElementLevelDifference();
        primaryInteractionVolumes_.reserve(numIvEstimate);
        secondaryInteractionVolumes_.reserve(numIvEstimate);
        primaryIvDataHandles_.reserve(numIvEstimate);
        secondaryIvDataHandles_.reserve(numIvEstimate);

        // helper class to fill flux variables caches
        FluxVariablesCacheFiller filler(problem);

        // resize the cache container
        fluxVarsCache_.resize(globalScvfIndices_.size());

        // go through the caches and fill them
        i = 0;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            auto& scvfCache = fluxVarsCache_[i++];
            if (!scvfCache.isUpdated())
                filler.fill(*this, scvfCache, element, fvGeometry, elemVolVars, scvf, true);
        }

        for (const auto& dataJ : assemblyMapI)
        {
            const auto elementJ = fvGridGeometry.element(dataJ.globalJ);
            for (const auto scvfIdx : dataJ.scvfsJ)
            {
                auto& scvfCache = fluxVarsCache_[i++];
                if (!scvfCache.isUpdated())
                    filler.fill(*this, scvfCache, elementJ, fvGeometry, elemVolVars, fvGeometry.scvf(scvfIdx), true);
            }
        }
    }

    /*!
     * \brief Prepares the transmissibilities of a single scv face
     * \note the fvGeometry is assumed to be bound to the same element
     * \note this function has to be called prior to flux calculations on this scvf.
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        // For mpfa schemes we will have to prepare the caches of all
        // scvfs that are embedded in the interaction volumes this scvf is embedded
        DUNE_THROW(Dune::NotImplemented, "Scvf-local binding of the flux variables cache in mpfa schemes");
    }

    /*!
     * \brief Update the transmissibilities if the volume variables have changed
     * \note Results in undefined behaviour if called before bind() or with a different element
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void update(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        // Update only if the filler puts
        // solution-dependent stuff into the caches
        if (FluxVariablesCacheFiller::isSolDependent)
        {
            const auto& problem = gridFluxVarsCache().problem();
            const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
            const auto& assemblyMapI = fvGridGeometry.connectivityMap()[fvGridGeometry.elementMapper().index(element)];

            // helper class to fill flux variables caches
            FluxVariablesCacheFiller filler(problem);

            // set all the caches to "outdated"
            for (auto& cache : fluxVarsCache_)
                cache.setUpdateStatus(false);

            // go through the caches maybe update them
            unsigned int i = 0;
            for (const auto& scvf : scvfs(fvGeometry))
            {
                auto& scvfCache = fluxVarsCache_[i++];
                if (!scvfCache.isUpdated())
                    filler.fill(*this, scvfCache, element, fvGeometry, elemVolVars, scvf);
            }

            for (const auto& dataJ : assemblyMapI)
            {
                const auto elementJ = fvGridGeometry.element(dataJ.globalJ);
                for (const auto scvfIdx : dataJ.scvfsJ)
                {
                    auto& scvfCache = fluxVarsCache_[i++];
                    if (!scvfCache.isUpdated())
                        filler.fill(*this, scvfCache, elementJ, fvGeometry, elemVolVars, fvGeometry.scvf(scvfIdx));
                }
            }
        }
    }

    //! access operators in the case of no caching
    template<class SubControlVolumeFace, typename std::enable_if_t<!std::is_integral<SubControlVolumeFace>::value, int> = 0>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    //! access operators in the case of no caching
    const FluxVariablesCache& operator [](const std::size_t scvfIdx) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvfIdx)]; }

    //! access operators in the case of no caching
    template<class SubControlVolumeFace, typename std::enable_if_t<!std::is_integral<SubControlVolumeFace>::value, int> = 0>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    //! access operators in the case of no caching
    FluxVariablesCache& operator [](const std::size_t scvfIdx)
    { return fluxVarsCache_[getLocalScvfIdx_(scvfIdx)]; }

    //! access to the stored interaction volumes
    const std::vector<PrimaryInteractionVolume>& primaryInteractionVolumes() const
    { return primaryInteractionVolumes_; }

    //! access to the stored interaction volumes
    std::vector<PrimaryInteractionVolume>& primaryInteractionVolumes()
    { return primaryInteractionVolumes_; }

    //! access to the stored data handles
    const std::vector<PrimaryIvDataHandle>& primaryDataHandles() const
    { return primaryIvDataHandles_; }

    //! access to the stored data handles
    std::vector<PrimaryIvDataHandle>& primaryDataHandles()
    { return primaryIvDataHandles_; }

    //! access to the stored interaction volumes
    const std::vector<SecondaryInteractionVolume>& secondaryInteractionVolumes() const
    { return secondaryInteractionVolumes_; }

    //! access to the stored interaction volumes
    std::vector<SecondaryInteractionVolume>& secondaryInteractionVolumes()
    { return secondaryInteractionVolumes_; }

    //! access to the stored data handles
    const std::vector<SecondaryIvDataHandle>& secondaryDataHandles() const
    { return secondaryIvDataHandles_; }

    //! access to the stored data handles
    std::vector<SecondaryIvDataHandle>& secondaryDataHandles()
    { return secondaryIvDataHandles_; }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    const GridFluxVariablesCache* gridFluxVarsCachePtr_;

    //! clears all containers
    void clear_()
    {
        fluxVarsCache_.clear();
        globalScvfIndices_.clear();
        primaryInteractionVolumes_.clear();
        secondaryInteractionVolumes_.clear();
        primaryIvDataHandles_.clear();
        secondaryIvDataHandles_.clear();
    }

    //! get index of an scvf in the local container
    unsigned int getLocalScvfIdx_(const int scvfIdx) const
    {
        auto it = std::find(globalScvfIndices_.begin(), globalScvfIndices_.end(), scvfIdx);
        assert(it != globalScvfIndices_.end() && "Could not find the flux vars cache for scvfIdx");
        return std::distance(globalScvfIndices_.begin(), it);
    }

    // the local flux vars caches and corresponding indices
    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<std::size_t> globalScvfIndices_;

    // stored interaction volumes and handles
    std::vector<PrimaryInteractionVolume> primaryInteractionVolumes_;
    std::vector<SecondaryInteractionVolume> secondaryInteractionVolumes_;
    std::vector<PrimaryIvDataHandle> primaryIvDataHandles_;
    std::vector<SecondaryIvDataHandle> secondaryIvDataHandles_;
};

} // end namespace

#endif
