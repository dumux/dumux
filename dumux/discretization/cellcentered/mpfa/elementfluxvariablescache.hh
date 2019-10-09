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
 * \brief Structure to store interaction volumes and data handles
 */
template<class PrimaryIV, class PrimaryIVDataHandle,
         class SecondaryIV, class SecondaryIVDataHandle>
struct InteractionVolumeDataStorage
{
    std::vector<PrimaryIV> primaryInteractionVolumes;
    std::vector<SecondaryIV> secondaryInteractionVolumes;

    std::vector<PrimaryIVDataHandle> primaryDataHandles;
    std::vector<SecondaryIVDataHandle> secondaryDataHandles;
};

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

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the grid flux variables
    using GridFluxVariablesCache = GFVC;

private:
    //! the flux variable cache filler type
    using FluxVariablesCacheFiller = typename GFVC::Traits::FluxVariablesCacheFiller;

    //! Class to store the flux variables caches related to boundary interaction volumes
    class BoundaryCacheData
    {
        // allow the element flux variables class access to private members
        friend CCMpfaElementFluxVariablesCache<GFVC, true>;

    public:
        //! access operators
        template<class SubControlVolumeFace>
        const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
        { return fluxVarCaches_[getLocalIdx_(scvf.index())]; }

        template<class SubControlVolumeFace>
        FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
        { return fluxVarCaches_[getLocalIdx_(scvf.index())]; }

        //! clear all containers
        void clear()
        {
            fluxVarCaches_.clear();
            cacheScvfIndices_.clear();
            ivDataStorage_.primaryInteractionVolumes.clear();
            ivDataStorage_.secondaryInteractionVolumes.clear();
            ivDataStorage_.primaryDataHandles.clear();
            ivDataStorage_.secondaryDataHandles.clear();
        }

    public:
        //! map a global scvf index to the local storage index
        int getLocalIdx_(const int scvfIdx) const
        {
            auto it = std::find(cacheScvfIndices_.begin(), cacheScvfIndices_.end(), scvfIdx);
            assert(it != cacheScvfIndices_.end() && "Could not find the local idx for the given scvf idx!");
            return std::distance(cacheScvfIndices_.begin(), it);
        }

        std::vector<std::size_t> cacheScvfIndices_;
        std::vector<FluxVariablesCache> fluxVarCaches_;

        // stored boundary interaction volumes and handles
        using IVDataStorage = InteractionVolumeDataStorage<PrimaryInteractionVolume,
                                                           PrimaryIvDataHandle,
                                                           SecondaryInteractionVolume,
                                                           SecondaryIvDataHandle>;
        IVDataStorage ivDataStorage_;
    };

public:
    //! The constructor
    CCMpfaElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global)
    {}

    //! Bind the flux var caches for scvfs inside the element only
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars)
    { DUNE_THROW(Dune::NotImplemented, "Local element binding of the flux variables cache in mpfa schemes"); }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars)
    {
        boundaryCacheData_.clear();

        // find out how much memory needs to be reserved
        std::size_t numPrimaryIv;   numPrimaryIv = 0;
        std::size_t numSecondaryIv; numSecondaryIv = 0;
        std::size_t numCaches;      numCaches = 0;

        const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto& gridIvIndexSets = gridGeometry.gridInteractionVolumeIndexSets();

        // lambda to check if a scvf was handled already
        auto scvfHandled = [&] (auto idx)
        {
            return std::find(boundaryCacheData_.cacheScvfIndices_.begin(),
                             boundaryCacheData_.cacheScvfIndices_.end(),
                             idx) != boundaryCacheData_.cacheScvfIndices_.end();
        };

        // lambda to increase counters for a given scvf
        auto handleScvf = [&] (const auto& scvf, const auto& indexSet, bool isSecondary)
        {
            const auto& scvfIndices = indexSet.gridScvfIndices();
            if ( indexSet.nodalIndexSet().numBoundaryScvfs() > 0
                 && !std::any_of(scvfIndices.begin(), scvfIndices.end(), scvfHandled) )
            {
                boundaryCacheData_.cacheScvfIndices_.insert(boundaryCacheData_.cacheScvfIndices_.end(),
                                                            scvfIndices.begin(),
                                                            scvfIndices.end());
                numCaches += scvfIndices.size();
                if (isSecondary) numSecondaryIv++;
                else numPrimaryIv++;
            }
        };

        // search for ivs at boundary vertices
        for (const auto& scvf : scvfs(fvGeometry))
            gridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()) ?
                    handleScvf(scvf, gridIvIndexSets.secondaryIndexSet(scvf), true) :
                    handleScvf(scvf, gridIvIndexSets.primaryIndexSet(scvf),  false) ;

        // skip the rest if there are no boundary caches to be created
        if (numCaches > 0)
        {
            const auto& assemblyMapI = gridGeometry.connectivityMap()[gridGeometry.elementMapper().index(element)];

            for (const auto& dataJ : assemblyMapI)
            {
                for (const auto& scvfJIdx : dataJ.scvfsJ)
                {
                    const auto& scvfJ = fvGeometry.scvf(scvfJIdx);
                    if (gridGeometry.vertexUsesSecondaryInteractionVolume(scvfJ.vertexIndex()))
                        handleScvf(scvfJ, gridIvIndexSets.secondaryIndexSet(scvfJ), true);
                    else
                        handleScvf(scvfJ, gridIvIndexSets.primaryIndexSet(scvfJ), false);
                }
            }

            // now prepare the caches of all scvfs that have been found to be handled
            boundaryCacheData_.ivDataStorage_.primaryInteractionVolumes.reserve(numPrimaryIv);
            boundaryCacheData_.ivDataStorage_.secondaryInteractionVolumes.reserve(numSecondaryIv);
            boundaryCacheData_.ivDataStorage_.primaryDataHandles.reserve(numPrimaryIv);
            boundaryCacheData_.ivDataStorage_.secondaryDataHandles.reserve(numSecondaryIv);

            boundaryCacheData_.fluxVarCaches_.resize(numCaches);
            for (auto& cache : boundaryCacheData_.fluxVarCaches_)
                cache.setUpdateStatus(false);

            FluxVariablesCacheFiller filler(gridFluxVarsCachePtr_->problem());
            for (auto scvfIdx : boundaryCacheData_.cacheScvfIndices_)
            {
                const auto& scvf = fvGeometry.scvf(scvfIdx);
                auto& cache = boundaryCacheData_[scvf];
                if (!cache.isUpdated())
                    filler.fill(boundaryCacheData_, cache, boundaryCacheData_.ivDataStorage_,
                                element, fvGeometry, elemVolVars, scvf, /*forceUpdate*/true);
            }
        }
    }

    //! Bind the flux var caches for an individual scvf
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf)
    { DUNE_THROW(Dune::NotImplemented, "Scvf-local binding of the flux variables cache in mpfa schemes"); }

    //! Specialization for the global caching being enabled - do nothing here
    template<class FVElementGeometry, class ElementVolumeVariables>
    void update(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        // Update only if the filler puts solution-dependent stuff into the caches
        if (FluxVariablesCacheFiller::isSolDependent)
        {
            // helper class to fill flux variables caches
            FluxVariablesCacheFiller filler(gridFluxVarsCachePtr_->problem());

            // first, set all the caches to "outdated"
            for (auto& cache : boundaryCacheData_.fluxVarCaches_)
                cache.setUpdateStatus(false);

            // go through the caches maybe update them
            std::size_t cacheIdx = 0;
            for (auto scvfIdx : boundaryCacheData_.cacheScvfIndices_)
            {
                auto& scvfCache = boundaryCacheData_.fluxVarCaches_[cacheIdx++];
                if (!scvfCache.isUpdated())
                    filler.fill(boundaryCacheData_, scvfCache, boundaryCacheData_.ivDataStorage_,
                                element, fvGeometry, elemVolVars, fvGeometry.scvf(scvfIdx));
            }
        }
    }

    //! access operators in the case of caching
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return !isEmbeddedInBoundaryIV_(scvf) ? (*gridFluxVarsCachePtr_)[scvf] : boundaryCacheData_[scvf]; }

    //! access to the interaction volume an scvf is embedded in
    template<class SubControlVolumeFace>
    const PrimaryInteractionVolume& primaryInteractionVolume(const SubControlVolumeFace& scvf) const
    {
        return isEmbeddedInBoundaryIV_(scvf)
               ? boundaryCacheData_.ivDataStorage_.primaryInteractionVolumes[ (*this)[scvf].ivIndexInContainer() ]
               : gridFluxVarsCachePtr_->primaryInteractionVolume(scvf);
    }

    //! access to the data handle of an interaction volume an scvf is embedded in
    template<class SubControlVolumeFace>
    const PrimaryIvDataHandle& primaryDataHandle(const SubControlVolumeFace& scvf) const
    {
        return isEmbeddedInBoundaryIV_(scvf)
               ? boundaryCacheData_.ivDataStorage_.primaryDataHandles[ (*this)[scvf].ivIndexInContainer() ]
               : gridFluxVarsCachePtr_->primaryDataHandle(scvf);
    }

    //! access to the interaction volume an scvf is embedded in
    template<class SubControlVolumeFace>
    const SecondaryInteractionVolume& secondaryInteractionVolume(const SubControlVolumeFace& scvf) const
    {
        return isEmbeddedInBoundaryIV_(scvf)
               ? boundaryCacheData_.ivDataStorage_.secondaryInteractionVolumes[ (*this)[scvf].ivIndexInContainer() ]
               : gridFluxVarsCachePtr_->secondaryInteractionVolume(scvf);
    }

    //! access to the data handle of an interaction volume an scvf is embedded in
    template<class SubControlVolumeFace>
    const SecondaryIvDataHandle& secondaryDataHandle(const SubControlVolumeFace& scvf) const
    {
        return isEmbeddedInBoundaryIV_(scvf)
               ? boundaryCacheData_.ivDataStorage_.secondaryDataHandles[ (*this)[scvf].ivIndexInContainer() ]
               : gridFluxVarsCachePtr_->secondaryDataHandle(scvf);
    }

    //! The global object we are a restriction of
    const GridFluxVariablesCache& gridFluxVarsCache() const
    {  return *gridFluxVarsCachePtr_; }

private:
    //! returns true if an scvf is contained in an interaction volume that touches the boundary
    template<class SubControlVolumeFace>
    bool isEmbeddedInBoundaryIV_(const SubControlVolumeFace& scvf) const
    {
        const auto& gridGeometry = gridFluxVarsCachePtr_->problem().gridGeometry();
        const auto& gridIvIndexSets = gridGeometry.gridInteractionVolumeIndexSets();
        if (gridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
            return gridIvIndexSets.secondaryIndexSet(scvf).nodalIndexSet().numBoundaryScvfs() > 0;
        else
            return gridIvIndexSets.primaryIndexSet(scvf).nodalIndexSet().numBoundaryScvfs() > 0;
    }

    const GridFluxVariablesCache* gridFluxVarsCachePtr_;

    // we store those caches that touch the boundary locally here
    // for the case that the boundary conditions change, which would
    // leave the grid-wide cache outdated.
    BoundaryCacheData boundaryCacheData_;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The flux variables caches for an element with caching disabled
 */
template<class GFVC>
class CCMpfaElementFluxVariablesCache<GFVC, false>
{
    //! the flux variable cache filler type
    using FluxVariablesCacheFiller = typename GFVC::Traits::FluxVariablesCacheFiller;

public:
    //! export the interaction volume types
    using PrimaryInteractionVolume = typename GFVC::PrimaryInteractionVolume;
    using SecondaryInteractionVolume = typename GFVC::SecondaryInteractionVolume;

    //! export the data handle types used
    using PrimaryIvDataHandle = typename GFVC::PrimaryIvDataHandle;
    using SecondaryIvDataHandle = typename GFVC::SecondaryIvDataHandle;

    //! export the flux variable cache type
    using FluxVariablesCache = typename GFVC::FluxVariablesCache;


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
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
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
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars)
    {
        clear_();

        // some references for convenience
        const auto& problem = gridFluxVarsCache().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();

        // the assembly map of the given element
        const auto& assemblyMapI = gridGeometry.connectivityMap()[gridGeometry.elementMapper().index(element)];

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
        ivDataStorage_.primaryInteractionVolumes.reserve(numIvEstimate);
        ivDataStorage_.secondaryInteractionVolumes.reserve(numIvEstimate);
        ivDataStorage_.primaryDataHandles.reserve(numIvEstimate);
        ivDataStorage_.secondaryDataHandles.reserve(numIvEstimate);

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
                filler.fill(*this, scvfCache, ivDataStorage_, element, fvGeometry, elemVolVars, scvf, true);
        }

        for (const auto& dataJ : assemblyMapI)
        {
            const auto elementJ = gridGeometry.element(dataJ.globalJ);
            for (const auto scvfIdx : dataJ.scvfsJ)
            {
                auto& scvfCache = fluxVarsCache_[i++];
                if (!scvfCache.isUpdated())
                    filler.fill(*this, scvfCache, ivDataStorage_, elementJ, fvGeometry, elemVolVars, fvGeometry.scvf(scvfIdx), true);
            }
        }
    }

    /*!
     * \brief Prepares the transmissibilities of a single scv face
     * \note the fvGeometry is assumed to be bound to the same element
     * \note this function has to be called prior to flux calculations on this scvf.
     */
    template<class FVElementGeometry, class ElementVolumeVariables>
    void bindScvf(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
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
    void update(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        // Update only if the filler puts
        // solution-dependent stuff into the caches
        if (FluxVariablesCacheFiller::isSolDependent)
        {
            const auto& problem = gridFluxVarsCache().problem();
            const auto& gridGeometry = fvGeometry.gridGeometry();
            const auto& assemblyMapI = gridGeometry.connectivityMap()[gridGeometry.elementMapper().index(element)];

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
                    filler.fill(*this, scvfCache, ivDataStorage_, element, fvGeometry, elemVolVars, scvf);
            }

            for (const auto& dataJ : assemblyMapI)
            {
                const auto elementJ = gridGeometry.element(dataJ.globalJ);
                for (const auto scvfIdx : dataJ.scvfsJ)
                {
                    auto& scvfCache = fluxVarsCache_[i++];
                    if (!scvfCache.isUpdated())
                        filler.fill(*this, scvfCache, ivDataStorage_, elementJ, fvGeometry, elemVolVars, fvGeometry.scvf(scvfIdx));
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

    //! access to the interaction volume an scvf is embedded in
    template<class SubControlVolumeFace>
    const PrimaryInteractionVolume& primaryInteractionVolume(const SubControlVolumeFace& scvf) const
    { return ivDataStorage_.primaryInteractionVolumes[ (*this)[scvf].ivIndexInContainer() ]; }

    //! access to the data handle of an interaction volume an scvf is embedded in
    template<class SubControlVolumeFace>
    const PrimaryIvDataHandle& primaryDataHandle(const SubControlVolumeFace& scvf) const
    { return ivDataStorage_.primaryDataHandles[ (*this)[scvf].ivIndexInContainer() ]; }

    //! access to the interaction volume an scvf is embedded in
    template<class SubControlVolumeFace>
    const SecondaryInteractionVolume& secondaryInteractionVolume(const SubControlVolumeFace& scvf) const
    { return ivDataStorage_.secondaryInteractionVolumes[ (*this)[scvf].ivIndexInContainer() ]; }

    //! access to the data handle of an interaction volume an scvf is embedded in
    template<class SubControlVolumeFace>
    const SecondaryIvDataHandle& secondaryDataHandle(const SubControlVolumeFace& scvf) const
    { return ivDataStorage_.secondaryDataHandles[ (*this)[scvf].ivIndexInContainer() ]; }

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
        ivDataStorage_.primaryInteractionVolumes.clear();
        ivDataStorage_.secondaryInteractionVolumes.clear();
        ivDataStorage_.primaryDataHandles.clear();
        ivDataStorage_.secondaryDataHandles.clear();
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
    using IVDataStorage = InteractionVolumeDataStorage<PrimaryInteractionVolume,
                                                       PrimaryIvDataHandle,
                                                       SecondaryInteractionVolume,
                                                       SecondaryIvDataHandle>;
    IVDataStorage ivDataStorage_;
};

} // end namespace

#endif
