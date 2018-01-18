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

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include "fluxvariablescachefiller.hh"
#include "methods.hh"

namespace Dumux
{

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The flux variables caches for an element
 * \note The class is specialized for a version with and without caching.
 *       If grid caching is enabled the flux caches are stored for the whole
 *       gridview in the corresponding GridFluxVariablesCache. This is memory
 *       intensive but faster. For caching disabled, the flux caches are locally
 *       computed for each element whenever needed.
 */
template<class TypeTag, bool EnableGridFluxVariablesCache>
class CCMpfaElementFluxVariablesCache;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The flux variables caches for an element with caching enabled
 */
template<class TypeTag>
class CCMpfaElementFluxVariablesCache<TypeTag, true>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);


public:
    //! The constructor
    CCMpfaElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    //! Specialization for the global caching being enabled - do nothing here
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) {}

    //! Specialization for the global caching being enabled - do nothing here
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) {}

    //! Specialization for the global caching being enabled - do nothing here
    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf) {}

    //! Specialization for the global caching being enabled - do nothing here
    void update(const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars) {}

    //! access operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return (*gridFluxVarsCachePtr_)[scvf]; }

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
template<class TypeTag>
class CCMpfaElementFluxVariablesCache<TypeTag, false>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GridIndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using FluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFiller<TypeTag>;
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using PrimaryIvDataHandle = typename PrimaryInteractionVolume::Traits::DataHandle;
    using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
    using SecondaryIvDataHandle = typename SecondaryInteractionVolume::Traits::DataHandle;

public:
    CCMpfaElementFluxVariablesCache(const GridFluxVariablesCache& global)
    : gridFluxVarsCachePtr_(&global) {}

    /*!
     * \brief Prepares the transmissibilities of the scv faces in an element
     * \note the fvGeometry is assumed to be bound to the same element
     * \note this function has to be called prior to flux calculations on the element.
     */
    void bindElement(const Element& element,
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
    void bind(const Element& element,
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
        // TODO should this be a property? Statically allocated memory might be an issue for local adaptivity in general
        static const std::size_t maxDiff = getParamFromGroup<std::size_t>(GET_PROP_VALUE(TypeTag, ModelParameterGroup),
                                                                          "Grid.MaxLocalElementLevelDifference", 2);
        const auto numIvEstimate = FVElementGeometry::maxNumElementScvfs*maxDiff;
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
    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        // For mpfa schemes we will have to prepare the caches of all
        // scvfs that are embedded in the interaction volumes this scvf is embedded
        DUNE_THROW(Dune::NotImplemented, "Scvf-local binding of the flux variables cache in mpfa schemes");
    }

    /*!
     * \brief Update the transmissibilities if the volume variables have changed
     * \note Results in undefined behaviour if called before bind() or with a different element
     */
    void update(const Element& element,
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
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    //! access operators in the case of no caching
    const FluxVariablesCache& operator [](const GridIndexType scvfIdx) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvfIdx)]; }

    //! access operators in the case of no caching
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    //! access operators in the case of no caching
    FluxVariablesCache& operator [](const GridIndexType scvfIdx)
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
    std::vector<GridIndexType> globalScvfIndices_;

    // stored interaction volumes and handles
    std::vector<PrimaryInteractionVolume> primaryInteractionVolumes_;
    std::vector<SecondaryInteractionVolume> secondaryInteractionVolumes_;
    std::vector<PrimaryIvDataHandle> primaryIvDataHandles_;
    std::vector<SecondaryIvDataHandle> secondaryIvDataHandles_;
};

} // end namespace

#endif
