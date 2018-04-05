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
 * \brief Flux variable caches on a gridview
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCMPFA_GRID_FLUXVARSCACHE_HH

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cellcentered/mpfa/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Data handle physics traits
 */
template<class ModelTraits>
struct IvDataHandlePhysicsTraits
{
    static constexpr bool enableAdvection = ModelTraits::enableAdvection();
    static constexpr bool enableMolecularDiffusion = ModelTraits::enableMolecularDiffusion();
    static constexpr bool enableHeatConduction = ModelTraits::enableEnergyBalance();

    static constexpr int numPhases = ModelTraits::numPhases();
    static constexpr int numComponents = ModelTraits::numComponents();
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Data handle physics traits
 */
template<class P,
         class FVC, class FVCF,
         class PIV, class SIV,
         class PDH, class SDH>
struct CCMpfaDefaultGridFluxVariablesCacheTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;
    using FluxVariablesCacheFiller = FVCF;

    using PrimaryInteractionVolume = PIV;
    using SecondaryInteractionVolume = SIV;
    using PrimaryIvDataHandle = PDH;
    using SecondaryIvDataHandle = SDH;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = CCMpfaElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;

    // Reserve memory (over-) estimate for interaction volumes and corresponding data.
    // The overestimate doesn't hurt as we are not in a memory-limited configuration.
    // We need to avoid reallocation because in the caches we store pointers to the data handles.
    // Default -> each facet has two neighbors (local adaption) and all scvfs belongs to different ivs.
    // If you want to use higher local differences change the parameter below.
    static constexpr std::size_t maxLocalElementLevelDifference()
    { return 2; };
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Traits, bool cachingEnabled>
class CCMpfaGridFluxVariablesCache;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class TheTraits>
class CCMpfaGridFluxVariablesCache<TheTraits, true>
{
    using Problem = typename TheTraits::Problem;
    using ThisType = CCMpfaGridFluxVariablesCache<TheTraits, true>;
public:
    //! export the Traits
    using Traits = TheTraits;

    //! export the interaction volume types
    using PrimaryInteractionVolume = typename Traits::PrimaryInteractionVolume;
    using SecondaryInteractionVolume = typename Traits::SecondaryInteractionVolume;

    //! export the data handle types used
    using PrimaryIvDataHandle = typename Traits::PrimaryIvDataHandle;
    using SecondaryIvDataHandle = typename Traits::SecondaryIvDataHandle;

    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! export the flux variable cache filler type
    using FluxVariablesCacheFiller = typename Traits::FluxVariablesCacheFiller;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! The constructor
    CCMpfaGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    //! When global caching is enabled, precompute transmissibilities for all scv faces
    template<class FVGridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // Update only if the filler puts solution-dependent
        // stuff into the caches or if update is enforced
        if (FluxVariablesCacheFiller::isSolDependent || forceUpdate)
        {
            // clear previous data if forced update is desired
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

                // reserve memory estimate for caches, interaction volumes and corresponding data
                fluxVarsCache_.resize(fvGridGeometry.numScvf());
            }

            // instantiate helper class to fill the caches
            FluxVariablesCacheFiller filler(problem());

            // set all the caches to "outdated"
            for (auto& cache : fluxVarsCache_)
                cache.setUpdateStatus(false);

            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
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

    template<class FVElementGeometry, class ElementVolumeVariables>
    void updateElement(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars)
    {
        // Update only if the filler puts
        // solution-dependent stuff into the caches
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
                        const auto& scvf = fvGeometry.scvf(scvfIdx);
                        filler.fill(*this, scvfCache, elementJ, fvGeometry, elemVolVars, scvf);
                    }
                }
            }
        }
    }

    //! access operators in the case of caching
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    //! access operators in the case of caching
    template<class SubControlVolumeFace>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

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

    const Problem& problem() const
    { return *problemPtr_; }

private:

    //! clear all containers
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
    std::vector<PrimaryIvDataHandle> primaryIvDataHandles_;
    std::vector<SecondaryIvDataHandle> secondaryIvDataHandles_;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class TheTraits>
class CCMpfaGridFluxVariablesCache<TheTraits, false>
{
    using Problem = typename TheTraits::Problem;
    using ThisType = CCMpfaGridFluxVariablesCache<TheTraits, false>;
public:
    //! export the Traits
    using Traits = TheTraits;

    //! export the interaction volume types
    using PrimaryInteractionVolume = typename Traits::PrimaryInteractionVolume;
    using SecondaryInteractionVolume = typename Traits::SecondaryInteractionVolume;

    //! export the data handle types used
    using PrimaryIvDataHandle = typename Traits::PrimaryIvDataHandle;
    using SecondaryIvDataHandle = typename Traits::SecondaryIvDataHandle;

    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! export the flux variable cache filler type
    using FluxVariablesCacheFiller = typename Traits::FluxVariablesCacheFiller;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! The constructor
    CCMpfaGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    //! When global flux variables caching is disabled, we don't need to update the cache
    template<class FVGridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false) {}

    //! When global flux variables caching is disabled, we don't need to update the cache
    template<class FVElementGeometry, class ElementVolumeVariables>
    void updateElement(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
