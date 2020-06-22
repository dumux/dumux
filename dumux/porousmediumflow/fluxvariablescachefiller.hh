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
 * \ingroup PorousmediumflowModels
 * \brief A helper class to fill the flux variables cache
 */
#ifndef DUMUX_POROUSMEDIUM_FLUXVARIABLESCACHE_FILLER_HH
#define DUMUX_POROUSMEDIUM_FLUXVARIABLESCACHE_FILLER_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class PorousMediumFluxVariablesCacheFillerImplementation;

/*!
 * \ingroup PorousmediumflowModels
 * \brief The flux variables cache filler class for porous media
 *
 * Helps filling the flux variables cache depending several policies
 */
template<class TypeTag>
using PorousMediumFluxVariablesCacheFiller = PorousMediumFluxVariablesCacheFillerImplementation<TypeTag, GetPropType<TypeTag, Properties::GridGeometry>::discMethod>;

//! Specialization of the flux variables cache filler for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheFillerImplementation<TypeTag, DiscretizationMethod::cctpfa>
{
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool advectionEnabled = ModelTraits::enableAdvection();
    static constexpr bool diffusionEnabled = ModelTraits::enableMolecularDiffusion();
    static constexpr bool heatConductionEnabled = ModelTraits::enableEnergyBalance();

    static constexpr bool advectionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentAdvection>();
    static constexpr bool diffusionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentMolecularDiffusion>();
    static constexpr bool heatConductionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentHeatConduction>();


public:
    static constexpr bool isSolDependent = (advectionEnabled && advectionIsSolDependent) ||
                                           (diffusionEnabled && diffusionIsSolDependent) ||
                                           (heatConductionEnabled && heatConductionIsSolDependent);

    //! The constructor. Sets the problem pointer
    PorousMediumFluxVariablesCacheFillerImplementation(const Problem& problem)
    : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheContainer Either the element or global flux variables cache
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables
     * \param scvf The corresponding sub-control volume face
     * \param forceUpdateAll if true, forces all caches to be updated (even the solution-independent ones)
     */
    template<class FluxVariablesCacheContainer, class FluxVariablesCache>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              bool forceUpdateAll = false)
    {
        // fill the physics-related quantities of the caches
        if (forceUpdateAll)
        {
            if constexpr (advectionEnabled)
                fillAdvection_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
            if constexpr (diffusionEnabled)
                fillDiffusion_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
            if constexpr (heatConductionEnabled)
                fillHeatConduction_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
        }
        else
        {
            if constexpr (advectionEnabled && advectionIsSolDependent)
                fillAdvection_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
            if constexpr (diffusionEnabled && diffusionIsSolDependent)
                fillDiffusion_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
            if constexpr (heatConductionEnabled && heatConductionIsSolDependent)
                fillHeatConduction_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
        }
    }

private:

    const Problem& problem() const
    { return *problemPtr_; }

    //! method to fill the advective quantities
    template<class FluxVariablesCache>
    void fillAdvection_(FluxVariablesCache& scvfFluxVarsCache,
                        const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf)
    {
        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
        using AdvectionFiller = typename AdvectionType::Cache::Filler;

        // forward to the filler for the advective quantities
        AdvectionFiller::fill(scvfFluxVarsCache, problem(), element, fvGeometry, elemVolVars, scvf, *this);
    }

    //! method to fill the diffusive quantities
    template<class FluxVariablesCache>
    void fillDiffusion_(FluxVariablesCache& scvfFluxVarsCache,
                        const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf)
    {
        using DiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
        using DiffusionFiller = typename DiffusionType::Cache::Filler;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        static constexpr int numPhases = ModelTraits::numFluidPhases();
        static constexpr int numComponents = ModelTraits::numFluidComponents();

        // forward to the filler of the diffusive quantities
        if constexpr (FluidSystem::isTracerFluidSystem())
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
                    DiffusionFiller::fill(scvfFluxVarsCache, phaseIdx, compIdx, problem(), element, fvGeometry, elemVolVars, scvf, *this);
        else
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
                    if (compIdx != FluidSystem::getMainComponent(phaseIdx))
                        DiffusionFiller::fill(scvfFluxVarsCache, phaseIdx, compIdx, problem(), element, fvGeometry, elemVolVars, scvf, *this);
    }

    //! method to fill the quantities related to heat conduction
    template<class FluxVariablesCache>
    void fillHeatConduction_(FluxVariablesCache& scvfFluxVarsCache,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf)
    {
        using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;
        using HeatConductionFiller = typename HeatConductionType::Cache::Filler;

        // forward to the filler of the diffusive quantities
        HeatConductionFiller::fill(scvfFluxVarsCache, problem(), element, fvGeometry, elemVolVars, scvf, *this);
    }

    const Problem* problemPtr_;
};

//! Specialization of the flux variables cache filler for the cell centered mpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheFillerImplementation<TypeTag, DiscretizationMethod::ccmpfa>
{
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using MpfaHelper = typename GridGeometry::MpfaHelper;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;

    using PrimaryInteractionVolume = GetPropType<TypeTag, Properties::PrimaryInteractionVolume>;
    using PrimaryDataHandle = typename ElementFluxVariablesCache::PrimaryIvDataHandle;
    using PrimaryLocalFaceData = typename PrimaryInteractionVolume::Traits::LocalFaceData;
    using SecondaryInteractionVolume = GetPropType<TypeTag, Properties::SecondaryInteractionVolume>;
    using SecondaryDataHandle = typename ElementFluxVariablesCache::SecondaryIvDataHandle;
    using SecondaryLocalFaceData = typename SecondaryInteractionVolume::Traits::LocalFaceData;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static constexpr bool advectionEnabled = ModelTraits::enableAdvection();
    static constexpr bool diffusionEnabled = ModelTraits::enableMolecularDiffusion();
    static constexpr bool heatConductionEnabled = ModelTraits::enableEnergyBalance();

    static constexpr bool advectionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentAdvection>();
    static constexpr bool diffusionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentMolecularDiffusion>();
    static constexpr bool heatConductionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentHeatConduction>();

public:
    //! This cache filler is always solution-dependent, as it updates the
    //! vectors of cell unknowns with which the transmissibilities have to be
    //! multiplied in order to obtain the fluxes.
    static constexpr bool isSolDependent = true;

    //! The constructor. Sets problem pointer.
    PorousMediumFluxVariablesCacheFillerImplementation(const Problem& problem)
    : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheStorage Class that holds the scvf flux vars caches
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param ivDataStorage Class that stores the interaction volumes & handles
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables (primary/secondary variables)
     * \param scvf The corresponding sub-control volume face
     * \param forceUpdateAll if true, forces all caches to be updated (even the solution-independent ones)
     */
    template<class FluxVarsCacheStorage, class FluxVariablesCache, class IVDataStorage>
    void fill(FluxVarsCacheStorage& fluxVarsCacheStorage,
              FluxVariablesCache& scvfFluxVarsCache,
              IVDataStorage& ivDataStorage,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              bool forceUpdateAll = false)
    {
        // Set pointers
        elementPtr_ = &element;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
        const auto& gridGeometry = fvGeometry.gridGeometry();

        // 1. prepare interaction volume (iv)
        // 2. solve for all transmissibilities and store them in data handles
        // 3. set pointers to transmissibilities in caches of all the scvfs of the iv
        if (gridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
        {
            if (forceUpdateAll)
            {
                // create new interaction volume
                const auto ivIndexInContainer = ivDataStorage.secondaryInteractionVolumes.size();
                const auto& indexSet = gridGeometry.gridInteractionVolumeIndexSets().secondaryIndexSet(scvf);
                ivDataStorage.secondaryInteractionVolumes.emplace_back();
                secondaryIv_ = &ivDataStorage.secondaryInteractionVolumes.back();
                secondaryIv_->bind(indexSet, problem(), fvGeometry);

                // create the corresponding data handle
                ivDataStorage.secondaryDataHandles.emplace_back();
                secondaryIvDataHandle_ = &ivDataStorage.secondaryDataHandles.back();
                prepareDataHandle_(*secondaryIv_, *secondaryIvDataHandle_, forceUpdateAll);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_<FluxVariablesCache>(fluxVarsCacheStorage, *secondaryIv_, ivIndexInContainer);
            }
            else
            {
                // get previously created interaction volume/handle
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                secondaryIv_ = &ivDataStorage.secondaryInteractionVolumes[ivIndexInContainer];
                secondaryIvDataHandle_ = &ivDataStorage.secondaryDataHandles[ivIndexInContainer];
                prepareDataHandle_(*secondaryIv_, *secondaryIvDataHandle_, forceUpdateAll);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_<FluxVariablesCache>(fluxVarsCacheStorage, *secondaryIv_, ivIndexInContainer);
            }
        }

        // primary interaction volume type
        else
        {
            if (forceUpdateAll)
            {
                // create new interaction volume
                const auto ivIndexInContainer = ivDataStorage.primaryInteractionVolumes.size();
                const auto& indexSet = gridGeometry.gridInteractionVolumeIndexSets().primaryIndexSet(scvf);
                ivDataStorage.primaryInteractionVolumes.emplace_back();
                primaryIv_ = &ivDataStorage.primaryInteractionVolumes.back();
                primaryIv_->bind(indexSet, problem(), fvGeometry);

                // create the corresponding data handle
                ivDataStorage.primaryDataHandles.emplace_back();
                primaryIvDataHandle_ = &ivDataStorage.primaryDataHandles.back();
                prepareDataHandle_(*primaryIv_, *primaryIvDataHandle_, forceUpdateAll);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_<FluxVariablesCache>(fluxVarsCacheStorage, *primaryIv_, ivIndexInContainer);
            }
            else
            {
                // get previously created interaction volume/handle
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                primaryIv_ = &ivDataStorage.primaryInteractionVolumes[ivIndexInContainer];
                primaryIvDataHandle_ = &ivDataStorage.primaryDataHandles[ivIndexInContainer];
                prepareDataHandle_(*primaryIv_, *primaryIvDataHandle_, forceUpdateAll);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_<FluxVariablesCache>(fluxVarsCacheStorage, *primaryIv_, ivIndexInContainer);
            }
        }
    }

    //! returns the stored interaction volume pointer
    const PrimaryInteractionVolume& primaryInteractionVolume() const
    { return *primaryIv_; }

    //! returns the stored interaction volume pointer
    const SecondaryInteractionVolume& secondaryInteractionVolume() const
    { return *secondaryIv_; }

    //! returns the stored data handle pointer
    const PrimaryDataHandle& primaryIvDataHandle() const
    { return *primaryIvDataHandle_; }

    //! returns the stored data handle pointer
    const SecondaryDataHandle& secondaryIvDataHandle() const
    { return *secondaryIvDataHandle_; }

    //! returns the currently stored iv-local face data object
    const PrimaryLocalFaceData& primaryIvLocalFaceData() const
    { return *primaryLocalFaceData_; }

    //! returns the currently stored iv-local face data object
    const SecondaryLocalFaceData& secondaryIvLocalFaceData() const
    { return *secondaryLocalFaceData_; }

private:

    const Problem& problem() const { return *problemPtr_; }
    const Element& element() const { return *elementPtr_; }
    const FVElementGeometry& fvGeometry() const { return *fvGeometryPtr_; }
    const ElementVolumeVariables& elemVolVars() const { return *elemVolVarsPtr_; }

    //! Method to fill the flux var caches within an interaction volume
    template<class FluxVariablesCache, class FluxVarsCacheStorage, class InteractionVolume>
    void fillCachesInInteractionVolume_(FluxVarsCacheStorage& fluxVarsCacheStorage,
                                        InteractionVolume& iv,
                                        unsigned int ivIndexInContainer)
    {
        // determine if secondary interaction volumes are used here
        static constexpr bool isSecondary = MpfaHelper::considerSecondaryIVs()
                                            && std::is_same_v<InteractionVolume, SecondaryInteractionVolume>;

        // First we upate data which are not dependent on the physical processes.
        // We store pointers to the other flux var caches, so that we have to obtain
        // this data only once and can use it again in the sub-cache fillers.
        const auto numGlobalScvfs = iv.localFaceData().size();
        std::vector<const SubControlVolumeFace*> ivScvfs(numGlobalScvfs);
        std::vector<FluxVariablesCache*> ivFluxVarCaches(numGlobalScvfs);

        unsigned int i = 0;
        for (const auto& d : iv.localFaceData())
        {
            // obtain the scvf
            const auto& scvfJ = fvGeometry().scvf(d.gridScvfIndex());
            ivScvfs[i] = &scvfJ;
            ivFluxVarCaches[i] = &fluxVarsCacheStorage[scvfJ];
            ivFluxVarCaches[i]->setIvIndexInContainer(ivIndexInContainer);
            ivFluxVarCaches[i]->setUpdateStatus(true);
            ivFluxVarCaches[i]->setSecondaryIvUsage(isSecondary);
            ivFluxVarCaches[i]->setIvLocalFaceIndex(d.ivLocalScvfIndex());
            if (dim < dimWorld)
                if (d.isOutsideFace())
                    ivFluxVarCaches[i]->setIndexInOutsideFaces(d.scvfLocalOutsideScvfIndex());
            i++;
        }

        if constexpr (advectionEnabled)
            fillAdvection_(iv, ivScvfs, ivFluxVarCaches);
        if constexpr (diffusionEnabled)
            fillDiffusion_(iv, ivScvfs, ivFluxVarCaches);
        if constexpr (heatConductionEnabled)
            fillHeatConduction_(iv, ivScvfs, ivFluxVarCaches);
    }

    //! fills the advective quantities (enabled advection)
    template<class InteractionVolume, class FluxVariablesCache>
    void fillAdvection_(InteractionVolume& iv,
                        const std::vector<const SubControlVolumeFace*>& ivScvfs,
                        const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
        using AdvectionFiller = typename AdvectionType::Cache::Filler;

        // fill advection caches
        for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
        {
            // set pointer to current local face data object
            // ifs are evaluated at compile time and are optimized away
            if (std::is_same_v<PrimaryInteractionVolume, SecondaryInteractionVolume>)
            {
                // we cannot make a disctinction, thus we set both pointers
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
            }
            else if (std::is_same_v<InteractionVolume, PrimaryInteractionVolume>)
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
            else
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);

            // fill this scvfs cache
            AdvectionFiller::fill(*ivFluxVarCaches[i],
                                  problem(),
                                  iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                  fvGeometry(),
                                  elemVolVars(),
                                  *ivScvfs[i],
                                  *this);
        }
    }

    //! fills the diffusive quantities (diffusion enabled)
    template<class InteractionVolume, class FluxVariablesCache>
    void fillDiffusion_(InteractionVolume& iv,
                        const std::vector<const SubControlVolumeFace*>& ivScvfs,
                        const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using DiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
        using DiffusionFiller = typename DiffusionType::Cache::Filler;

        static constexpr int numPhases = ModelTraits::numFluidPhases();
        static constexpr int numComponents = ModelTraits::numFluidComponents();

        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
                if constexpr (!FluidSystem::isTracerFluidSystem())
                    if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                        continue;

                // fill diffusion caches
                for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
                {
                    // set pointer to current local face data object
                    // ifs are evaluated at compile time and are optimized away
                    if constexpr (std::is_same_v<PrimaryInteractionVolume, SecondaryInteractionVolume>)
                    {
                        // we cannot make a disctinction, thus we set both pointers
                        primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                        secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
                    }
                    else if constexpr (std::is_same_v<InteractionVolume, PrimaryInteractionVolume>)
                        primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                    else
                        secondaryLocalFaceData_ = &(iv.localFaceData()[i]);

                    // fill this scvfs cache
                    DiffusionFiller::fill(*ivFluxVarCaches[i],
                                          phaseIdx,
                                          compIdx,
                                          problem(),
                                          iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                          fvGeometry(),
                                          elemVolVars(),
                                          *ivScvfs[i],
                                          *this);
                }
            }
        }
    }

    //! fills the quantities related to heat conduction (heat conduction enabled)
    template<class InteractionVolume, class FluxVariablesCache>
    void fillHeatConduction_(InteractionVolume& iv,
                             const std::vector<const SubControlVolumeFace*>& ivScvfs,
                             const std::vector<FluxVariablesCache*>& ivFluxVarCaches)
    {
        using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;
        using HeatConductionFiller = typename HeatConductionType::Cache::Filler;

        // fill heat conduction caches
        for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
        {
            // set pointer to current local face data object
            // ifs are evaluated at compile time and are optimized away
            if constexpr (std::is_same_v<PrimaryInteractionVolume, SecondaryInteractionVolume>)
            {
                // we cannot make a disctinction, thus we set both pointers
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
            }
            else if constexpr (std::is_same_v<InteractionVolume, PrimaryInteractionVolume>)
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
            else
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);

            // fill this scvfs cache
            HeatConductionFiller::fill(*ivFluxVarCaches[i],
                                       problem(),
                                       iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                       fvGeometry(),
                                       elemVolVars(),
                                       *ivScvfs[i],
                                       *this);
        }
    }

    //! Solves the local systems and stores the result in the handles
    template< class InteractionVolume, class DataHandle>
    void prepareDataHandle_([[maybe_unused]] InteractionVolume& iv, [[maybe_unused]] DataHandle& handle, [[maybe_unused]] bool forceUpdate)
    {
        // (maybe) solve system subject to intrinsic permeability
        if constexpr (advectionEnabled)
        {
            using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
            if constexpr (AdvectionType::discMethod == DiscretizationMethod::ccmpfa)
                prepareAdvectionHandle_(iv, handle, forceUpdate);
        }

        // (maybe) solve system subject to diffusion tensors
        if constexpr (diffusionEnabled)
        {
            using DiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
            if constexpr (DiffusionType::discMethod == DiscretizationMethod::ccmpfa)
                prepareDiffusionHandles_(iv, handle, forceUpdate);
        }

        // (maybe) solve system subject to thermal conductivity
        if constexpr (heatConductionEnabled)
        {
            using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;
            if constexpr (HeatConductionType::discMethod == DiscretizationMethod::ccmpfa)
                prepareHeatConductionHandle_(iv, handle, forceUpdate);
        }
    }

    //! prepares the quantities necessary for advective fluxes in the handle
    template<class InteractionVolume, class DataHandle>
    void prepareAdvectionHandle_(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll)
    {
        // get instance of the interaction volume-local assembler
        using Traits = typename InteractionVolume::Traits;
        using IvLocalAssembler = typename Traits::template LocalAssembler<Problem, FVElementGeometry, ElementVolumeVariables>;
        IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

        // lambda to obtain the permeability tensor
        auto getK = [] (const auto& volVars) { return volVars.permeability(); };

        // Assemble T only if permeability is sol-dependent or if update is forced
        if (forceUpdateAll || advectionIsSolDependent)
            localAssembler.assembleMatrices(handle.advectionHandle(), iv, getK);

        // assemble pressure vectors
        for (unsigned int pIdx = 0; pIdx < ModelTraits::numFluidPhases(); ++pIdx)
        {
            // set context in handle
            handle.advectionHandle().setPhaseIndex(pIdx);

            // maybe (re-)assemble gravity contribution vector
            auto getRho = [pIdx] (const auto& volVars) { return volVars.density(pIdx); };
            static const bool enableGravity = getParamFromGroup<bool>(problem().paramGroup(), "Problem.EnableGravity");
            if (enableGravity)
                localAssembler.assembleGravity(handle.advectionHandle(), iv, getRho);

            // reassemble pressure vector
            auto getPressure = [pIdx] (const auto& volVars) { return volVars.pressure(pIdx); };
            localAssembler.assembleU(handle.advectionHandle(), iv, getPressure);
        }
    }

    //! prepares the quantities necessary for diffusive fluxes in the handle
    template<class InteractionVolume, class DataHandle>
    void prepareDiffusionHandles_(InteractionVolume& iv,
                                  DataHandle& handle,
                                  bool forceUpdateAll)
    {
        for (unsigned int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++phaseIdx)
        {
            for (unsigned int compIdx = 0; compIdx < ModelTraits::numFluidComponents(); ++compIdx)
            {
                // skip main component
                using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
                if constexpr (!FluidSystem::isTracerFluidSystem())
                    if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                        continue;

                // fill data in the handle
                handle.diffusionHandle().setPhaseIndex(phaseIdx);
                handle.diffusionHandle().setComponentIndex(compIdx);

                using DiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;

                // get instance of the interaction volume-local assembler
                using Traits = typename InteractionVolume::Traits;
                using IvLocalAssembler = typename Traits::template LocalAssembler<Problem, FVElementGeometry, ElementVolumeVariables>;
                IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

                // maybe (re-)assemble matrices
                if (forceUpdateAll || diffusionIsSolDependent)
                {
                    // lambda to obtain diffusion coefficient
                    const auto getD = [&](const auto& volVars)
                    {
                        if constexpr (FluidSystem::isTracerFluidSystem())
                            return volVars.effectiveDiffusionCoefficient(0, 0, compIdx);
                        else
                            return volVars.effectiveDiffusionCoefficient(phaseIdx, FluidSystem::getMainComponent(phaseIdx), compIdx);
                    };

                    // Effective diffusion coefficients might get zero if saturation = 0.
                    // Compute epsilon to detect obsolete rows in the iv-local matrices during assembly
                    const auto& scv = *scvs(fvGeometry()).begin();
                    const auto& scvf = *scvfs(fvGeometry()).begin();
                    const auto& vv = elemVolVars()[scv];
                    const auto D = [&] ()
                    {
                        // diffusion coefficients below 1e-20 are treated as zeroes!!
                        using std::max;
                        if constexpr (!FluidSystem::isTracerFluidSystem())
                            return max(1e-20, vv.diffusionCoefficient(phaseIdx, FluidSystem::getMainComponent(phaseIdx), compIdx));
                        else
                            return max(1e-20, vv.diffusionCoefficient(0, 0, compIdx));
                    } ();

                    auto eps = 1e-7*computeTpfaTransmissibility(scvf, scv, D, vv.extrusionFactor())*Extrusion::area(scvf);
                    localAssembler.assembleMatrices(handle.diffusionHandle(), iv, getD, eps);
                }

                // assemble vector of mole fractions
                auto getMassOrMoleFraction = [phaseIdx, compIdx] (const auto& volVars)
                {
                    return (DiffusionType::referenceSystemFormulation() == ReferenceSystemFormulation::massAveraged) ? volVars.massFraction(phaseIdx, compIdx) :
                                                                                                                       volVars.moleFraction(phaseIdx, compIdx);
                };

                localAssembler.assembleU(handle.diffusionHandle(), iv, getMassOrMoleFraction);
            }
        }
    }

    //! prepares the quantities necessary for conductive fluxes in the handle
    template<class InteractionVolume, class DataHandle>
    void prepareHeatConductionHandle_(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll)
    {
        // get instance of the interaction volume-local assembler
        using Traits = typename InteractionVolume::Traits;
        using IvLocalAssembler = typename Traits::template LocalAssembler<Problem, FVElementGeometry, ElementVolumeVariables>;
        IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

        // lambda to obtain the effective thermal conductivity
        auto getLambda = [] (const auto& volVars) { return volVars.effectiveThermalConductivity(); };

        // maybe (re-)assemble matrices
        if (forceUpdateAll || heatConductionIsSolDependent)
            localAssembler.assembleMatrices(handle.heatConductionHandle(), iv, getLambda);

        // assemble vector of temperatures
        auto getTemperature = [] (const auto& volVars) { return volVars.temperature(); };
        localAssembler.assembleU(handle.heatConductionHandle(), iv, getTemperature);
    }

    const Problem* problemPtr_;
    const Element* elementPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;

    // We store pointers to an inner and a boundary interaction volume.
    // These are updated during the filling of the caches and the
    // physics-related caches have access to them
    PrimaryInteractionVolume* primaryIv_;
    SecondaryInteractionVolume* secondaryIv_;

    // pointer to the current interaction volume data handle
    PrimaryDataHandle* primaryIvDataHandle_;
    SecondaryDataHandle* secondaryIvDataHandle_;

    // We do an interaction volume-wise filling of the caches
    // While filling, we store a pointer to the current localScvf
    // face data object of the IV so that the individual caches
    // can access it and don't have to retrieve it again
    const PrimaryLocalFaceData* primaryLocalFaceData_;
    const SecondaryLocalFaceData* secondaryLocalFaceData_;
};

} // end namespace Dumux

#endif
