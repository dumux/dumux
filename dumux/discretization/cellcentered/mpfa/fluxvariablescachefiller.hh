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
 * \brief A helper class to fill the flux variable caches used in the flux constitutive laws
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/tensorlambdafactory.hh>
#include <dumux/discretization/cellcentered/mpfa/localassembler.hh>

namespace Dumux
{

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Helper class to fill the flux variables caches within
 *        the interaction volume around a given sub-control volume face.
 */
template<class TypeTag>
class CCMpfaFluxVariablesCacheFiller
{
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using MpfaHelper = typename FVGridGeometry::MpfaHelper;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using PrimaryDataHandle = typename ElementFluxVariablesCache::PrimaryIvDataHandle;
    using PrimaryLocalFaceData = typename PrimaryInteractionVolume::Traits::LocalFaceData;
    using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
    using SecondaryDataHandle = typename ElementFluxVariablesCache::SecondaryIvDataHandle;
    using SecondaryLocalFaceData = typename SecondaryInteractionVolume::Traits::LocalFaceData;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static constexpr bool doAdvection = ModelTraits::enableAdvection();
    static constexpr bool doDiffusion = ModelTraits::enableMolecularDiffusion();
    static constexpr bool doHeatConduction = ModelTraits::enableEnergyBalance();

    static constexpr bool soldependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool soldependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static constexpr bool soldependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);

public:
    //! This cache filler is always solution-dependent, as it updates the
    //! vectors of cell unknowns with which the transmissibilities have to be
    //! multiplied in order to obtain the fluxes.
    static constexpr bool isSolDependent = true;

    //! The constructor. Sets problem pointer.
    CCMpfaFluxVariablesCacheFiller(const Problem& problem) : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheContainer Either the element or global flux variables cache
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables (primary/secondary variables)
     * \param scvf The corresponding sub-control volume face
     * \param forceUpdateAll if true, forces all caches to be updated (even the solution-independent ones)
     */
    template<class FluxVariablesCacheContainer>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
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

        // prepare interaction volume and fill caches of all the scvfs connected to it
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        if (fvGridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
        {
            if (forceUpdateAll)
            {
                // the local index of the interaction volume to be created in its container
                const auto ivIndexInContainer = fluxVarsCacheContainer.secondaryInteractionVolumes().size();

                // prepare the locally cached boundary interaction volume
                const auto& indexSet = fvGridGeometry.gridInteractionVolumeIndexSets().secondaryIndexSet(scvf);
                fluxVarsCacheContainer.secondaryInteractionVolumes().emplace_back();
                secondaryIv_ = &fluxVarsCacheContainer.secondaryInteractionVolumes().back();
                secondaryIv_->setUpLocalScope(indexSet, problem(), fvGeometry);

                // prepare the corresponding data handle
                fluxVarsCacheContainer.secondaryDataHandles().emplace_back();
                secondaryIvDataHandle_ = &fluxVarsCacheContainer.secondaryDataHandles().back();
                secondaryIvDataHandle_->resize(*secondaryIv_);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, *secondaryIv_, *secondaryIvDataHandle_, ivIndexInContainer, true);
            }
            else
            {
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                secondaryIv_ = &fluxVarsCacheContainer.secondaryInteractionVolumes()[ivIndexInContainer];
                secondaryIvDataHandle_ = &fluxVarsCacheContainer.secondaryDataHandles()[ivIndexInContainer];

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, *secondaryIv_, *secondaryIvDataHandle_, ivIndexInContainer);
            }
        }
        else
        {
            if (forceUpdateAll)
            {
                // the local index of the interaction volume to be created in its container
                const auto ivIndexInContainer = fluxVarsCacheContainer.primaryInteractionVolumes().size();

                // prepare the locally cached boundary interaction volume
                const auto& indexSet = fvGridGeometry.gridInteractionVolumeIndexSets().primaryIndexSet(scvf);
                fluxVarsCacheContainer.primaryInteractionVolumes().emplace_back();
                primaryIv_ = &fluxVarsCacheContainer.primaryInteractionVolumes().back();
                primaryIv_->setUpLocalScope(indexSet, problem(), fvGeometry);

                // prepare the corresponding data handle
                fluxVarsCacheContainer.primaryDataHandles().emplace_back();
                primaryIvDataHandle_ = &fluxVarsCacheContainer.primaryDataHandles().back();
                primaryIvDataHandle_->resize(*primaryIv_);

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, *primaryIv_, *primaryIvDataHandle_, ivIndexInContainer, true);
            }
            else
            {
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                primaryIv_ = &fluxVarsCacheContainer.primaryInteractionVolumes()[ivIndexInContainer];
                primaryIvDataHandle_ = &fluxVarsCacheContainer.primaryDataHandles()[ivIndexInContainer];

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheContainer, *primaryIv_, *primaryIvDataHandle_, ivIndexInContainer);
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
    template<class FluxVarsCacheContainer, class InteractionVolume, class DataHandle>
    void fillCachesInInteractionVolume_(FluxVarsCacheContainer& fluxVarsCacheContainer,
                                        InteractionVolume& iv,
                                        DataHandle& handle,
                                        unsigned int ivIndexInContainer,
                                        bool forceUpdateAll = false)
    {
        // determine if secondary interaction volumes are used here
        static constexpr bool isSecondary = MpfaHelper::considerSecondaryIVs()
                                            && std::is_same<InteractionVolume, SecondaryInteractionVolume>::value;

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
            const auto& scvfJ = fvGeometry().scvf(d.globalScvfIndex());
            ivScvfs[i] = &scvfJ;
            ivFluxVarCaches[i] = &fluxVarsCacheContainer[scvfJ];
            ivFluxVarCaches[i]->setIvIndexInContainer(ivIndexInContainer);
            ivFluxVarCaches[i]->setUpdateStatus(true);
            ivFluxVarCaches[i]->setSecondaryIvUsage(isSecondary);
            i++;
        }

        fillAdvection(iv, handle, ivScvfs, ivFluxVarCaches, forceUpdateAll);
        fillDiffusion(iv, handle, ivScvfs, ivFluxVarCaches, forceUpdateAll);
        fillHeatConduction(iv, handle, ivScvfs, ivFluxVarCaches, forceUpdateAll);
    }

    //! fills the advective quantities (enabled advection)
    template< class InteractionVolume,
              class DataHandle,
              bool enableAdvection = doAdvection,
              typename std::enable_if_t<enableAdvection, int> = 0 >
    void fillAdvection(InteractionVolume& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                       bool forceUpdateAll = false)
    {
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
        using AdvectionFiller = typename AdvectionType::Cache::Filler;

        // fill data in the handle
        fillAdvectionHandle(iv, handle, forceUpdateAll);

        // fill advection caches
        for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
        {
            // set pointer to current local face data object
            // ifs are evaluated at compile time and are optimized away
            if (std::is_same<PrimaryInteractionVolume, SecondaryInteractionVolume>::value)
            {
                // we cannot make a disctinction, thus we set both pointers
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
            }
            else if (std::is_same<InteractionVolume, PrimaryInteractionVolume>::value)
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

    //! do nothing if advection is not enabled
    template< class InteractionVolume,
              class DataHandle,
              bool enableAdvection = doAdvection,
              typename std::enable_if_t<!enableAdvection, int> = 0 >
    void fillAdvection(InteractionVolume& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                       bool forceUpdateAll = false)
    {}

    //! fills the diffusive quantities (diffusion enabled)
    template< class InteractionVolume,
              class DataHandle,
              bool enableDiffusion = doDiffusion,
              typename std::enable_if_t<enableDiffusion, int> = 0 >
    void fillDiffusion(InteractionVolume& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                       bool forceUpdateAll = false)
    {
        using DiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
        using DiffusionFiller = typename DiffusionType::Cache::Filler;

        static constexpr int numPhases = ModelTraits::numPhases();
        static constexpr int numComponents = ModelTraits::numComponents();

        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

                // fill data in the handle
                fillDiffusionHandle(iv, handle, forceUpdateAll, phaseIdx, compIdx);

                // fill diffusion caches
                for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
                {
                    // set pointer to current local face data object
                    // ifs are evaluated at compile time and are optimized away
                    if (std::is_same<PrimaryInteractionVolume, SecondaryInteractionVolume>::value)
                    {
                        // we cannot make a disctinction, thus we set both pointers
                        primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                        secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
                    }
                    else if (std::is_same<InteractionVolume, PrimaryInteractionVolume>::value)
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

    //! do nothing if diffusion is not enabled
    template< class InteractionVolume,
              class DataHandle,
              bool enableDiffusion = doDiffusion,
              typename std::enable_if_t<!enableDiffusion, int> = 0 >
    void fillDiffusion(InteractionVolume& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                       bool forceUpdateAll = false)
    {}

    //! fills the quantities related to heat conduction (heat conduction enabled)
    template< class InteractionVolume,
              class DataHandle,
              bool enableHeatConduction = doHeatConduction,
              typename std::enable_if_t<enableHeatConduction, int> = 0 >
    void fillHeatConduction(InteractionVolume& iv,
                            DataHandle& handle,
                            const std::vector<const SubControlVolumeFace*>& ivScvfs,
                            const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                            bool forceUpdateAll = false)
    {
        using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);
        using HeatConductionFiller = typename HeatConductionType::Cache::Filler;

        // prepare data in handle
        fillHeatConductionHandle(iv, handle, forceUpdateAll);

        // fill heat conduction caches
        for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
        {
            // set pointer to current local face data object
            // ifs are evaluated at compile time and are optimized away
            if (std::is_same<PrimaryInteractionVolume, SecondaryInteractionVolume>::value)
            {
                // we cannot make a disctinction, thus we set both pointers
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
            }
            else if (std::is_same<InteractionVolume, PrimaryInteractionVolume>::value)
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

    //! do nothing if heat conduction is disabled
    template< class InteractionVolume,
              class DataHandle,
              bool enableHeatConduction = doHeatConduction,
              typename std::enable_if_t<!enableHeatConduction, int> = 0 >
    void fillHeatConduction(InteractionVolume& iv,
                            DataHandle& handle,
                            const std::vector<const SubControlVolumeFace*>& ivScvfs,
                            const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                            bool forceUpdateAll = false)
    {}

    //! prepares the quantities necessary for advective fluxes in the handle
    template< class InteractionVolume,
              class DataHandle,
              class AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType),
              typename std::enable_if_t<AdvectionType::discMethod == DiscretizationMethod::ccmpfa, int> = 0 >
    void fillAdvectionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll)
    {
        using LambdaFactory = TensorLambdaFactory<DiscretizationMethod::ccmpfa>;

        // get instance of the interaction volume-local assembler
        static constexpr MpfaMethods M = InteractionVolume::MpfaMethod;
        using IvLocalAssembler = InteractionVolumeAssembler< Problem, FVElementGeometry, ElementVolumeVariables, M >;
        IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

        // Use different assembly if gravity is enabled
        static const bool enableGravity = getParamFromGroup<bool>(problem().paramGroup(), "Problem.EnableGravity");

        // Assemble T only if permeability is sol-dependent or if update is forced
        if (forceUpdateAll || soldependentAdvection)
        {
            // distinguish between normal/surface grids (optimized away by compiler)
            if (dim < dimWorld)
            {
                if (enableGravity)
                    localAssembler.assembleWithGravity( handle.advectionTout(),
                                                        handle.advectionT(),
                                                        handle.gravityOutside(),
                                                        handle.gravity(),
                                                        handle.advectionCA(),
                                                        handle.advectionA(),
                                                        iv,
                                                        LambdaFactory::getAdvectionLambda() );

                else
                    localAssembler.assemble( handle.advectionTout(),
                                             handle.advectionT(),
                                             iv,
                                             LambdaFactory::getAdvectionLambda() );
            }

            // normal grids
            else
            {
                if (enableGravity)
                    localAssembler.assembleWithGravity( handle.advectionT(),
                                                        handle.gravity(),
                                                        handle.advectionCA(),
                                                        handle.advectionA(),
                                                        handle.advectionAB(),
                                                        handle.advectionN(),
                                                        iv,
                                                        LambdaFactory::getAdvectionLambda() );
                else
                    localAssembler.assemble( handle.advectionA(),
                                             handle.advectionAB(),
                                             handle.advectionCA(),
                                             handle.advectionT(),
                                             handle.advectionN(),
                                             iv,
                                             LambdaFactory::getAdvectionLambda() );
            }
        }

        // (maybe) only reassemble gravity vector
        else if (enableGravity)
        {
            if (dim == dimWorld)
                localAssembler.assembleGravity( handle.gravity(),
                                                iv,
                                                handle.advectionCA(),
                                                LambdaFactory::getAdvectionLambda() );
            else
                localAssembler.assembleGravity( handle.gravity(),
                                                handle.gravityOutside(),
                                                iv,
                                                handle.advectionCA(),
                                                handle.advectionA(),
                                                LambdaFactory::getAdvectionLambda() );
        }

        // assemble pressure vectors
        for (unsigned int pIdx = 0; pIdx < ModelTraits::numPhases(); ++pIdx)
        {
            const auto& evv = &elemVolVars();
            auto getPressure = [&evv, pIdx] (auto volVarIdx) { return (evv->operator[](volVarIdx)).pressure(pIdx); };
            localAssembler.assemble(handle.pressures(pIdx), iv, getPressure);
        }
    }

    //! prepares the quantities necessary for diffusive fluxes in the handle
    template< class InteractionVolume,
              class DataHandle,
              class DiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType),
              typename std::enable_if_t<DiffusionType::discMethod == DiscretizationMethod::ccmpfa, int> = 0 >
    void fillDiffusionHandle(InteractionVolume& iv,
                             DataHandle& handle,
                             bool forceUpdateAll,
                             int phaseIdx, int compIdx)
    {
        using LambdaFactory = TensorLambdaFactory<DiscretizationMethod::ccmpfa>;

        // get instance of the interaction volume-local assembler
        static constexpr MpfaMethods M = InteractionVolume::MpfaMethod;
        using IvLocalAssembler = InteractionVolumeAssembler< Problem, FVElementGeometry, ElementVolumeVariables, M >;
        IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

        // solve the local system subject to the tensor and update the handle
        handle.setPhaseIndex(phaseIdx);
        handle.setComponentIndex(compIdx);

        // assemble T
        if (forceUpdateAll || soldependentDiffusion)
        {
            if (dim < dimWorld)
                localAssembler.assemble( handle.diffusionTout(),
                                         handle.diffusionT(),
                                         iv,
                                         LambdaFactory::getDiffusionLambda(phaseIdx, compIdx) );
            else
                localAssembler. assemble( handle.diffusionT(),
                                          iv,
                                          LambdaFactory::getDiffusionLambda(phaseIdx, compIdx) );
        }

        // assemble vector of mole fractions
        const auto& evv = &elemVolVars();
        auto getMoleFraction = [&evv, phaseIdx, compIdx] (auto volVarIdx)
                               { return (evv->operator[](volVarIdx)).moleFraction(phaseIdx, compIdx); };
        localAssembler.assemble(handle.moleFractions(phaseIdx, compIdx), iv, getMoleFraction);
    }

    //! prepares the quantities necessary for conductive fluxes in the handle
    template< class InteractionVolume,
              class DataHandle,
              class HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType),
              typename std::enable_if_t<HeatConductionType::discMethod == DiscretizationMethod::ccmpfa, int> = 0 >
    void fillHeatConductionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll)
    {
        using LambdaFactory = TensorLambdaFactory<DiscretizationMethod::ccmpfa>;

        // get instance of the interaction volume-local assembler
        static constexpr MpfaMethods M = InteractionVolume::MpfaMethod;
        using IvLocalAssembler = InteractionVolumeAssembler< Problem, FVElementGeometry, ElementVolumeVariables, M >;
        IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

        if (forceUpdateAll || soldependentAdvection)
        {
            if (dim < dimWorld)
                localAssembler.assemble( handle.heatConductionTout(),
                                         handle.heatConductionT(),
                                         iv,
                                         LambdaFactory::template getHeatConductionLambda<typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel)>() );
            else
                localAssembler.assemble( handle.heatConductionT(),
                                         iv,
                                         LambdaFactory::template getHeatConductionLambda<typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel)>() );
        }

        // assemble vector of temperatures
        const auto& evv = &elemVolVars();
        auto getMoleFraction = [&evv] (auto volVarIdx) { return (evv->operator[](volVarIdx)).temperature(); };
        localAssembler.assemble(handle.temperatures(), iv, getMoleFraction);
    }

    //! fill handle only when advection uses mpfa
    template< class InteractionVolume,
              class DataHandle,
              class AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType),
              typename std::enable_if_t<AdvectionType::discMethod != DiscretizationMethod::ccmpfa, int> = 0 >
    void fillAdvectionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll) {}

    //! fill handle only when diffusion uses mpfa
    template< class InteractionVolume,
              class DataHandle,
              class DiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType),
              typename std::enable_if_t<DiffusionType::discMethod != DiscretizationMethod::ccmpfa, int> = 0 >
    void fillDiffusionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll, int phaseIdx, int compIdx) {}

    //! fill handle only when heat conduction uses mpfa
    template< class InteractionVolume,
              class DataHandle,
              class HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType),
              typename std::enable_if_t<HeatConductionType::discMethod != DiscretizationMethod::ccmpfa, int> = 0 >
    void fillHeatConductionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll) {}

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
