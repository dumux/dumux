// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
#include <dumux/discretization/cellcentered/mpfa/interactionvolume.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class PorousMediumFluxVariablesCacheFillerImplementation;

/*!
 * \ingroup PorousmediumflowModels
 * \brief The flux variables cache filler class for porous media
 *
 * Helps filling the flux variables cache depending several policies
 */
template<class TypeTag>
using PorousMediumFluxVariablesCacheFiller = PorousMediumFluxVariablesCacheFillerImplementation<TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

//! Specialization of the flux variables cache filler for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheFillerImplementation<TypeTag, DiscretizationMethods::CCTpfa>
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
class PorousMediumFluxVariablesCacheFillerImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using MpfaHelper = typename GridGeometry::MpfaHelper;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using InteractionVolume = CCMpfaInteractionVolume<GridGeometry, Scalar>;
    using Fluxes = typename InteractionVolume::Fluxes;
    using FluxId = typename Fluxes::FluxId;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static constexpr bool advectionEnabled = ModelTraits::enableAdvection();
    static constexpr bool diffusionEnabled = ModelTraits::enableMolecularDiffusion();
    static constexpr bool heatConductionEnabled = ModelTraits::enableEnergyBalance();

    static constexpr bool advectionUsesMpfa = GetPropType<TypeTag, Properties::AdvectionType>::discMethod == DiscretizationMethods::ccmpfa;
    static constexpr bool diffusionUsesMpfa = GetPropType<TypeTag, Properties::MolecularDiffusionType>::discMethod == DiscretizationMethods::ccmpfa;
    static constexpr bool heatConductionUsesMpfa = GetPropType<TypeTag, Properties::HeatConductionType>::discMethod == DiscretizationMethods::ccmpfa;

    static constexpr bool advectionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentAdvection>();
    static constexpr bool diffusionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentMolecularDiffusion>();
    static constexpr bool heatConductionIsSolDependent = getPropValue<TypeTag, Properties::SolutionDependentHeatConduction>();

    static constexpr unsigned int numFluidPhases = ModelTraits::numFluidPhases();
    static constexpr unsigned int numFluidComponents = ModelTraits::numFluidComponents();

    // These represent the number of tensors/fluxes handled with mpfa
    static constexpr unsigned int numDiffusiveFluxesPerPhase = FluidSystem::isTracerFluidSystem() ? numFluidComponents : numFluidComponents - 1;
    static constexpr unsigned int numAdvectionTensors = advectionUsesMpfa ? 1 : 0;
    static constexpr unsigned int numAdvectionFluxes = advectionUsesMpfa ? numFluidPhases : 0;
    static constexpr unsigned int numDiffusionTensors = diffusionUsesMpfa ? numFluidPhases*numDiffusiveFluxesPerPhase : 0;
    static constexpr unsigned int numDiffusionFluxes = numDiffusionTensors;

    static constexpr unsigned int advectionTensorIdOffset = 0;
    static constexpr unsigned int advectionFluxIdOffset = 0;

    static constexpr unsigned int firstDiffusionTensorIdOffset = numAdvectionTensors;
    static constexpr unsigned int firstDiffusionFluxIdOffset = numAdvectionFluxes;

    static constexpr unsigned int heatConductionTensorIdOffset = numAdvectionTensors + numDiffusionTensors;
    static constexpr unsigned int heatConductionFluxIdOffset = numAdvectionFluxes + numDiffusionFluxes;

public:
    //! Modes of creating/updating the caches
    enum Mode { create, update, forcedUpdate };

    //! This cache filler is always solution-dependent, as it updates the
    //! values of cell unknowns with which the transmissibilities have to be
    //! multiplied in order to obtain the fluxes.
    static constexpr bool isSolDependent = true;

    //! The constructor. Sets problem pointer.
    PorousMediumFluxVariablesCacheFillerImplementation(const Problem& problem)
    : problemPtr_(&problem)
    {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheStorage Class that holds the scvf flux vars caches
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param ivDataStorage Class that stores the interaction volumes & fluxes
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables (primary/secondary variables)
     * \param scvf The corresponding sub-control volume face
     */
    template<class FluxVarsCacheStorage, class FluxVariablesCache, class IVDataStorage>
    void fill(FluxVarsCacheStorage& fluxVarsCacheStorage,
              FluxVariablesCache& scvfFluxVarsCache,
              IVDataStorage& ivDataStorage,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              Mode updateMode = Mode::update)
    {
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
        const auto ivIndexInStorage = [&] () {
            if (updateMode == Mode::create)
            {
                const auto& iv = fvGeometry.gridGeometry().gridInteractionVolumes().get(scvf);
                const auto isDirichletScvf = [pptr=problemPtr_] (const auto& element, const auto& scvf) {
                    return pptr->boundaryTypes(element, scvf).hasOnlyDirichlet();
                };
                ivDataStorage.emplace_back(typename IVDataStorage::value_type{
                    &iv, iv.fluxes(fvGeometry, isDirichletScvf)
                });
                return static_cast<unsigned int>(ivDataStorage.size() - 1);
            }
            else
                return static_cast<unsigned int>(scvfFluxVarsCache.ivIndexInContainer());
        } ();

        fluxes_ = ivDataStorage[ivIndexInStorage].fluxes.get();
        fillCachesInInteractionVolume_(fluxVarsCacheStorage,
                                       ivDataStorage[ivIndexInStorage],
                                       ivIndexInStorage,
                                       updateMode);
    }

    const auto& advectionFluxIds() const { return advectionFluxIds_; }
    const auto& diffusionFluxIds() const { return diffusionFluxIds_; }
    const auto& heatConductionFluxId() const { return heatConductionFluxId_; }

    const Fluxes* fluxesPtr() const
    { return fluxes_; }

private:
    const Problem& problem_() const { return *problemPtr_; }
    const FVElementGeometry& fvGeometry_() const { return *fvGeometryPtr_; }
    const ElementVolumeVariables& elemVolVars_() const { return *elemVolVarsPtr_; }

    //! Method to fill the flux var caches within an interaction volume
    template<class FluxVarsCacheStorage, class IVDataStorage>
    void fillCachesInInteractionVolume_(FluxVarsCacheStorage& fluxVarsCacheStorage,
                                        IVDataStorage& ivDataStorage,
                                        unsigned int ivIndexInContainer,
                                        Mode updateMode)
    {
        if constexpr (advectionEnabled)
            fillAdvection_(fluxVarsCacheStorage, ivDataStorage, ivIndexInContainer, updateMode);
        if constexpr (diffusionEnabled)
            fillDiffusion_(fluxVarsCacheStorage, ivDataStorage, ivIndexInContainer, updateMode);
        if constexpr (heatConductionEnabled)
            fillHeatConduction_(fluxVarsCacheStorage, ivDataStorage, ivIndexInContainer, updateMode);

        // after updating the fluxes in the interaction volume, all embedded scvfs are up-to-date
        ivDataStorage.iv->visitFluxGridScvfIndices([&] (const auto scvfIdx) {
            const auto& scvf = fvGeometry_().scvf(scvfIdx);
            auto& cache = fluxVarsCacheStorage[scvf];
            cache.setIvIndexInContainer(ivIndexInContainer);
            cache.setUpdateStatus(true);
        });
    }

    template<class FluxVarsCacheStorage, class IVDataStorage>
    void fillAdvection_(FluxVarsCacheStorage& fluxVarsCacheStorage,
                        IVDataStorage& ivDataStorage,
                        unsigned int ivIndexInContainer,
                        Mode updateMode)
    {
        static const bool gravity = getParamFromGroup<bool>(problem_().paramGroup(), "Problem.EnableGravity");

        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
        using AdvectionFiller = typename AdvectionType::Cache::Filler;

        updateTransmissibilities_<advectionIsSolDependent>(
            ivDataStorage,
            advectionTensorIdOffset,
            updateMode,
            [evv=elemVolVarsPtr_] (const auto& scv) { return (*evv)[scv].permeability(); }
        );

        for (unsigned int pIdx = 0; pIdx < ModelTraits::numFluidPhases(); ++pIdx)
            advectionFluxIds_[pIdx] = updateValues_(
                ivDataStorage,
                advectionTensorIdOffset,
                advectionFluxIdOffset + pIdx,
                updateMode,
                [pIdx] (const auto& volVars) { return volVars.pressure(pIdx); },
                gravity ? std::optional<typename Fluxes::ForceAccessor>{
                            [pIdx, p=problemPtr_, evv=elemVolVarsPtr_] (const auto& scv) {
                                auto g = p->spatialParams().gravity(scv.center());
                                g *= (*evv)[scv].density(pIdx);
                                return g;
                        }} : std::optional<typename Fluxes::ForceAccessor>{}
            );

        ivDataStorage.iv->visitFluxGridScvfIndices([&] (const auto scvfIdx) {
            const auto& scvf = fvGeometry_().scvf(scvfIdx);
            const auto& element = fvGeometry_().gridGeometry().element(scvf.insideScvIdx());
            AdvectionFiller::fill(fluxVarsCacheStorage[scvf],
                                  problem_(),
                                  element,
                                  fvGeometry_(),
                                  elemVolVars_(),
                                  scvf,
                                  *this);
        });
    }

    template<class FluxVarsCacheStorage, class IVDataStorage>
    void fillDiffusion_(FluxVarsCacheStorage& fluxVarsCacheStorage,
                        IVDataStorage& ivDataStorage,
                        unsigned int ivIndexInContainer,
                        Mode updateMode)
    {
        using DiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
        using DiffusionFiller = typename DiffusionType::Cache::Filler;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        unsigned int idOffset = 0;
        for (unsigned int pIdx = 0; pIdx < ModelTraits::numFluidPhases(); ++pIdx)
            for (unsigned int cIdx = 0; cIdx < ModelTraits::numFluidComponents(); ++cIdx)
            {
                // skip main component
                if constexpr (!FluidSystem::isTracerFluidSystem())
                    if (cIdx == FluidSystem::getMainComponent(pIdx))
                        continue;

                // Effective diffusion coefficients might be zero if saturation = 0.
                // Compute epsilon to detect obsolete rows in the iv-local matrices during assembly
                // TODO: pass/accept zero threshold!!
                static const auto zeroD = getParamFromGroup<Scalar>(
                    problem_().paramGroup(),
                    "Mpfa.ZeroEffectiveDiffusionCoefficientThreshold",
                    1e-16
                );

                updateTransmissibilities_<diffusionIsSolDependent>(
                    ivDataStorage,
                    firstDiffusionTensorIdOffset + idOffset,
                    updateMode,
                    [evv=elemVolVarsPtr_, pIdx, cIdx] (const auto& scv) {
                        if constexpr (FluidSystem::isTracerFluidSystem())
                            return (*evv)[scv].effectiveDiffusionCoefficient(0, 0, cIdx);
                        else
                            return (*evv)[scv].effectiveDiffusionCoefficient(pIdx, FluidSystem::getMainComponent(pIdx), cIdx);
                    }
                );

                diffusionFluxIds_[pIdx][cIdx] = updateValues_(
                    ivDataStorage,
                    firstDiffusionTensorIdOffset + idOffset,
                    firstDiffusionFluxIdOffset + idOffset,
                    updateMode,
                    [pIdx, cIdx] (const auto& volVars) {
                        return (DiffusionType::referenceSystemFormulation() == ReferenceSystemFormulation::massAveraged)
                            ? volVars.massFraction(pIdx, cIdx)
                            : volVars.moleFraction(pIdx, cIdx);
                    }
                );

                ivDataStorage.iv->visitFluxGridScvfIndices([&] (const auto scvfIdx) {
                    const auto& scvf = fvGeometry_().scvf(scvfIdx);
                    const auto& element = fvGeometry_().gridGeometry().element(scvf.insideScvIdx());
                    DiffusionFiller::fill(fluxVarsCacheStorage[scvf],
                                          pIdx,
                                          cIdx,
                                          problem_(),
                                          element,
                                          fvGeometry_(),
                                          elemVolVars_(),
                                          scvf,
                                          *this);
                });

                idOffset++;
            }
    }

    template<class FluxVarsCacheStorage, class IVDataStorage>
    void fillHeatConduction_(FluxVarsCacheStorage& fluxVarsCacheStorage,
                             IVDataStorage& ivDataStorage,
                             unsigned int ivIndexInContainer,
                             Mode updateMode)
    {
        using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;
        using HeatConductionFiller = typename HeatConductionType::Cache::Filler;

        updateTransmissibilities_<heatConductionIsSolDependent>(
            ivDataStorage,
            heatConductionTensorIdOffset,
            updateMode,
            [evv=elemVolVarsPtr_] (const auto& scv) { return (*evv)[scv].effectiveThermalConductivity(); }
        );

        heatConductionFluxId_ = updateValues_(
            ivDataStorage,
            heatConductionTensorIdOffset,
            heatConductionFluxIdOffset,
            updateMode,
            [] (const auto& volVars) { return volVars.temperature(); }
        );

        ivDataStorage.iv->visitFluxGridScvfIndices([&] (const auto scvfIdx) {
            const auto& scvf = fvGeometry_().scvf(scvfIdx);
            const auto& element = fvGeometry_().gridGeometry().element(scvf.insideScvIdx());
            HeatConductionFiller::fill(fluxVarsCacheStorage[scvf],
                                       problem_(),
                                       element,
                                       fvGeometry_(),
                                       elemVolVars_(),
                                       scvf,
                                       *this);
        });
    }

    template<bool isSolDependent, class IVDataStorage, class Tensor>
    void updateTransmissibilities_(IVDataStorage& ivDataStorage,
                                   const unsigned int expectedTensorOffset,
                                   const Mode mode,
                                   Tensor&& tensor)
    {
        switch (mode)
        {
            case Mode::create:
            {
                ivDataStorage.tensorIds.push_back(ivDataStorage.fluxes->registerTensor(std::move(tensor)));
                if (ivDataStorage.tensorIds.size() != expectedTensorOffset + 1)
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected number of tensor ids");
                break;
            }
            case Mode::forcedUpdate:
                ivDataStorage.fluxes->updateTensor(ivDataStorage.tensorIds.at(expectedTensorOffset), std::move(tensor));
                break;
            case Mode::update:
            {
                if constexpr (isSolDependent)
                    ivDataStorage.fluxes->updateTensor(ivDataStorage.tensorIds.at(expectedTensorOffset), std::move(tensor));
                break;
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Unknown update mode");
                break;
        }
    }

    template<class IVDataStorage, class Values>
    FluxId updateValues_(IVDataStorage& ivDataStorage,
                         const unsigned int tensorIdOffset,
                         const unsigned int fluxIdOffset,
                         const Mode mode,
                         const Values& volVarValues,
                         const std::optional<typename Fluxes::ForceAccessor>& forces = {})
    {
        auto values = [evv=elemVolVarsPtr_, v=volVarValues] (const auto& scv) { return v((*evv)[scv]); };
        auto boundaryValues = [evv=elemVolVarsPtr_, v=volVarValues] (const auto& scvf) { return v((*evv)[scvf.outsideScvIdx()]); };

        if (mode == Mode::update)
        {
            const auto fluxId = ivDataStorage.fluxIds.at(fluxIdOffset);
            ivDataStorage.fluxes->updateValuesFor(fluxId, std::move(values), std::move(boundaryValues), forces);
            return fluxId;
        }
        else
        {
            ivDataStorage.fluxIds.push_back(
                ivDataStorage.fluxes->registerValuesFor(
                    ivDataStorage.tensorIds.at(tensorIdOffset), std::move(values), std::move(boundaryValues), forces
                )
            );
            if (ivDataStorage.fluxIds.size() != fluxIdOffset + 1)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected flux id");
            return ivDataStorage.fluxIds.back();
        }
    }

    const Problem* problemPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;

    // during filling, we set this pointer s.t. the individual caches can get and store it
    const Fluxes* fluxes_;

    // during filling, keep track of registered flux ids
    std::array<FluxId, numFluidPhases> advectionFluxIds_;
    std::array<std::array<FluxId, numFluidComponents>, numFluidPhases> diffusionFluxIds_;
    FluxId heatConductionFluxId_;
};

} // end namespace Dumux

#endif
