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
 * \brief The flux variables cache filler class
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_GLOBAL_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_GLOBAL_FLUXVARSCACHE_FILLER_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include "fluxvariablescachefillerbase.hh"

namespace Dumux
{
//! Forward declaration of the actual implementation
template<class TypeTag, bool advection, bool diffusion, bool energy>
class CCMpfaFluxVariablesCacheFillerImplementation;

/*!
 * \ingroup ImplicitModel
 * \brief Helper class to fill the flux var caches
 */
template<class TypeTag>
using CCMpfaFluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFillerImplementation<TypeTag,
                                                                                    GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                                                    GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                                                    GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

//! Implementation for purely advective problems
template<class TypeTag>
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, false, false> : public CCMpfaAdvectionCacheFiller<TypeTag>
{
    using AdvectionFiller = CCMpfaAdvectionCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer>
    static void fillFluxVarCache(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        // Instantiate interaction volume depending on if scvf touches the boundary or not
        if (problem.model().globalFvGeometry().touchesInteriorOrDomainBoundary(scvf))
        {
            BoundaryInteractionVolume iv(problem.model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf),
                                         problem,
                                         fvGeometry,
                                         elemVolVars);

            // forward to the filler for the advective quantities
            AdvectionFiller::fillCaches(problem,
                                        element,
                                        fvGeometry,
                                        elemVolVars,
                                        scvf,
                                        iv,
                                        fluxVarsCacheContainer);

            //set update status and  maybe update data on interior boundaries
            for (const auto scvfIdxJ : iv.globalScvfs())
            {
                if (enableInteriorBoundaries)
                    fluxVarsCacheContainer[scvfIdxJ].updateInteriorBoundaryData(iv, fvGeometry.scvf(scvfIdxJ));
                fluxVarsCacheContainer[scvfIdxJ].setUpdateStatus(true);
            }
        }
        else
        {
            InteractionVolume iv(problem.model().globalFvGeometry().interactionVolumeSeed(scvf),
                                 problem,
                                 fvGeometry,
                                 elemVolVars);

            // forward to the filler for the advective quantities
            AdvectionFiller::fillCaches(problem,
                                        element,
                                        fvGeometry,
                                        elemVolVars,
                                        scvf,
                                        iv,
                                        fluxVarsCacheContainer);

            // set update status
            for (const auto scvfIdxJ : iv.globalScvfs())
                fluxVarsCacheContainer[scvfIdxJ].setUpdateStatus(true);
        }
    }

    //! function to update the flux var caches during derivative calculation
    template<class FluxVarsCacheContainer>
    static void updateFluxVarCache(const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        //! If advection is solution-independent, only update the caches for scvfs that touch
        //! a boundary as we have to update the boundary contributions (possibly solution-dependent)
        if (!solDependentAdvection)
        {
            const bool touchesBoundary = problem.model().globalFvGeometry().touchesInteriorOrDomainBoundary(scvf);

            //! Do the whole update where a boundary is touched
            if (!useTpfaBoundary && touchesBoundary)
                fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
            //! Simply tell the cache that it is up to date
            else
                fluxVarsCacheContainer[scvf.index()].setUpdateStatus(true);
        }
        else
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
    }
};

//! Implementation for problems considering advection & diffusion
template<class TypeTag>
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, true, false> : public CCMpfaAdvectionCacheFiller<TypeTag>,
                                                                                 public CCMpfaDiffusionCacheFiller<TypeTag>
{
    using AdvectionFiller = CCMpfaAdvectionCacheFiller<TypeTag>;
    using DiffusionFiller = CCMpfaDiffusionCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    static constexpr bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool advectionUsesMpfa = AdvectionType::myDiscretizationMethod == DiscretizationMethods::CCMpfa;

    static constexpr bool solDependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static constexpr bool diffusionUsesMpfa = MolecularDiffusionType::myDiscretizationMethod == DiscretizationMethods::CCMpfa;

    //! Whether or not the caches have to be updated after solution deflection depends on
    //!     - the solution dependency of the parameters
    //!     - if tpfa is used on the boundaries (advection has to be updated because of neumann boundaries)
    //!     - if single-phasic compositional model is used
    //!       (diffusion has to be updated because of neumann boundaries for the transported component)
    static constexpr bool doAdvectionUpdate = (solDependentAdvection || !useTpfaBoundary) && advectionUsesMpfa ;
    static constexpr bool doDiffusionUpdate = (solDependentDiffusion || (numPhases == 1 && !useTpfaBoundary)) && diffusionUsesMpfa;

    using BoolPair = std::pair<bool, bool>;

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer>
    static void fillFluxVarCache(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 FluxVarsCacheContainer& fluxVarsCacheContainer,
                                 const BoolPair doAdvectionOrDiffusion = BoolPair(true, true))
    {
        // Instantiate interaction volume depending on if scvf touches the boundary or not
        if (problem.model().globalFvGeometry().touchesInteriorOrDomainBoundary(scvf))
        {
            BoundaryInteractionVolume iv(problem.model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf),
                                         problem,
                                         fvGeometry,
                                         elemVolVars);

            // forward to the filler for the advective quantities
            if (doAdvectionOrDiffusion.first)
                AdvectionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // forward to the filler for the diffusive quantities
            if (doAdvectionOrDiffusion.second)
                DiffusionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            //set update status and  maybe update data on interior boundaries
            for (const auto scvfIdxJ : iv.globalScvfs())
            {
                if (enableInteriorBoundaries)
                    fluxVarsCacheContainer[scvfIdxJ].updateInteriorBoundaryData(iv, fvGeometry.scvf(scvfIdxJ));
                fluxVarsCacheContainer[scvfIdxJ].setUpdateStatus(true);
            }
        }
        else
        {
            InteractionVolume iv(problem.model().globalFvGeometry().interactionVolumeSeed(scvf),
                                 problem,
                                 fvGeometry,
                                 elemVolVars);

            // forward to the filler for the advective quantities
            if (doAdvectionOrDiffusion.first)
                AdvectionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // forward to the filler for the diffusive quantities
            if (doAdvectionOrDiffusion.second)
                DiffusionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // set update status
            for (const auto scvfIdxJ : iv.globalScvfs())
                fluxVarsCacheContainer[scvfIdxJ].setUpdateStatus(true);
        }
    }

    //! function to update the flux var caches during derivative calculation
    template<class FluxVarsCacheContainer>
    static void updateFluxVarCache(const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        const bool touchesDomainBoundary = problem.model().globalFvGeometry().touchesDomainBoundary(scvf);

        // maybe update Advection or diffusion or both
        const bool updateAdvection = doAdvectionUpdate && touchesDomainBoundary;
        const bool updateDiffusion = doDiffusionUpdate && touchesDomainBoundary;

        if (updateAdvection && updateDiffusion)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolPair(true, true));
        else if (updateAdvection && !updateDiffusion)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolPair(true, false));
        else if (!updateAdvection && updateDiffusion)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolPair(false, true));

        // tell the cache that it has been updated in any case
        fluxVarsCacheContainer[scvf.index()].setUpdateStatus(true);
    }
};

//! Implementation for problems considering advection & heat conduction
template<class TypeTag>
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, false, true> : public CCMpfaAdvectionCacheFiller<TypeTag>,
                                                                                 public CCMpfaHeatConductionCacheFiller<TypeTag>
{
    using AdvectionFiller = CCMpfaAdvectionCacheFiller<TypeTag>;
    using HeatConductionFiller = CCMpfaHeatConductionCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    static constexpr bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool advectionUsesMpfa = AdvectionType::myDiscretizationMethod == DiscretizationMethods::CCMpfa;

    static constexpr bool solDependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);
    static constexpr bool heatConductionUsesMpfa = HeatConductionType::myDiscretizationMethod == DiscretizationMethods::CCMpfa;

    //! Whether or not the caches have to be updated after solution deflection depends on
    //!     - the solution dependency of the parameters
    //!     - if tpfa is used on the boundaries (advection has to be updated because of neumann boundaries)
    static constexpr bool doAdvectionUpdate = (solDependentAdvection || !useTpfaBoundary) && advectionUsesMpfa ;
    static constexpr bool doHeatConductionUpdate = (solDependentHeatConduction || !useTpfaBoundary) && heatConductionUsesMpfa;

    using BoolPair = std::pair<bool, bool>;

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer>
    static void fillFluxVarCache(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 FluxVarsCacheContainer& fluxVarsCacheContainer,
                                 const BoolPair doAdvectionOrHeatConduction = BoolPair(true, true))
    {
        // Instantiate interaction volume depending on if scvf touches the boundary or not
        if (problem.model().globalFvGeometry().touchesInteriorOrDomainBoundary(scvf))
        {
            BoundaryInteractionVolume iv(problem.model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf),
                                         problem,
                                         fvGeometry,
                                         elemVolVars);

            // forward to the filler for the advective quantities
            if (doAdvectionOrHeatConduction.first)
                AdvectionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // forward to the filler for the heat conduction quantities
            if (doAdvectionOrHeatConduction.second)
                HeatConductionFiller::fillCaches(problem,
                                                 element,
                                                 fvGeometry,
                                                 elemVolVars,
                                                 scvf,
                                                 iv,
                                                 fluxVarsCacheContainer);

            //set update status and  maybe update data on interior boundaries
            for (const auto scvfIdxJ : iv.globalScvfs())
            {
                if (enableInteriorBoundaries)
                    fluxVarsCacheContainer[scvfIdxJ].updateInteriorBoundaryData(iv, fvGeometry.scvf(scvfIdxJ));
                fluxVarsCacheContainer[scvfIdxJ].setUpdateStatus(true);
            }
        }
        else
        {
            InteractionVolume iv(problem.model().globalFvGeometry().interactionVolumeSeed(scvf),
                                 problem,
                                 fvGeometry,
                                 elemVolVars);

            // forward to the filler for the advective quantities
            if (doAdvectionOrHeatConduction.first)
                AdvectionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // forward to the filler for the heat conduction quantities
            if (doAdvectionOrHeatConduction.second)
                HeatConductionFiller::fillCaches(problem,
                                                 element,
                                                 fvGeometry,
                                                 elemVolVars,
                                                 scvf,
                                                 iv,
                                                 fluxVarsCacheContainer);

            // set update status
            for (const auto scvfIdxJ : iv.globalScvfs())
                fluxVarsCacheContainer[scvfIdxJ].setUpdateStatus(true);
        }
    }

    //! function to update the flux var caches during derivative calculation
    template<class FluxVarsCacheContainer>
    static void updateFluxVarCache(const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        const bool touchesDomainBoundary = problem.model().globalFvGeometry().touchesDomainBoundary(scvf);

        // maybe update Advection or diffusion or both
        const bool updateAdvection = doAdvectionUpdate && touchesDomainBoundary;
        const bool updateHeatConduction = doHeatConductionUpdate && touchesDomainBoundary;

        if (updateAdvection && updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolPair(true, true));
        else if (updateAdvection && !updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolPair(true, false));
        else if (!updateAdvection && updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolPair(false, true));

        // tell the cache that it has been updated in any case
        fluxVarsCacheContainer[scvf.index()].setUpdateStatus(true);
    }
};

//! Implementation for problems considering advection, diffusion & heat conduction
template<class TypeTag>
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, true, true> : public CCMpfaAdvectionCacheFiller<TypeTag>,
                                                                                public CCMpfaDiffusionCacheFiller<TypeTag>,
                                                                                public CCMpfaHeatConductionCacheFiller<TypeTag>
{
    using AdvectionFiller = CCMpfaAdvectionCacheFiller<TypeTag>;
    using DiffusionFiller = CCMpfaDiffusionCacheFiller<TypeTag>;
    using HeatConductionFiller = CCMpfaHeatConductionCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    static constexpr bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool advectionUsesMpfa = AdvectionType::myDiscretizationMethod == DiscretizationMethods::CCMpfa;

    static constexpr bool solDependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static constexpr bool diffusionUsesMpfa = MolecularDiffusionType::myDiscretizationMethod == DiscretizationMethods::CCMpfa;

    static constexpr bool solDependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);
    static constexpr bool heatConductionUsesMpfa = HeatConductionType::myDiscretizationMethod == DiscretizationMethods::CCMpfa;

    //!     - if single-phasic compositional model is used
    //!       (diffusion has to be updated because of neumann boundaries for the transported component)
    static constexpr bool doAdvectionUpdate = (solDependentAdvection || !useTpfaBoundary) && advectionUsesMpfa ;
    static constexpr bool doDiffusionUpdate = (solDependentDiffusion || (numPhases == 1 && !useTpfaBoundary)) && diffusionUsesMpfa;
    static constexpr bool doHeatConductionUpdate = (solDependentHeatConduction || !useTpfaBoundary) && heatConductionUsesMpfa;

    using BoolTriplet = std::array<bool, 3>;

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer>
    static void fillFluxVarCache(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 FluxVarsCacheContainer& fluxVarsCacheContainer,
                                 const BoolTriplet doAdvectionOrDiffusionOrHeatConduction = BoolTriplet({true, true, true}))
    {
        // Instantiate interaction volume depending on if scvf touches the boundary or not
        if (problem.model().globalFvGeometry().touchesInteriorOrDomainBoundary(scvf))
        {
            BoundaryInteractionVolume iv(problem.model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf),
                                         problem,
                                         fvGeometry,
                                         elemVolVars);

            // forward to the filler for the advective quantities
            if (doAdvectionOrDiffusionOrHeatConduction[0])
                AdvectionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // forward to the filler for the diffusive quantities
            if (doAdvectionOrDiffusionOrHeatConduction[1])
                DiffusionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // forward to the filler for the diffusive quantities
            if (doAdvectionOrDiffusionOrHeatConduction[2])
                HeatConductionFiller::fillCaches(problem,
                                                 element,
                                                 fvGeometry,
                                                 elemVolVars,
                                                 scvf,
                                                 iv,
                                                 fluxVarsCacheContainer);

            //set update status and  maybe update data on interior boundaries
            for (const auto scvfIdxJ : iv.globalScvfs())
            {
                if (enableInteriorBoundaries)
                    fluxVarsCacheContainer[scvfIdxJ].updateInteriorBoundaryData(iv, fvGeometry.scvf(scvfIdxJ));
                fluxVarsCacheContainer[scvfIdxJ].setUpdateStatus(true);
            }
        }
        else
        {
            InteractionVolume iv(problem.model().globalFvGeometry().interactionVolumeSeed(scvf),
                                 problem,
                                 fvGeometry,
                                 elemVolVars);

            // forward to the filler for the advective quantities
            if (doAdvectionOrDiffusionOrHeatConduction[0])
                AdvectionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // forward to the filler for the diffusive quantities
            if (doAdvectionOrDiffusionOrHeatConduction[1])
                DiffusionFiller::fillCaches(problem,
                                            element,
                                            fvGeometry,
                                            elemVolVars,
                                            scvf,
                                            iv,
                                            fluxVarsCacheContainer);

            // forward to the filler for the diffusive quantities
            if (doAdvectionOrDiffusionOrHeatConduction[2])
                HeatConductionFiller::fillCaches(problem,
                                                 element,
                                                 fvGeometry,
                                                 elemVolVars,
                                                 scvf,
                                                 iv,
                                                 fluxVarsCacheContainer);

            // set update status
            for (const auto scvfIdxJ : iv.globalScvfs())
                fluxVarsCacheContainer[scvfIdxJ].setUpdateStatus(true);
        }
    }

    //! function to update the flux var caches during derivative calculation
    template<class FluxVarsCacheContainer>
    static void updateFluxVarCache(const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        const bool touchesDomainBoundary = problem.model().globalFvGeometry().touchesDomainBoundary(scvf);

        // maybe update Advection or diffusion or both
        const bool updateAdvection = doAdvectionUpdate && touchesDomainBoundary;
        const bool updateDiffusion = doDiffusionUpdate && touchesDomainBoundary;
        const bool updateHeatConduction = doHeatConductionUpdate && touchesDomainBoundary;

        if (updateAdvection && updateDiffusion && updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolTriplet({true, true, true}));
        else if (updateAdvection && updateDiffusion && !updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolTriplet({true, true, false}));
        else if (updateAdvection && !updateDiffusion && !updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolTriplet({true, false, false}));
        else if (updateAdvection && !updateDiffusion && updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolTriplet({true, false, true}));
        else if (!updateAdvection && updateDiffusion && updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolTriplet({false, true, true}));
        else if (!updateAdvection && updateDiffusion && !updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolTriplet({false, true, false}));
        else if (!updateAdvection && !updateDiffusion && !updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolTriplet({false, false, false}));
        else if (!updateAdvection && !updateDiffusion && updateHeatConduction)
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer, BoolTriplet({false, false, true}));


        // tell the cache that it has been updated in any case
        fluxVarsCacheContainer[scvf.index()].setUpdateStatus(true);
    }
};

} // end namespace

#endif
