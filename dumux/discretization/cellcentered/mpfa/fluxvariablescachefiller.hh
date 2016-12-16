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
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, false, false>
               : public CCMpfaAdvectionCacheFiller<TypeTag>
{
    using AdvectionFiller = CCMpfaAdvectionCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using Element = typename GridView::template Codim<0>::Entity;

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static const bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);

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
        // forward to the filler for the advective quantities
        AdvectionFiller::fillCaches(problem,
                                    element,
                                    fvGeometry,
                                    elemVolVars,
                                    scvf,
                                    fluxVarsCacheContainer);
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
        // Do basically nothing if advection is solution-independent. Although we have
        // to set the update status to true here as it has been set to false before.
        // This is for compatibility reasons with compositional models.
        if (!solDependentAdvection && useTpfaBoundary)
            fluxVarsCacheContainer[scvf.index()].setUpdateStatus(true);
        // If we don't use tpfa on the boundaries we have to update the whole thing anyway,
        // as we need the local matrices to assemble the neumann fluxes (which could be solution-dependent)
        else
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
    }
};

//! Implementation for problems considering advection & diffusion
template<class TypeTag>
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, true, false>
               : public CCMpfaAdvectionCacheFiller<TypeTag>,
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

    using Element = typename GridView::template Codim<0>::Entity;

    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static const bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static const bool solDependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);

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
        AdvectionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
        DiffusionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
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
        // Do basically nothing if the parameters are solution-independent. Although we have
        // to set the update status to true here as it has been set to false before.
        // This is for compatibility reasons with compositional models.

        // TODO: How to treat !useTpfaBoundary???
        if (!solDependentAdvection && !solDependentDiffusion)
        {
            if (useTpfaBoundary)
                fluxVarsCacheContainer[scvf.index()].setUpdateStatus(true);
            else
                DiffusionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);

            // we're done here
            return;
        }

        if (solDependentAdvection)
            AdvectionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);

        if (solDependentDiffusion || !useTpfaBoundary)
            DiffusionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
    }
};

//! Implementation for problems considering advection & heat conduction
template<class TypeTag>
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, false, true>
               : public CCMpfaAdvectionCacheFiller<TypeTag>,
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
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using Element = typename GridView::template Codim<0>::Entity;

    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static const bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static const bool solDependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);

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
        AdvectionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
        HeatConductionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
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
        // Do basically nothing if the parameters are solution-independent. Although we have
        // to set the update status to true here as it has been set to false before.
        // This is for compatibility reasons with compositional models.

        // TODO: How to treat !useTpfaBoundary???
        if (!solDependentAdvection && !solDependentHeatConduction)
        {
            if (useTpfaBoundary)
                fluxVarsCacheContainer[scvf.index()].setUpdateStatus(true);
            else
                fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);

            // we're done here
            return;
        }

        if (solDependentAdvection || !useTpfaBoundary)
            AdvectionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);

        if (solDependentHeatConduction || !useTpfaBoundary)
            HeatConductionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
    }
};

//! Implementation for problems considering advection, diffusion & heat conduction
template<class TypeTag>
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, true, true>
               : public CCMpfaAdvectionCacheFiller<TypeTag>,
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

    using Element = typename GridView::template Codim<0>::Entity;

    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static const bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static const bool solDependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static const bool solDependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);

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
        AdvectionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
        DiffusionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
        HeatConductionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
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
        // Do basically nothing if the parameters are solution-independent. Although we have
        // to set the update status to true here as it has been set to false before.
        // This is for compatibility reasons with compositional models.

        // TODO: How to treat !useTpfaBoundary???
        if (!solDependentAdvection && !solDependentDiffusion && !solDependentHeatConduction)
        {
            if (useTpfaBoundary)
                fluxVarsCacheContainer[scvf.index()].setUpdateStatus(true);
            else
            {
                DiffusionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
                HeatConductionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
            }

            // we're done here
            return;
        }

        if (solDependentAdvection)
            AdvectionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);

        if (solDependentDiffusion || !useTpfaBoundary)
            DiffusionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);

        if (solDependentHeatConduction || !useTpfaBoundary)
            HeatConductionFiller::fillCaches(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCacheContainer);
    }
};

} // end namespace

#endif
