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
 * \ingroup GodunovDiscretization
 * \brief A helper class to fill the flux variable caches used in the flux constitutive laws
 */
#ifndef DUMUX_DISCRETIZATION_GODUNOV_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_GODUNOV_FLUXVARSCACHE_FILLER_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

/*!
* \ingroup GodunovDiscretization
* \brief A helper class to fill the flux variable caches used in the flux constitutive laws
*/
template<class TypeTag>
class GodunovFluxVariablesCacheFiller
{
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    using Element = typename GridView::template Codim<0>::Entity;


    static constexpr bool doAdvection = false;
    static constexpr bool doDiffusion = false;
    static constexpr bool doHeatConduction = false;

    /*
    static constexpr bool soldependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool soldependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static constexpr bool soldependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);
    */

public:

    //! The constructor. Sets the problem pointer
    GodunovFluxVariablesCacheFiller(const Problem& problem) : problemPtr_(&problem) {}

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
    template<class FluxVariablesCacheContainer>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              bool forceUpdateAll = false)
    {
        // fill the physics-related quantities of the caches
        fillAdvection(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
        fillDiffusion(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
    }

private:

    const Problem& problem() const
    { return *problemPtr_; }

    //! method to fill the advective quantities
    template<bool advectionEnabled = doAdvection>
    typename std::enable_if<advectionEnabled>::type
    fillAdvection(FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        // not implemented so far
    }

    //! do nothing if advection is not enabled
    template<bool advectionEnabled = doAdvection>
    typename std::enable_if<!advectionEnabled>::type
    fillAdvection(FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {}

    //! method to fill the diffusive quantities
    template<bool diffusionEnabled = doDiffusion>
    typename std::enable_if<diffusionEnabled>::type
    fillDiffusion(FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        // not implemented so far
    }

    //! do nothing if diffusion is not enabled
    template<bool diffusionEnabled = doDiffusion>
    typename std::enable_if<!diffusionEnabled>::type
    fillDiffusion(FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {}

    //! method to fill the quantities related to heat conduction
    template<bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCache& scvfFluxVarsCache,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf)
    {
    }

    //! do nothing if heat conduction is disabled
    template<bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<!heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCache& scvfFluxVarsCache,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf)
    {}

    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
