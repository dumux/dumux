// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief This file contains the data which is required to calculate
 *        diffusive heat fluxes with Fourier's law.
 */
#ifndef DUMUX_TEST_MULTIDOMAIN_PORENETWORK_FLUX_ADVECTION_HH
#define DUMUX_TEST_MULTIDOMAIN_PORENETWORK_FLUX_ADVECTION_HH

#include <array>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkFlux
 * \brief Hagen–Poiseuille-type flux law to describe the advective flux, blocking wetting phase backflow at outlet
 */
template<class ScalarT, class... TransmissibilityLawTypes>
class CreepingFlowBlockingWPhaseOutletBack
{

public:
    //! Export the Scalar type
    using Scalar = ScalarT;

    //! Export the transmissibility law
    using Transmissibility = Detail::Transmissibility<TransmissibilityLawTypes...>;

    /*!
     * \brief Returns the advective flux of a fluid phase
     *        across the given sub-control volume face (corresponding to a pore throat).
     * \note The flux is given in N*m, and can be converted
     *       into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
     *       for the mobility (1/viscosity) or the product of density and mobility, respectively.
     */
    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace, class ElemFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElemFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];
        // calculate the pressure difference
        const Scalar deltaP = insideVolVars.pressure(phaseIdx) - outsideVolVars.pressure(phaseIdx);
        const Scalar transmissibility = fluxVarsCache.transmissibility(phaseIdx);
        using std::isfinite;
        assert(isfinite(transmissibility));
        Scalar volumeFlow = transmissibility*deltaP;

        static const bool blockWOutlet = getParam<bool>("Problem.DisWReflow", false);
        if (blockWOutlet && fvGeometry.gridGeometry().throatLabel(fvGeometry.gridGeometry().elementMapper().index(element)) == 3 && deltaP <= -1e-20 && phaseIdx == 0)
            return 0.0*deltaP;

        static const bool blockNinlet = getParam<bool>("Problem.DisNOutflow", false);
        if (blockNinlet && fvGeometry.gridGeometry().throatLabel(fvGeometry.gridGeometry().elementMapper().index(element)) == 2 && deltaP > 1e-20 && phaseIdx == 1)
            return 0.0*deltaP;

        // add gravity term
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (enableGravity)
        {
            const Scalar rho = 0.5*insideVolVars.density(phaseIdx) + 0.5*outsideVolVars.density(phaseIdx);
            const Scalar g = problem.spatialParams().gravity(scvf.center()) * scvf.unitOuterNormal();

            // The transmissibility is with respect to the effective throat length (potentially dropping the pore body radii).
            // For gravity, we need to consider the total throat length (i.e., the cell-center to cell-center distance).
            // This might cause some inconsistencies TODO: is there a better way?
            volumeFlow += transmissibility * fluxVarsCache.poreToPoreDistance() * rho * g;
        }

        return volumeFlow;
    }

    /*!
     * \brief Returns the throat conductivity
     */
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                            const ElementVolumeVariables& elemVolVars,
                                            const FluxVariablesCache& fluxVarsCache,
                                            const int phaseIdx)
    {
        if constexpr (ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1)
            return Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
        else
        {
            static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 2);

            const auto& spatialParams = problem.spatialParams();
            using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
            const int wPhaseIdx = spatialParams.template wettingPhase<FluidSystem>(element, elemVolVars);
            using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
            const bool invaded = fluxVarsCache.invaded();
            static const double thetaScalingFactor = getParam<double>("Constraint.Problem.ThetaScalingFactor", 1);
            if constexpr (Dumux::Detail::hasProblemThetaFunction<Problem, Element, FVElementGeometry, ElementVolumeVariables, FluxVariablesCache, SubControlVolumeFace>())
            {
                auto theta = problem.theta(element, fvGeometry, elemVolVars, fluxVarsCache, scvf);
                theta = std::clamp(theta, 0.0, thetaScalingFactor);
                if (phaseIdx == wPhaseIdx)
                {
                    const Scalar k1p = Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
                    const Scalar kw = invaded ? Transmissibility::wettingLayerTransmissibility(element, fvGeometry, scvf, fluxVarsCache)
                                              : Transmissibility::entryWettingLayerTransmissibility(element, fvGeometry, scvf, fluxVarsCache);

                    return kw/thetaScalingFactor * theta +  k1p * (1-theta/thetaScalingFactor) ;
                }
                else // non-wetting phase
                {
                    // auto entryKn = Transmissibility::entryNonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache);
                    return  Transmissibility::nonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache)/thetaScalingFactor* theta;
                }
            }
            else
            {
                if (phaseIdx == wPhaseIdx)
                {
                    return invaded ? Transmissibility::wettingLayerTransmissibility(element, fvGeometry, scvf, fluxVarsCache)
                                   : Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
                }
                else // non-wetting phase
                {
                    return invaded ? Transmissibility::nonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache)
                                   : 0.0;
                }
            }
        }
    }

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static std::array<Scalar, 2> calculateTransmissibilities(const Problem& problem,
                                                             const Element& element,
                                                             const FVElementGeometry& fvGeometry,
                                                             const ElementVolumeVariables& elemVolVars,
                                                             const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                             const FluxVariablesCache& fluxVarsCache)
    {
        static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1);
        const Scalar t = calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, 0);
        return {t, -t};
    }
};

} // end namespace Dumux::PoreNetwork

#endif
