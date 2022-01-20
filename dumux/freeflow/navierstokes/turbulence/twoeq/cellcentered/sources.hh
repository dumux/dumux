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
 * \ingroup TwoEqModel
 * \copydoc Dumux::TwoEqSources
 */
#ifndef DUMUX_CC_TWOEQ_SOURCES_HH
#define DUMUX_CC_TWOEQ_SOURCES_HH

#include <dune/common/hybridutilities.hh>
#include <dumux/common/properties.hh>

namespace Dumux::TwoEqSources {

template <class Problem, class VolumeVariables, class Element, class SubControlVolume>
auto wilcoxTKESource(const Problem& problem,
                     const VolumeVariables& volVars,
                     const Element& element,
                     const SubControlVolume scv)
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
    static const auto enableKOmegaProductionLimiter = getParamFromGroup<bool>(problem.paramGroup(), "KOmega.EnableProductionLimiter", false);
    using std::min;

    // production
    Scalar productionTerm = 2.0 * volVars.dynamicEddyViscosity() * problem.stressTensorScalarProduct(element, scv);
    if (enableKOmegaProductionLimiter)
    {
        Scalar productionAlternative = 20.0 * volVars.density() * volVars.betaK() * volVars.turbulentKineticEnergy() * volVars.dissipation();
        productionTerm = min(productionTerm, productionAlternative);
    }
    Scalar source = productionTerm;

    // destruction
    source -= volVars.betaK() * volVars.density() * volVars.turbulentKineticEnergy() * volVars.dissipation();

    return source;
}

template <class Problem, class VolumeVariables, class Element, class SubControlVolume>
auto wilcoxDissipationSource(const Problem& problem,
                             const VolumeVariables& volVars,
                             const Element& element,
                             const SubControlVolume scv)
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
    static const auto enableKOmegaProductionLimiter = getParamFromGroup<bool>(problem.paramGroup(), "KOmega.EnableProductionLimiter", false);
    static const auto enableKOmegaCrossDiffusion = getParamFromGroup<bool>(problem.paramGroup(), "KOmega.EnableCrossDiffusion", true);
    using std::min;

    // production
    Scalar productionTerm = 2.0 * volVars.dynamicEddyViscosity() * problem.stressTensorScalarProduct(element, scv);
    if (enableKOmegaProductionLimiter)
    {
        Scalar productionAlternative = 20.0 * volVars.density() * volVars.betaK() * volVars.turbulentKineticEnergy() * volVars.dissipation();
        productionTerm = min(productionTerm, productionAlternative);
    }
    Scalar source = volVars.alpha() * volVars.dissipation() / volVars.turbulentKineticEnergy() * productionTerm;

    // destruction
    source -= volVars.betaOmega() * volVars.density() * volVars.dissipation() * volVars.dissipation();

    // cross-diffusion term
    if (enableKOmegaCrossDiffusion)
    {
        Scalar gradientProduct = 0.0;
        for (unsigned int i = 0; i < volVars.storedTurbulentKineticEnergyGradient().size(); ++i)
            gradientProduct += volVars.storedTurbulentKineticEnergyGradient()[i]
                             * volVars.storedDissipationGradient()[i];
        if (gradientProduct > 0.0)
            source += 0.125 * volVars.density() / volVars.dissipation() * gradientProduct;
    }
    return source;
}

template <class Problem, class VolumeVariables, class Element, class SubControlVolume>
auto shearStressTransportTKESource(const Problem& problem,
                                   const VolumeVariables& volVars,
                                   const Element& element,
                                   const SubControlVolume scv)
{
    // production
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
    Scalar source = 2.0 * volVars.dynamicEddyViscosity() * problem.stressTensorScalarProduct(element, scv);

    // destruction
    source -= volVars.betaStar() * volVars.density() * volVars.turbulentKineticEnergy() * volVars.dissipation();

    return source;
}

template <class Problem, class VolumeVariables, class Element, class SubControlVolume>
auto shearStressTransportDissipationSource(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const Element& element,
                                           const SubControlVolume scv)
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
    // production
    Scalar productionTerm = 2.0 * volVars.dynamicEddyViscosity() * problem.stressTensorScalarProduct(element, scv);
    Scalar source =  volVars.gammaWeighted() / volVars.kinematicEddyViscosity() * productionTerm;

    // destruction
    source -= volVars.betaWeighted() * volVars.density() * volVars.dissipation() * volVars.dissipation();

    // cross-diffusion term
    Scalar gradientProduct = 0.0;
    for (unsigned int i = 0; i < volVars.storedTurbulentKineticEnergyGradient().size(); ++i)
        gradientProduct += volVars.storedTurbulentKineticEnergyGradient()[i]
                         * volVars.storedDissipationGradient()[i];
    source += 2.0 * volVars.density() * (1.0 - volVars.fInner()) * volVars.sigmaOmegaInner() / volVars.dissipation() * gradientProduct;

    return source;
}


} // end namespace Dumux

#endif
