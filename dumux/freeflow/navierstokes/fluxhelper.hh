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

#ifndef DUMUX_NAVIERSTOKES_FLUXHELPER_HH
#define DUMUX_NAVIERSTOKES_FLUXHELPER_HH

#include <dumux/common/math.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

namespace Dumux {
namespace NavierStokes {

template<class Traits>
class UpwindFluxHelper
{
    static constexpr bool enableEneryBalance = Traits::enableEnergyBalance();
    static constexpr bool isCompositional = (Traits::numFluidComponents() > 1);

public:
    template <class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static Scalar advectiveEnergyFlux(const VolumeVariables& insideVolVars,
                                      const VolumeVariables& outsideVolVars,
                                      const SubControlVolumeFace& scvf,
                                      const Scalar velocity,
                                      const Scalar upwindWeight)
    {
        auto upwindTerm = [](const auto& volVars) { return volVars.density() * volVars.enthalpy(); };
        return advectiveUpwindFlux(insideVolVars,
                                   outsideVolVars,
                                   scvf,
                                   velocity,
                                   upwindWeight,
                                   upwindTerm);
    }

    template <bool enable = isCompositional, std::enable_if_t<!enable, int> = 0,
              class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static Scalar massFlux(const VolumeVariables& insideVolVars,
                           const VolumeVariables& outsideVolVars,
                           const SubControlVolumeFace& scvf,
                           const Scalar velocity,
                           const Scalar upwindWeight)
    {
        auto upwindTerm = [](const auto& volVars) { return volVars.density(); };
        return advectiveUpwindFlux(insideVolVars,
                                   outsideVolVars,
                                   scvf,
                                   velocity,
                                   upwindWeight,
                                   upwindTerm);
    }

    template <bool enable = isCompositional, std::enable_if_t<enable, int> = 0,
              class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static Scalar advectiveComponentFlux(const VolumeVariables& insideVolVars,
                                         const VolumeVariables& outsideVolVars,
                                         const SubControlVolumeFace& scvf,
                                         const Scalar velocity,
                                         const Scalar upwindWeight,
                                         const int compIdx)
    {
        static constexpr bool useMoles = Traits::useMoles();

        auto upwindTerm = [compIdx](const auto& volVars)
        {
            const auto density = useMoles ? volVars.molarDensity() : volVars.density();
            const auto fraction = useMoles ? volVars.moleFraction(compIdx) : volVars.massFraction(compIdx);
            return density * fraction;
        };

        return advectiveUpwindFlux(insideVolVars,
                                   outsideVolVars,
                                   scvf,
                                   velocity,
                                   upwindWeight,
                                   upwindTerm);
    }

    template<class VolumeVariables, class SubControlVolumeFace, class Scalar, class UpwindTerm>
    static Scalar advectiveUpwindFlux(const VolumeVariables& insideVolVars,
                                      const VolumeVariables& outsideVolVars,
                                      const SubControlVolumeFace& scvf,
                                      const Scalar velocity,
                                      const Scalar upwindWeight,
                                      UpwindTerm upwindTerm)
    {
        const bool insideIsUpstream = scvf.directionSign() == sign(velocity);

        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? outsideVolVars : insideVolVars;

        const Scalar flux = (upwindWeight * upwindTerm(upstreamVolVars) +
                            (1.0 - upwindWeight) * upwindTerm(downstreamVolVars))
                            * velocity * scvf.directionSign();

        return flux;
    }
};

template<class Traits, class Vector, int offset = 0>
class BoundaryFluxHelper
{
    static constexpr bool enableEneryBalance = Traits::enableEnergyBalance();
    static constexpr bool isCompositional = (Traits::numFluidComponents() > 1);

    using FluxHelper = UpwindFluxHelper<Traits>;

public:
    template<class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static Vector outflowFlux(const VolumeVariables& insideVolVars,
                              const VolumeVariables& outsideVolVars,
                              const SubControlVolumeFace& scvf,
                              const Scalar velocity,
                              const Scalar upwindWeight = 1.0)
    {
        Vector flux(0.0);
        advectiveFlux(flux, insideVolVars, outsideVolVars, scvf, velocity, upwindWeight);

        return flux;
    }

    template<class Problem, class Element, class FVElementGeometry, class VolumeVariables, class Scalar>
    static Vector outflowFlux(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& insideVolVars,
                              typename VolumeVariables::PrimaryVariables&& boundaryPriVars,
                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                              const Scalar velocity,
                              const Scalar upwindWeight = 1.0)
    {
        Vector flux(0.0);
        VolumeVariables boundaryVolVars;
        boundaryVolVars.update(elementSolution<FVElementGeometry>(std::move(boundaryPriVars)),
                               problem,
                               element,
                               fvGeometry.scv(scvf.insideScvIdx()));

        advectiveFlux(flux, insideVolVars, boundaryVolVars, scvf, velocity, upwindWeight);

        return flux;
    }

    template<class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static void advectiveFlux(Vector& flux,
                              const VolumeVariables& insideVolVars,
                              const VolumeVariables& outsideVolVars,
                              const SubControlVolumeFace& scvf,
                              const Scalar velocity,
                              const Scalar upwindWeight)
    {
        if constexpr (!isCompositional)
            flux[Traits::Indices::conti0EqIdx - offset] = FluxHelper::massFlux(insideVolVars, outsideVolVars, scvf, velocity, upwindWeight);
        else
        {
            static constexpr auto numComponents = Traits::numFluidComponents();

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                flux[Traits::Indices::conti0EqIdx + compIdx - offset] = FluxHelper::advectiveComponentFlux(insideVolVars, outsideVolVars, scvf, velocity, upwindWeight, compIdx);
            }

            // in case one balance is substituted by the total mass balance
            if constexpr (Traits::replaceCompEqIdx() < numComponents)
            {
                flux[Traits::Indices::conti0EqIdx + Traits::replaceCompEqIdx() - offset] = std::accumulate(flux.begin(), flux.end(), 0.0);
            }
        }

        if constexpr (enableEneryBalance)
            flux[Traits::Indices::energyEqIdx - offset] = FluxHelper::advectiveEnergyFlux(insideVolVars, outsideVolVars, scvf, velocity, upwindWeight);
    }
};

} // end namespace NavierStokes
} // end namespace Dumux

#endif
