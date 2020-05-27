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

#ifndef DUMUX_RANS_FLUXHELPER_HH
#define DUMUX_RANS_FLUXHELPER_HH

#include <dumux/common/math.hh>
#include <dumux/freeflow/navierstokes/fluxhelper.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

namespace Dumux {
namespace Rans {

template<class Traits, class Vector, int offset = 0>
class TurbulenceFluxHelper
{
    static constexpr TurbulenceModel turbulenceModel = Traits::turbulenceModel();
    using UpwindFluxHelper = NavierStokes::UpwindFluxHelper<Traits>;

public:
    template <class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static void advectiveTurbulenceFlux(Vector& flux,
                                        const VolumeVariables& insideVolVars,
                                        const VolumeVariables& outsideVolVars,
                                        const SubControlVolumeFace& scvf,
                                        const Scalar velocity,
                                        const Scalar upwindWeight)
    {
        if constexpr (turbulenceModel == TurbulenceModel::oneeq)
        {
            // calculate advective flux
            auto upwindTerm = [](const auto& volVars)
            {
                return volVars.viscosityTilde();
            };
            flux[Traits::Indices::viscosityTildeEqIdx - offset] = UpwindFluxHelper::advectiveUpwindFlux(insideVolVars,
                                                                                                        outsideVolVars,
                                                                                                        scvf,
                                                                                                        velocity,
                                                                                                        upwindWeight,
                                                                                                        upwindTerm);
        }
        else if constexpr (turbulenceModel == TurbulenceModel::kepsilon ||
                           turbulenceModel == TurbulenceModel::komega ||
                           turbulenceModel == TurbulenceModel::lowrekepsilon)
        {
            auto upwindTermK = [](const auto& volVars)
            {
                return volVars.turbulentKineticEnergy();
            };
            flux[Traits::Indices::turbulentKineticEnergyEqIdx - offset] = UpwindFluxHelper::advectiveUpwindFlux(insideVolVars,
                                                                                                                outsideVolVars,
                                                                                                                scvf,
                                                                                                                velocity,
                                                                                                                upwindWeight,
                                                                                                                upwindTermK);
            auto upwindTermD = [](const auto& volVars)
            {
                if constexpr (turbulenceModel == TurbulenceModel::lowrekepsilon)
                    return volVars.dissipationTilde();
                else
                    return volVars.dissipation();
            };
            flux[Traits::Indices::dissipationEqIdx - offset] = UpwindFluxHelper::advectiveUpwindFlux(insideVolVars,
                                                                                                     outsideVolVars,
                                                                                                     scvf,
                                                                                                     velocity,
                                                                                                     upwindWeight,
                                                                                                     upwindTermD);
        }
    }
};

template<class Traits, class Vector, int offset = 0>
class BoundaryFluxHelper
{
    using OutflowFluxHelper = NavierStokes::BoundaryFluxHelper<Traits, Vector, offset>;
    using FluxHelper = TurbulenceFluxHelper<Traits, Vector, offset>;

public:
    template<class VolumeVariables, class SubControlVolumeFace, class Scalar>
    static Vector outflowFlux(const VolumeVariables& insideVolVars,
                              const VolumeVariables& outsideVolVars,
                              const SubControlVolumeFace& scvf,
                              const Scalar velocity,
                              const Scalar upwindWeight = 1.0)
    {
        Vector flux(0.0);
        // Account for advective fluxes of Navier-Stokes equations
        OutflowFluxHelper::advectiveFlux(flux, insideVolVars, outsideVolVars, scvf, velocity, upwindWeight);
        // Account for advective turbulence fluxes
        FluxHelper::advectiveTurbulenceFlux(flux, insideVolVars, outsideVolVars, scvf, velocity, upwindWeight);

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

        // Account for advective fluxes of Navier-Stokes equations
        OutflowFluxHelper::advectiveFlux(flux, insideVolVars, boundaryVolVars, scvf, velocity, upwindWeight);
        // Account for advective turbulence fluxes
        FluxHelper::advectiveTurbulenceFlux(flux, insideVolVars, boundaryVolVars, scvf, velocity, upwindWeight);

        return flux;
    }
};

} // end namespace Rans
} // end namespace Dumux

#endif
