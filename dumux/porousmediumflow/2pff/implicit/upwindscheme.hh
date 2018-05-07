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
 * \brief Base class for the upwind scheme
 */
#ifndef DUMUX_TWOP_FRACFLOW_UPWINDSCHEME_HH
#define DUMUX_TWOP_FRACFLOW_UPWINDSCHEME_HH

#include <type_traits>
#include <dune/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <math.h>

namespace Dumux {

//! Upwind scheme for the fractional flow formulation
class TwoPFractionalFlowUpwindScheme
{
    enum
    {
        wPhaseIdx = 0,
        nPhaseIdx = 1
    };

    enum {
        viscousFluxIdx = 0,
        gravityFluxIdx = 1,
        capillaryFluxIdx = 2
    };

public:
    // applies a simple upwind scheme to the precalculated advective flux
    template<class FluxVariables, class UpwindTermFunction, class Flux>
    static auto apply(const FluxVariables& fluxVars,
                      const UpwindTermFunction& upwindTerm,
                      const Flux& flux, int eqIdx)
    -> std::decay_t<decltype(flux[0])>
    {
        using Scalar = std::decay_t<decltype(flux[0])>;
        using MaterialLaw = typename std::decay_t<decltype(fluxVars.problem().spatialParams())>::MaterialLaw;
        static const bool useHybridUpwinding = getParamFromGroup<bool>(fluxVars.problem().paramGroup(), "Problem.UseHybridUpwinding");
        static const bool enableGravity = getParamFromGroup<bool>(fluxVars.problem().paramGroup(), "Problem.EnableGravity");
        static constexpr int transportEqIdx = 0;

        if (useHybridUpwinding)
        {
            const auto& insideVolVars = fluxVars.elemVolVars()[fluxVars.scvFace().insideScvIdx()];
            const auto& outsideVolVars = fluxVars.elemVolVars()[fluxVars.scvFace().outsideScvIdx()];

            if (eqIdx == transportEqIdx)
            {
                // viscous Flux
                const auto mobW_V = std::signbit(flux[viscousFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx)
                                                                       : insideVolVars.mobility(wPhaseIdx);
                const auto mobN_V = std::signbit(flux[viscousFluxIdx]) ? outsideVolVars.mobility(nPhaseIdx)
                                                                       : insideVolVars.mobility(nPhaseIdx);
                const auto mobT_V = std::max(mobW_V + mobN_V, 1.0e-30);

                const Scalar viscousFlux = mobW_V/mobT_V*flux[viscousFluxIdx];

                // Calculate gravity flux
                Scalar gravFlux = 0.0;
                if (enableGravity)
                {
                    const auto& scvf = fluxVars.scvFace();
                    // do averaging for the density over all neighboring elements
                    const auto rho = [&](int phaseIdx)
                    {
                        // boundaries
                        if (scvf.boundary())
                            return insideVolVars.density(phaseIdx);

                        // inner faces with two neighboring elements
                        else
                            return (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;
                    };

                    if(rho(nPhaseIdx) - rho(wPhaseIdx) < 0)
                    {
                        const auto mobW_G = std::signbit(flux[gravityFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx)
                                                                               : insideVolVars.mobility(wPhaseIdx);
                        const auto mobN_G = std::signbit(flux[gravityFluxIdx]) ? insideVolVars.mobility(nPhaseIdx)
                                                                               : outsideVolVars.mobility(nPhaseIdx);
                        const auto mobT_G = std::max(mobW_G + mobN_G, 1.0e-30);

                        gravFlux = mobW_G*mobN_G/mobT_G * flux[gravityFluxIdx];
                    }
                    else
                    {
                        const auto mobW_G = std::signbit(flux[gravityFluxIdx]) ? insideVolVars.mobility(wPhaseIdx)
                                                                               : outsideVolVars.mobility(wPhaseIdx);
                        const auto mobN_G = std::signbit(flux[gravityFluxIdx]) ? outsideVolVars.mobility(nPhaseIdx)
                                                                               : insideVolVars.mobility(nPhaseIdx);
                        const auto mobT_G = std::max(mobW_G + mobN_G, 1.0e-30);

                        gravFlux = mobW_G*mobN_G/mobT_G * flux[gravityFluxIdx];
                    }

                }

                // Calculate capillary flux
                Scalar S_min = std::min(insideVolVars.saturation(wPhaseIdx),outsideVolVars.saturation(wPhaseIdx));
                Scalar S_max = std::max(insideVolVars.saturation(wPhaseIdx),outsideVolVars.saturation(wPhaseIdx));

                Scalar D_max = 0.0;

                //Here, we assume constant viscosity and constant material laws
                const auto& fvGeometry = fluxVars.fvGeometry();
                const auto& scvf = fluxVars.scvFace();
                const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

                const auto& materialLaws = fluxVars.problem().spatialParams().materialLawParamsAtPos(scv.center());

                // compute D_max
                // TODO use the real maximum of D or a approximation


//naehrungsformel
/*Scalar Sw = 0.0 ;
                Scalar Swmax = 0.37+0.2*log10 (insideVolVars.viscosity(wPhaseIdx)/insideVolVars.viscosity(nPhaseIdx));
                if (S_max < Swmax)
                {
                     Sw = S_max ;
                }
                else
                {
                    if (S_min > Swmax)
                    {
                         Sw =S_min;
                    }
                    else
                    {
                         Sw = Swmax;
                    }
                }

                    Scalar mobW = MaterialLaw::krw(materialLaws, Sw)/insideVolVars.viscosity(wPhaseIdx);
                    Scalar mobN = MaterialLaw::krn(materialLaws, Sw)/insideVolVars.viscosity(nPhaseIdx);
                    Scalar dPc_dSw = MaterialLaw::dpc_dsw(materialLaws, Sw);
                    D_max =  -(mobW*mobN)/(mobW + mobN)*dPc_dSw;

*/


//numerische Methode
///////////////////////////////////////////////
                 unsigned int numIntervals = 10;
                 for(int k=0; k <= numIntervals; k++)
                 {
                     Scalar Sw = S_min + (S_max-S_min)*k/numIntervals;
                     Scalar mobW = MaterialLaw::krw(materialLaws, Sw)/insideVolVars.viscosity(wPhaseIdx);
                     Scalar mobN = MaterialLaw::krn(materialLaws, Sw)/insideVolVars.viscosity(nPhaseIdx);
                     Scalar dPc_dSw = MaterialLaw::dpc_dsw(materialLaws, Sw);
                     D_max = std::max(D_max, -(mobW*mobN)/(mobW + mobN)*dPc_dSw);

                 }
                ///////////////////////////////////////////////



//Mittelwert of Sw_max and Sw_min
/*                   Scalar Sw = (insideVolVars.saturation(wPhaseIdx) + outsideVolVars.saturation(wPhaseIdx))/2.0;
                    Scalar mobW = MaterialLaw::krw(materialLaws, Sw)/insideVolVars.viscosity(wPhaseIdx);
                    Scalar mobN = MaterialLaw::krn(materialLaws, Sw)/insideVolVars.viscosity(nPhaseIdx);
                    Scalar dPc_dSw = MaterialLaw::dpc_dsw(materialLaws, Sw);
                    D_max =  -(mobW*mobN)/(mobW + mobN)*dPc_dSw;

*/

//Exact sw_max
//                Scalar Sw = 0.0 ;
//                if (S_max < 0.5782)
//                {
//                        Sw = S_max ;
//                }
//                else
//                {
//                    Sw = (  S_min > 0.5782) ? S_min :  0.5782;
//                }
//                Scalar mobW = MaterialLaw::krw(materialLaws, Sw)/insideVolVars.viscosity(wPhaseIdx);
//                Scalar mobN = MaterialLaw::krn(materialLaws, Sw)/insideVolVars.viscosity(nPhaseIdx);
//                Scalar dPc_dSw = MaterialLaw::dpc_dsw(materialLaws, Sw);
//                D_max =  -(mobW*mobN)/(mobW + mobN)*dPc_dSw;


                return viscousFlux
                       + D_max*flux[capillaryFluxIdx]
                       + gravFlux;
            }
            else
            {
                DUNE_THROW(Dune::InvalidStateException, "Unknown equation index!");
            }
        } // end hybrid upwinding

        // we use phase potential upwinding as the default
        else
        {
            const auto& insideVolVars = fluxVars.elemVolVars()[fluxVars.scvFace().insideScvIdx()];
            const auto& outsideVolVars = fluxVars.elemVolVars()[fluxVars.scvFace().outsideScvIdx()];

            if (eqIdx == transportEqIdx)
            {
                // viscous Flux
                const auto mobW_V = std::signbit(flux[viscousFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx)
                                                                       : insideVolVars.mobility(wPhaseIdx);
                const auto mobN_V = std::signbit(flux[viscousFluxIdx]) ? outsideVolVars.mobility(nPhaseIdx)
                                                                       : insideVolVars.mobility(nPhaseIdx);
                const auto mobT_V = std::max(mobW_V + mobN_V, 1.0e-30);

                const Scalar viscousFlux = mobW_V/mobT_V*flux[viscousFluxIdx];

                // Calculate gravity flux
                Scalar gravFlux = 0.0;
                if (enableGravity)
                {
                    const auto& scvf = fluxVars.scvFace();
                    // do averaging for the density over all neighboring elements
                    const auto rho = [&](int phaseIdx)
                    {
                        // boundaries
                        if (scvf.boundary())
                            return insideVolVars.density(phaseIdx);

                        // inner faces with two neighboring elements
                        else
                            return (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;
                    };

                    if(rho(nPhaseIdx) - rho(wPhaseIdx) < 0)
                    {
                        const auto mobW_G = std::signbit(flux[gravityFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx)
                                                                               : insideVolVars.mobility(wPhaseIdx);
                        const auto mobN_G = std::signbit(flux[gravityFluxIdx]) ? insideVolVars.mobility(nPhaseIdx)
                                                                               : outsideVolVars.mobility(nPhaseIdx);
                        const auto mobT_G = std::max(mobW_G + mobN_G, 1.0e-30);

                        gravFlux = mobW_G*mobN_G/mobT_G * flux[gravityFluxIdx];
                    }
                    else
                    {
                        const auto mobW_G = std::signbit(flux[gravityFluxIdx]) ? insideVolVars.mobility(wPhaseIdx)
                                                                               : outsideVolVars.mobility(wPhaseIdx);
                        const auto mobN_G = std::signbit(flux[gravityFluxIdx]) ? outsideVolVars.mobility(nPhaseIdx)
                                                                               : insideVolVars.mobility(nPhaseIdx);
                        const auto mobT_G = std::max(mobW_G + mobN_G, 1.0e-30);

                        gravFlux = mobW_G*mobN_G/mobT_G * flux[gravityFluxIdx];
                    }

                }

                // capillary flux
                const auto mobW_C = std::signbit(flux[capillaryFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx)
                                                                         : insideVolVars.mobility(wPhaseIdx);
                const auto mobN_C = std::signbit(flux[capillaryFluxIdx]) ? insideVolVars.mobility(nPhaseIdx)
                                                                         : outsideVolVars.mobility(nPhaseIdx);
                const auto mobT_C = std::max(mobW_C + mobN_C, 1.0e-30);

                const Scalar capillaryFlux = mobW_C*mobN_C/mobT_C*flux[capillaryFluxIdx];

                return    viscousFlux
                        + capillaryFlux
                        + gravFlux;


                // for PPU the flux[capillaryFluxIdx] should contain tij*(pc_i - pc_j)
                // while for IHU it should contain tij*(s_i - s_j)
                // return mobW/mobT*flux[viscousFluxIdx]
                //        + mobW*mobN/mobT*flux[capillaryFluxIdx]
                //        + mobW*mobN/mobT*flux[gravityFluxIdx];
            }
            else
            {
                DUNE_THROW(Dune::InvalidStateException, "Unknown equation index!");
            }
        } // end phase potential upwinding
    }
};

} // end namespace Dumux

#endif
