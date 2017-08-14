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

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{

//! Upwind scheme for the fractional flow formulation
template<class TypeTag>
class TwoPFractionalFlowUpwindScheme
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    enum {
        viscousFluxIdx = 0,
        gravityFluxIdx = 1,
        capillaryFluxIdx = 2
    };

    enum {
        transportEqIdx = Indices::transportEqIdx
    };

public:
    // applies a simple upwind scheme to the precalculated advective flux
    template<class FluxVariables, class UpwindTermFunction, class Flux>
    static Scalar apply(const FluxVariables& fluxVars,
                        const UpwindTermFunction& upwindTerm,
                        const Flux& flux, int eqIdx)
    {
        static const bool useHybridUpwinding = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, UseHybridUpwinding);

        if (useHybridUpwinding)
        {
            const auto& insideVolVars = fluxVars.elemVolVars()[fluxVars.scvFace().insideScvIdx()];
            const auto& outsideVolVars = fluxVars.elemVolVars()[fluxVars.scvFace().outsideScvIdx()];

            if (eqIdx == transportEqIdx)
            {
                // Calculate viscous flux
                const auto mobW_v = std::signbit(flux[viscousFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx)
                                                                     : insideVolVars.mobility(wPhaseIdx);
                const auto mobN_v = std::signbit(flux[viscousFluxIdx]) ? outsideVolVars.mobility(nPhaseIdx)
                                                                     : insideVolVars.mobility(nPhaseIdx);
                const auto mobT_v = mobW_v + mobN_v;

                const Scalar viscousFlux = mobW_v/mobT_v*flux[viscousFluxIdx];

                // Calculate gravity flux
                const auto mobW_mobN_G =
                        std::signbit(flux[gravityFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx) * insideVolVars.mobility(nPhaseIdx)
                                                           : outsideVolVars.mobility(nPhaseIdx) * insideVolVars.mobility(wPhaseIdx);
                const auto mobT_G =
                        std::signbit(flux[gravityFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx) + insideVolVars.mobility(nPhaseIdx)
                                                           : outsideVolVars.mobility(nPhaseIdx) + insideVolVars.mobility(wPhaseIdx);

                const Scalar gravFlux = mobW_mobN_G/mobT_G * flux[gravityFluxIdx];

                // Calculate capillary flux
                Scalar S_min = std::min(insideVolVars.saturation(wPhaseIdx),outsideVolVars.saturation(wPhaseIdx));
                Scalar S_max = std::max(insideVolVars.saturation(wPhaseIdx),outsideVolVars.saturation(wPhaseIdx));

                Scalar D_max = 0.0;
                //Here, we assume constant viscosity and constant material laws
                //ToDo generalize for non-constant material laws, etc.
                const auto& fvGeometry = fluxVars.fvGeometry();
                const auto& scvf = fluxVars.scvFace();
                const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

                auto materialLaws = fluxVars.problem().spatialParams().materialLawParamsAtPos(scv.center());
                unsigned int numIntervals = 10;
                for(int k=0; k <= numIntervals; k++)
                {
                    Scalar Sw = S_min + (S_max-S_min)*k/numIntervals;
                    Scalar mobW = MaterialLaw::krw(materialLaws, Sw)/insideVolVars.viscosity(wPhaseIdx);
                    Scalar mobN = MaterialLaw::krn(materialLaws, Sw)/insideVolVars.viscosity(nPhaseIdx);
                    Scalar dPc_dSw = MaterialLaw::dpc_dsw(materialLaws, Sw);

                    D_max = std::max(D_max,-(mobW*mobN)/(mobW + mobN)*dPc_dSw);
                }

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
                const auto mobT_V = mobW_V + mobN_V;

                Scalar viscousFlux = mobW_V/mobT_V*flux[viscousFluxIdx];

                // Calculate gravity flux
                const auto mobW_mobN_G =
                        std::signbit(flux[gravityFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx) * insideVolVars.mobility(nPhaseIdx)
                                                           : outsideVolVars.mobility(nPhaseIdx) * insideVolVars.mobility(wPhaseIdx);
                const auto mobT_G =
                        std::signbit(flux[gravityFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx) + insideVolVars.mobility(nPhaseIdx)
                                                           : outsideVolVars.mobility(nPhaseIdx) + insideVolVars.mobility(wPhaseIdx);

                const Scalar gravFlux = mobW_mobN_G/mobT_G * flux[gravityFluxIdx];

                // capillary Flux
                const auto mobW_C = std::signbit(flux[capillaryFluxIdx]) ? outsideVolVars.mobility(wPhaseIdx)
                                                                     : insideVolVars.mobility(wPhaseIdx);
                const auto mobN_C = std::signbit(flux[capillaryFluxIdx]) ? outsideVolVars.mobility(nPhaseIdx)
                                                                     : insideVolVars.mobility(nPhaseIdx);
                const auto mobT_C = mobW_C + mobN_C;

                Scalar capillaryFlux = mobW_C*mobN_C/mobT_C*flux[capillaryFluxIdx];

                return   viscousFlux
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
