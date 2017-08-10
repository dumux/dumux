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
        conti0EqIdx = Indices::conti0EqIdx,
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
            // TODO: Implement hybrid upwinding here
            return 0.0;
        } // end hybrid upwinding

        // we use phase potential upwinding as the default
        else
        {
            const auto& insideVolVars = fluxVars.elemVolVars()[fluxVars.scvFace().insideScvIdx()];
            const auto& outsideVolVars = fluxVars.elemVolVars()[fluxVars.scvFace().outsideScvIdx()];

            if (eqIdx == transportEqIdx)
            {

                // get the potential fluxes
                const auto potFluxes = AdvectionType::flux(fluxVars.problem(),
                                                           fluxVars.element(),
                                                           fluxVars.fvGeometry(),
                                                           fluxVars.elemVolVars(),
                                                           fluxVars.scvFace(),
                                                           conti0EqIdx,
                                                           fluxVars.elemFluxVarsCache());

                // Decide which mobilities to use
                const auto mobW = std::signbit(potFluxes[wPhaseIdx]) ? outsideVolVars.mobility(wPhaseIdx)
                                                                     : insideVolVars.mobility(wPhaseIdx);
                const auto mobN = std::signbit(potFluxes[nPhaseIdx]) ? outsideVolVars.mobility(nPhaseIdx)
                                                                     : insideVolVars.mobility(nPhaseIdx);
                const auto mobT = mobW + mobN;




                const Scalar upwindW = std::signbit(potFluxes[wPhaseIdx]) ? upwindTerm(outsideVolVars, wPhaseIdx)
                                                                     : upwindTerm(insideVolVars, wPhaseIdx);

                const Scalar upwindN = std::signbit(potFluxes[nPhaseIdx]) ? upwindTerm(outsideVolVars, nPhaseIdx)
                                                                     : upwindTerm(insideVolVars, nPhaseIdx);


                return mobW/mobT*(upwindW*flux[wPhaseIdx] + upwindN*flux[nPhaseIdx])
                       + mobW*mobN/mobT*flux[capillaryFluxIdx]
                       + mobW*mobN/mobT*flux[gravityFluxIdx];



                // for PPU the flux[capillaryFluxIdx] should contain tij*(pc_i - pc_j)
                // while for IHU it should contain tij*(s_i - s_j)
                // return mobW/mobT*flux[viscousFluxIdx]
                //        + mobW*mobN/mobT*flux[capillaryFluxIdx]
                //        + mobW*mobN/mobT*flux[gravityFluxIdx];
            }

            else if (eqIdx == conti0EqIdx)
            {
                const Scalar upwindW = std::signbit(flux[wPhaseIdx]) ? upwindTerm(outsideVolVars, wPhaseIdx)
                                                                     : upwindTerm(insideVolVars, wPhaseIdx);

                const Scalar upwindN = std::signbit(flux[nPhaseIdx]) ? upwindTerm(outsideVolVars, nPhaseIdx)
                                                                     : upwindTerm(insideVolVars, nPhaseIdx);

                return upwindW*flux[wPhaseIdx] + upwindN*flux[nPhaseIdx];
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
