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
 * \ingroup Fluidmatrixinteractions
 * \brief Makes the twophase capillary pressure-saturation relations
 *        available under the M-phase API for material laws.
 * Also use the temperature dependent version of the material laws.
 */
#ifndef DUMUX_MP_2P_OFT_ADAPTER_HH
#define DUMUX_MP_2P_OFT_ADAPTER_HH

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Adapts the interface of the MpNc material law to the standard-Dumux material law.
 * Also use the temperature dependent version of the material laws.
 */
template <int wPhaseIdx, class TwoPLaw>
class TwoPOfTAdapter
{
    enum { nPhaseIdx = (wPhaseIdx == 0)?1:0 };

public:
    using Params = typename TwoPLaw::Params;
    using Scalar = typename Params::Scalar;
    enum { numPhases = 2 };

    /*!
     * \brief The capillary pressure-saturation curve.
     * \param pc Container for capillary pressure in \f$\mathrm{[Pa]}\f$
     * \param params Array of parameters
     * \param fluidState Fluidstate
     * \param wPhaseIdx the phase index of the wetting phase
     */
    template <class pcContainerT, class FluidState>
    static void capillaryPressures(pcContainerT &pc,
                                   const Params &params,
                                   const FluidState &fluidState,
                                   int wPhaseIdx = 0)
    {
        // non-wetting phase gets the capillary pressure added
        pc[nPhaseIdx] = 0;

        // wetting phase does not get anything added
        pc[wPhaseIdx] = - TwoPLaw::pc(params, fluidState.saturation(wPhaseIdx), fluidState.temperature(wPhaseIdx));
    }

    /*!
     * \brief The relative permeability of all phases.
     * \param kr Container for relative permeability
     * \param params Array of parameters
     * \param fluidState Fluidstate
     */
    template <class krContainerT, class FluidState>
    static void relativePermeabilities(krContainerT &kr,
                   const Params &params,
                   const FluidState &fluidState)
    {
        kr[wPhaseIdx] = TwoPLaw::krw(params, fluidState.saturation(wPhaseIdx));
        kr[nPhaseIdx] = TwoPLaw::krn(params, fluidState.saturation(wPhaseIdx));
    }
};
}

#endif
