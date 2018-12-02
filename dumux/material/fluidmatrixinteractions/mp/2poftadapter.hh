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
#ifndef DUMUX_MP_OFT_ADAPTER_HH
#define DUMUX_MP_OFT_ADAPTER_HH

#include <algorithm>
#include <cassert>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Adapts the interface of the MpNc material law to the standard-Dumux material law.
 * Also use the temperature dependent version of the material laws.
 */
template <class MaterialLaw, int numPhases>
class MPOfTAdapter
{
    static_assert(AlwaysFalse<MaterialLaw>::value, "Adapter not implemented for the specified number of phases");
};


/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Adapts the interface of the MpNc material law to the standard-Dumux material law.
 * Also use the temperature dependent version of the material laws.
 */
template <class MaterialLaw>
class MPOfTAdapter<MaterialLaw, 2>
{

public:
    using Params = typename MaterialLaw::Params;

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
                                   int wPhaseIdx)
    {
        assert(pc.size() == 2);
        const int nPhaseIdx = 1 - wPhaseIdx;

        // non-wetting phase gets the capillary pressure added
        pc[nPhaseIdx] = 0;

        // wetting phase does not get anything added
        pc[wPhaseIdx] = - MaterialLaw::pc(params, fluidState.saturation(wPhaseIdx), fluidState.temperature(wPhaseIdx));
    }

    /*!
     * \brief The relative permeability of all phases.
     * \param kr Container for relative permeability
     * \param params Array of parameters
     * \param fluidState Fluidstate
     * \param wPhaseIdx the phase index of the wetting phase
     */
    template <class krContainerT, class FluidState>
    static void relativePermeabilities(krContainerT &kr,
                                       const Params &params,
                                       const FluidState &fluidState,
                                       int wPhaseIdx)
    {
        assert(kr.size() == 2);
        const int nPhaseIdx = 1 - wPhaseIdx;
        kr[wPhaseIdx] = MaterialLaw::krw(params, fluidState.saturation(wPhaseIdx));
        kr[nPhaseIdx] = MaterialLaw::krn(params, fluidState.saturation(wPhaseIdx));
    }
};
}

#endif
