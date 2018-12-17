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
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * Implements a linear saturation-capillary pressure relation for
 * M-phase fluid systems.
 */
#ifndef DUMUX_MP_LINEAR_MATERIAL_HH
#define DUMUX_MP_LINEAR_MATERIAL_HH

#include "mplinearmaterialparams.hh"

#include <algorithm>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * Implements a linear saturation-capillary pressure relation for
 * M-phase fluid systems.
 *
 * \sa MpLinearMaterialParams
 */
template <int numPhasesV, class ScalarT, class ParamsT = MpLinearMaterialParams<numPhasesV, ScalarT> >
class MpLinearMaterial
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;
    enum { numPhases = numPhasesV };

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
     *
     * \param values Container for the return values
     * \param params Array of Parameters
     * \param state The fluid state
     * \param wPhaseIdx The phase index of the wetting phase
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &state,
                                   int wPhaseIdx = 0)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar S = state.saturation(phaseIdx);
            values[phaseIdx] =
                S*params.pcMaxSat(phaseIdx) +
                (1 - S)*params.pcMinSat(phaseIdx);
        }
    }

    /*!
     * \brief The relative permeability of all phases.
     * \param values Container for the return values
     * \param params Array of Parameters
     * \param state The fluid state
     * \param wPhaseIdx the phase index of the wetting phase
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &state,
                                       int wPhaseIdx = 0)
    {
        using std::max;
        using std::min;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            values[phaseIdx] = max(min(state.saturation(phaseIdx),1.0),0.0);
    }
};
} // end namespace Dumux

#endif
