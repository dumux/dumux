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
 * \ingroup NavierStokesModel
 * \brief Navier Stokes scalar boundary flux helper
 */
#ifndef DUMUX_MPNC_INITIALCONDITION_HELPER_HH
#define DUMUX_MPNC_INITIALCONDITION_HELPER_HH

#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
namespace Dumux {

/*!
 * \ingroup MPNCModel
 * \brief A helper function to get the correct initial conditions by updating the fluidstate for the MPNC model
 */
template <class Scalar, class FluidSystem>
struct MPNCInitialConditionHelper
{

    /*!
     * \brief Return the area-specific outflow fluxes for all scalar balance equations.
     *        This should only be used of flow reversal does never occur.
     *        A (deactivable) warning is emitted otherwise.
     */
    template<class FluidState, class ParameterCache>
    static void solveFluidStateForMPNCInitialCondition(FluidState &fluidState,
                                                       ParameterCache &paramCache,
                                                       int refPhaseIdx)
    {
        if (fluidState.saturation(0) < 1.0 && fluidState.saturation(0) > 0)
        {
            // make the fluid state consistent with local thermodynamic
            // equilibrium
            using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem>;

            ParameterCache paramCache;
            MiscibleMultiPhaseComposition::solve(fluidState, paramCache);
        }
        else
        {
            using ComputeFromReferencePhase = ComputeFromReferencePhase<Scalar, FluidSystem>;

            ParameterCache paramCache;
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             refPhaseIdx);
        }
    }
};

} // end namespace Dumux

#endif
