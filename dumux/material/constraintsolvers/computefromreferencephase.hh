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
 * \ingroup ConstraintSolvers
 * \brief Computes all quantities of a generic fluid state if a
 *        reference phase has been specified.
 *
 * This makes it is possible to specify just one phase and let the
 * remaining ones be calculated by the constraint solver. This
 * constraint solver assumes thermodynamic equilibrium
 */
#ifndef DUMUX_COMPUTE_FROM_REFERENCE_PHASE_HH
#define DUMUX_COMPUTE_FROM_REFERENCE_PHASE_HH

#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup ConstraintSolvers
 * \brief Computes all quantities of a generic fluid state if a
 *        reference phase has been specified.
 *
 * This makes it is possible to specify just one phase and let the
 * remaining ones be calculated by the constraint solver. This
 * constraint solver assumes thermodynamic equilibrium. It assumes the
 * following quantities to be set:
 *
 * - composition (mole+mass fractions) of the *reference* phase \f$x^\kappa_\beta\f$
 * - temperature of the *reference* phase \f$T_\beta\f$
 * - saturations of *all* phases \f$S_\alpha\f$, \f$S_\beta\f$
 * - pressures of *all* phases \f$p_\alpha\f$, \f$p_\beta\f$
 *
 *  \f$ f^\kappa_\beta = f^\kappa_\alpha = \Phi^\kappa_\alpha(\{x^\lambda_\alpha \}, T_\alpha, p_\alpha)  p_\alpha x^\kappa_\alpha\; \f$,
 *
 *  \f$ p_\alpha = p_\beta + p_{c\beta\alpha}\; \f$,
 *
 * after calling the solve() method the following quantities are
 * calculated in addition:
 *
 * - temperature of *all* phases \f$T_\alpha\f$, \f$T_\beta\f$
 * - density, molar density, molar volume of *all* phases
 *   \f$\rho_\alpha\f$, \f$\rho_\beta\f$, \f$\rho_{mol, \alpha}\f$, \f$\rho_{mol, \beta}\f$,
 *   \f$V_{mol, \alpha}\f$, \f$V_{mol, \beta}\f$
 * - composition in mole and mass fractions and molarities of *all* phases
 *   \f$x^\kappa_\alpha\f$, \f$x^\kappa_\beta\f$, \f$X^\kappa_\alpha\f$, \f$X^\kappa_\beta\f$,
 *   \f$c^\kappa_\alpha\f$, \f$c^\kappa_\beta\f$
 * - mean molar masses of *all* phases \f$M_\alpha\f$, \f$M_\beta\f$
 * - fugacity coefficients of *all* components in *all* phases
 *   \f$\Phi^\kappa_\alpha\f$, \f$\Phi^\kappa_\beta\f$
 */
template <class Scalar, class FluidSystem>
class ComputeFromReferencePhase
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    using CompositionFromFugacities = Dumux::CompositionFromFugacities<Scalar, FluidSystem>;
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

public:
    /*!
     * \brief Computes all quantities of a generic fluid state if a
     *        reference phase has been specified.
     *
     * This makes it is possible to specify just one phase and let the
     * remaining ones be calculated by the constraint solver. This
     * constraint solver assumes thermodynamic equilibrium. It assumes the
     * following quantities to be set:
     *
     * - composition (mole+mass fractions) of the *reference* phase
     * - temperature of the *all* phases
     * - saturations of *all* phases
     * - pressures of *all* phases
     *
     * after calling the solve() method the following quantities are
     * calculated in addition:
     *
     * - temperature of *all* phases
     * - density, molar density, molar volume of *all* phases
     * - composition in mole and mass fractions and molaries of *all* phases
     * - mean molar masses of *all* phases
     * - fugacity coefficients of *all* components in *all* phases
     *
     * \param fluidState Thermodynamic state of the fluids
     * \param paramCache  Container for cache parameters
     * \param refPhaseIdx The phase index of the reference phase
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      int refPhaseIdx)
    {
        ComponentVector fugVec;

        // compute the density and enthalpy of the
        // reference phase
        paramCache.updatePhase(fluidState, refPhaseIdx);
        fluidState.setDensity(refPhaseIdx,
                              FluidSystem::density(fluidState,
                                                   paramCache,
                                                   refPhaseIdx));
        fluidState.setMolarDensity(refPhaseIdx,
                                   FluidSystem::molarDensity(fluidState,
                                                             paramCache,
                                                             refPhaseIdx));

        // compute the fugacities of all components in the reference phase
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            fluidState.setFugacityCoefficient(refPhaseIdx,
                                              compIdx,
                                              FluidSystem::fugacityCoefficient(fluidState,
                                                                               paramCache,
                                                                               refPhaseIdx,
                                                                               compIdx));
            fugVec[compIdx] = fluidState.fugacity(refPhaseIdx, compIdx);
        }

        // compute all quantities for the non-reference phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            if (phaseIdx == refPhaseIdx)
                continue; // reference phase is already calculated

            CompositionFromFugacities::guessInitial(fluidState, paramCache, phaseIdx, fugVec);
            CompositionFromFugacities::solve(fluidState, paramCache, phaseIdx, fugVec);
        }
    }
};

} // end namespace Dumux

#endif
