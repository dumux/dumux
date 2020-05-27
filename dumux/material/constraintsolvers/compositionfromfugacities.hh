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
 * \brief Determines the fluid composition given the component
 *        fugacities and an arbitary equation of state.
 */
#ifndef DUMUX_COMPOSITION_FROM_FUGACITIES_HH
#define DUMUX_COMPOSITION_FROM_FUGACITIES_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup ConstraintSolvers
 * \brief Calculates the chemical equilibrium from the component
 *        fugacities \f$ f^\kappa \f$ in the phase \f$ \alpha \f$.
 *
 * This constraint solver takes the component fugacity \f$f^\kappa\f$ of
 * of component \f$\kappa\f$, the temperature \f$ T_\alpha \f$, the pressure
 * \f$p_\alpha\f$ and the composition \f$x^\lambda_\alpha\f$ of a phase
 * \f$\alpha\f$ as input and calculates the mole fraction of component
 * \f$\kappa\f$ in that fluid phase \f$x^\kappa_\alpha\f$. This means
 * that the thermodynamic constraints used by this solver are
 *
 * \f$ f^\kappa = \Phi^\kappa_\alpha(\{x^\lambda_\alpha \}, T_\alpha, p_\alpha)  p_\alpha x^\kappa_\alpha\; \f$,
 *
 * where \f${f^\kappa}\f$, \f$ T_\alpha \f$ and \f$ p_\alpha \f$ are fixed values.
 */
template <class Scalar, class FluidSystem>
class CompositionFromFugacities
{
    enum { numComponents = FluidSystem::numComponents };

    using ParameterCache = typename FluidSystem::ParameterCache;

public:
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    /*!
     * \brief Guess an initial value for the composition of the phase.
     * \param fluidState Thermodynamic state of the fluids
     * \param paramCache  Container for cache parameters
     * \param phaseIdx The phase index
     * \param fugVec fugacity vector of the component
     */
    template <class FluidState>
    static void guessInitial(FluidState &fluidState,
                             ParameterCache &paramCache,
                             int phaseIdx,
                             const ComponentVector &fugVec)
    {
        if (!FluidSystem::isIdealMixture(phaseIdx))
        {
            // Pure component fugacities
            for (unsigned int i = 0; i < numComponents; ++ i)
            {
                fluidState.setMoleFraction(phaseIdx,
                                           i,
                                           1.0/numComponents);
            }
        }
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     * \param fluidState Thermodynamic state of the fluids
     * \param paramCache  Container for cache parameters
     * \param phaseIdx The phase index
     * \param targetFug target fugacity
     *
     * The phase's fugacities must already be set.
     */
    template <class FluidState>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      int phaseIdx,
                      const ComponentVector &targetFug)
    {
        // use a much more efficient method in case the phase is an
        // ideal mixture
        if (FluidSystem::isIdealMixture(phaseIdx)) {
            solveIdealMix_(fluidState, paramCache, phaseIdx, targetFug);
            return;
        }

        // save initial composition in case something goes wrong
        Dune::FieldVector<Scalar, numComponents> xInit;
        for (int i = 0; i < numComponents; ++i) {
            xInit[i] = fluidState.moleFraction(phaseIdx, i);
        }

        /////////////////////////
        // Newton method
        /////////////////////////

        // Jacobian matrix
        Dune::FieldMatrix<Scalar, numComponents, numComponents> J;
        // solution, i.e. phase composition
        Dune::FieldVector<Scalar, numComponents> x;
        // right hand side
        Dune::FieldVector<Scalar, numComponents> b;

        paramCache.updatePhase(fluidState, phaseIdx);

        // maximum number of iterations
        const int nMax = 25;
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            linearize_(J, b, fluidState, paramCache, phaseIdx, targetFug);

            // Solve J*x = b
            x = 0;
            try { J.solve(x, b); }
            catch (Dune::FMatrixError &e)
            { throw NumericalProblem(e.what()); }

            // update the fluid composition. b is also used to store
            // the defect for the next iteration.
            Scalar relError = update_(fluidState, paramCache, x, b, phaseIdx, targetFug);

            if (relError < 1e-9) {
                Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);
                fluidState.setMolarDensity(phaseIdx, rhoMolar);

                //std::cout << "num iterations: " << nIdx << "\n";
                return;
            }
        }

        DUNE_THROW(NumericalProblem,
                   "Calculating the " << FluidSystem::phaseName(phaseIdx)
                   << "Phase composition failed. Initial {x} = {"
                   << xInit
                   << "}, {fug_t} = {" << targetFug << "}, p = " << fluidState.pressure(phaseIdx)
                   << ", T = " << fluidState.temperature(phaseIdx));
    }


protected:
    // update the phase composition in case the phase is an ideal
    // mixture, i.e. the component's fugacity coefficients are
    // independent of the phase's composition.
    template <class FluidState>
    static void solveIdealMix_(FluidState &fluidState,
                               ParameterCache &paramCache,
                               int phaseIdx,
                               const ComponentVector &fugacities)
    {
        for (int i = 0; i < numComponents; ++ i) {
            Scalar phi = FluidSystem::fugacityCoefficient(fluidState,
                                                          paramCache,
                                                          phaseIdx,
                                                          i);
            Scalar gamma = phi * fluidState.pressure(phaseIdx);
            fluidState.setFugacityCoefficient(phaseIdx, i, phi);
            fluidState.setMoleFraction(phaseIdx, i, fugacities[i]/gamma);
        }

        paramCache.updatePhase(fluidState, phaseIdx);

        Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, rho);
        Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
        fluidState.setMolarDensity(phaseIdx, rhoMolar);
    }

    template <class FluidState>
    static void linearize_(Dune::FieldMatrix<Scalar, numComponents, numComponents> &J,
                           Dune::FieldVector<Scalar, numComponents> &defect,
                           FluidState &fluidState,
                           ParameterCache &paramCache,
                           int phaseIdx,
                           const ComponentVector &targetFug)
    {
        // reset jacobian
        J = 0;

        // calculate the defect (deviation of the current fugacities
        // from the target fugacities)
        for (int i = 0; i < numComponents; ++ i) {
            Scalar phi = FluidSystem::fugacityCoefficient(fluidState,
                                                          paramCache,
                                                          phaseIdx,
                                                          i);
            Scalar f = phi*fluidState.pressure(phaseIdx)*fluidState.moleFraction(phaseIdx, i);
            fluidState.setFugacityCoefficient(phaseIdx, i, phi);

            defect[i] = targetFug[i] - f;
        }

        // assemble jacobian matrix of the constraints for the composition
        for (int i = 0; i < numComponents; ++ i) {
            const Scalar eps = 1e-11; //std::max(1e-16, std::abs(x_i)*1e-9);

            ////////
            // approximately calculate partial derivatives of the
            // fugacity defect of all components in regard to the mole
            // fraction of the i-th component. This is done via
            // forward differences

            // deviate the mole fraction of the i-th component
            Scalar x_i = fluidState.moleFraction(phaseIdx, i);
            fluidState.setMoleFraction(phaseIdx, i, x_i + eps);
            paramCache.updateSingleMoleFraction(fluidState, phaseIdx, i);

            // compute new defect and derivative for all component
            // fugacities
            for (int j = 0; j < numComponents; ++j) {
                // compute the j-th component's fugacity coefficient ...
                Scalar phi = FluidSystem::fugacityCoefficient(fluidState,
                                                              paramCache,
                                                              phaseIdx,
                                                              j);
                // ... and its fugacity ...
                Scalar f =
                    phi *
                    fluidState.pressure(phaseIdx) *
                    fluidState.moleFraction(phaseIdx, j);
                // as well as the defect for this fugacity
                Scalar defJPlusEps = targetFug[j] - f;

                // use forward differences to calculate the defect's
                // derivative
                J[j][i] = (defJPlusEps - defect[j])/eps;
            }

            // reset composition to original value
            fluidState.setMoleFraction(phaseIdx, i, x_i);
            paramCache.updateSingleMoleFraction(fluidState, phaseIdx, i);

            // end forward differences
            ////////
        }
    }

    template <class FluidState>
    static Scalar update_(FluidState &fluidState,
                          ParameterCache &paramCache,
                          Dune::FieldVector<Scalar, numComponents> &x,
                          Dune::FieldVector<Scalar, numComponents> &b,
                          int phaseIdx,
                          const Dune::FieldVector<Scalar, numComponents> &targetFug)
    {
        // store original composition and calculate relative error
        Dune::FieldVector<Scalar, numComponents> origComp;
        Scalar relError = 0;
        Scalar sumDelta = 0;
        using std::max;
        using std::abs;
        using std::min;
        for (unsigned int i = 0; i < numComponents; ++i)
        {
            origComp[i] = fluidState.moleFraction(phaseIdx, i);
            relError = max(relError, abs(x[i]));
            sumDelta += abs(x[i]);
        }

        // chop update to at most 20% change in composition
        const Scalar maxDelta = 0.2;
        if (sumDelta > maxDelta)
            x /= (sumDelta/maxDelta);

        // change composition
        for (unsigned int i = 0; i < numComponents; ++i)
        {
            Scalar newx = origComp[i] - x[i];

            // only allow negative mole fractions if the target fugacity is negative
            if (targetFug[i] > 0)
                newx = max(0.0, newx);
            // only allow positive mole fractions if the target fugacity is positive
            else if (targetFug[i] < 0)
                newx = min(0.0, newx);
            // if the target fugacity is zero, the mole fraction must also be zero
            else
                newx = 0;

            fluidState.setMoleFraction(phaseIdx, i, newx);
            //sumMoleFrac += abs(newx);
        }

        paramCache.updateComposition(fluidState, phaseIdx);

        return relError;
    }

    template <class FluidState>
    static Scalar calculateDefect_(const FluidState &params,
                                   int phaseIdx,
                                   const ComponentVector &targetFug)
    {
        Scalar result = 0.0;
        using std::abs;
        for (int i = 0; i < numComponents; ++i) {
            // sum of the fugacity defect weighted by the inverse
            // fugacity coefficient
            result += abs(
                (targetFug[i] - params.fugacity(phaseIdx, i))
                /
                params.fugacityCoeff(phaseIdx, i) );
        }
        return result;
    }
};

} // end namespace Dumux

#endif
