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
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present
 */
#ifndef DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_HH
#define DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>

namespace Dumux {
/*!
 * \ingroup ConstraintSolvers
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present
 *
 * The constraint solver assumes the following quantities to be set:
 *
 * - temperatures of *all* phases
 * - saturations of *all* phases
 * - pressures of *all* phases
 *
 * It also assumes that the mole/mass fractions of all phases sum up
 * to 1. After calling the solve() method the following quantities
 * are calculated in addition:
 *
 * - temperature of *all* phases
 * - density, molar density, molar volume of *all* phases
 * - composition in mole and mass fractions and molarities of *all* phases
 * - mean molar masses of *all* phases
 * - fugacity coefficients of *all* components in *all* phases
 */
template <class Scalar, class FluidSystem>
class MiscibleMultiPhaseComposition
{
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;
    static const int numMajorComponents = FluidSystem::numPhases;

public:
    /*!
     * \brief @copybrief Dumux::MiscibleMultiPhaseComposition
     *
     * This function additionally considers a lowering of the saturation vapor pressure
     * of the wetting phase by the Kelvin equation:
     * \f[
     * p^\textrm{w}_\textrm{sat,Kelvin}
     * = p^\textrm{w}_\textrm{sat}
     *   \exp \left( -\frac{p_\textrm{c}}{\varrho_\textrm{w} R_\textrm{w} T} \right)
     * \f]
     *
     * \param fluidState A container with the current (physical) state of the fluid
     * \param paramCache A container for iterative calculation of fluid composition
     * \param knownPhaseIdx The index of the phase with known properties
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      int knownPhaseIdx = 0)
    {
#ifndef NDEBUG
        // currently this solver can only handle fluid systems which
        // assume ideal mixtures of all fluids. TODO: relax this
        // (requires solving a non-linear system of equations, i.e. using
        // Newton method.)
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));

        }
#endif

        //get the known mole fractions from the fluidState
        //in a 2pnc system the n>2 mole fractions are primary variables and are already set in the fluidstate
        Dune::FieldVector<Scalar, numComponents-numMajorComponents> xKnown(0.0);
        for (int knownCompIdx = 0; knownCompIdx < numComponents-numMajorComponents; ++knownCompIdx)
        {
            xKnown[knownCompIdx] = fluidState.moleFraction(knownPhaseIdx, knownCompIdx + numMajorComponents);
        }

        // compute all fugacity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updatePhase(fluidState, phaseIdx);

            // since we assume ideal mixtures, the fugacity
            // coefficients of the components cannot depend on
            // composition, i.e. the parameters in the cache are valid
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }


        // create the linear system of equations which defines the
        // mole fractions
        Dune::FieldMatrix<Scalar, numComponents*numPhases, numComponents*numPhases> M(0.0);
        Dune::FieldVector<Scalar, numComponents*numPhases> x(0.0);
        Dune::FieldVector<Scalar, numComponents*numPhases> b(0.0);

        // assemble the equations expressing the assumption that the
        // sum of all mole fractions in each phase must be 1
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            int rowIdx = numComponents*(numPhases - 1) + phaseIdx;
            b[rowIdx] = 1.0;

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int colIdx = phaseIdx*numComponents + compIdx;

                M[rowIdx][colIdx] = 1.0;
            }
        }

        // set the additional equations for the numComponents-numMajorComponents
        // this is only relevant if numphases != numcomponents e.g. in a 2pnc system
        // Components, of which the mole fractions are known, set to molefraction(knownCompIdx)=xKnown
        for(int knownCompIdx = 0; knownCompIdx < numComponents-numMajorComponents; ++knownCompIdx)
        {
            int rowIdx = numComponents + numPhases + knownCompIdx;
            int colIdx = knownPhaseIdx*numComponents + knownCompIdx + numMajorComponents;
            M[rowIdx][colIdx] = 1.0;
            b[rowIdx] = xKnown[knownCompIdx];
        }

        // assemble the equations expressing the fact that the
        // fugacities of each component are equal in all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            int col1Idx = compIdx;
            const auto entryPhase0 = fluidState.fugacityCoefficient(0, compIdx)*fluidState.pressure(0);

            for (unsigned int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            {
                int rowIdx = (phaseIdx - 1)*numComponents + compIdx;
                int col2Idx = phaseIdx*numComponents + compIdx;
                M[rowIdx][col1Idx] = entryPhase0;
                M[rowIdx][col2Idx] = -fluidState.fugacityCoefficient(phaseIdx, compIdx)*fluidState.pressure(phaseIdx);
            }
        }

        // preconditioning of M to reduce condition number
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            // Multiply row of main component (Raoult's Law) with 10e-5 (order of magn. of pressure)
            if (compIdx < numMajorComponents)
                M[compIdx] *= 10e-5;

            // Multiply row of sec. components (Henry's Law) with 10e-9 (order of magn. of Henry constant)
            else
                M[compIdx] *= 10e-9;
        }

        // solve for all mole fractions
        try { M.solve(x, b); }
        catch (const Dune::FMatrixError& e) {
            DUNE_THROW(NumericalProblem,
                "MiscibleMultiPhaseComposition: Failed to solve the linear equation system for the phase composition.");
        }

        // set all mole fractions and the the additional quantities in
        // the fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int rowIdx = phaseIdx*numComponents + compIdx;
                fluidState.setMoleFraction(phaseIdx, compIdx, x[rowIdx]);
            }
            paramCache.updateComposition(fluidState, phaseIdx);

            Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, value);

            value = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
            fluidState.setMolarDensity(phaseIdx, value);
        }
    }
};

} // end namespace Dumux

#endif
