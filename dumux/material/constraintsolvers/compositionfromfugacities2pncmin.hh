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
 *
 * \brief Determines the fluid composition given the component
 *        fugacities and an arbitary equation of state.
 */
#ifndef DUMUX_COMPOSITION_FROM_FUGACITIES_2PNCMIN_HH
#define DUMUX_COMPOSITION_FROM_FUGACITIES_2PNCMIN_HH

#include <dune/common/deprecated.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {
/*!
 * \ingroup ConstraintSolver
 * \brief Calculates the chemical equilibrium from the component
 *        fugacities in a phase.
 */
template <class Scalar, class FluidSystem>
class DUNE_DEPRECATED_MSG("CompositionFromFugacities2pncmin is deprecated. Use CompositionFromFugacities2pnc instead.")
  CompositionFromFugacities2pncmin
  : public CompositionFromFugacities2pnc<Scalar, FluidSystem>
{ };
// template <class Scalar, class FluidSystem>
// class CompositionFromFugacities2pncmin
// {
//     enum {
//             numComponents = FluidSystem::numComponents,
//             numMajorComponents = FluidSystem::numPhases
//          };
//
//     typedef typename FluidSystem::ParameterCache ParameterCache;
//
//
// public:
//     typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
//
//     /*!
//      * \brief Guess an initial value for the composition of the phase.
//      * \param fluidState Thermodynamic state of the fluids
//      * \param paramCache  Container for cache parameters
//      * \param phaseIdx The phase index
//      * \param phasePresence The presence index of the reference phase
//      * \param fugVec fugacity vector of the component
//      */
//     template <class FluidState>
//     static void guessInitial(FluidState &fluidState,
//                              ParameterCache &paramCache,
//                              int phaseIdx,
//                              int phasePresence,
//                              const ComponentVector &fugVec)
//     {
//         if (FluidSystem::isIdealMixture(phaseIdx))
//             return;
//
//         // Pure component fugacities
//         for (int i = 0; i < numComponents; ++ i)
//         {
//             fluidState.setMoleFraction(phaseIdx,i, 1.0/numComponents);
//         }
//     }
//
//     /*!
//      * \brief Calculates the chemical equilibrium from the component
//      *        fugacities in a phase. This constraint solver is developed for drying scenarios where
//      *        salt component is restricted to liquid phase and still for the sake for equilibrium
//      *        calculation some residual salt must be considered in the gas phase. In such cases for
//      *        existence of gas phase only, in the theoretical liquid phase, we set the mole fraction
//      *        of salt to  1e-10.
//      * \param fluidState Thermodynamic state of the fluids
//      * \param paramCache  Container for cache parameters
//      * \param phaseIdx The phase index
//      * \param targetFug target fugacity
//      * \param phasePresence Presence of the phase
//      *
//      * The phase's fugacities must already be set.
//      */
//     template <class FluidState>
//     static void solve(FluidState &fluidState,
//                       ParameterCache &paramCache,
//                       int phaseIdx,
//                       int phasePresence,
//                       const ComponentVector &targetFug)
//     {
//         // use a much more efficient method in case the phase is an
//         // ideal mixture
//         if (FluidSystem::isIdealMixture(phaseIdx))
//         {
//             solveIdealMix_(fluidState, paramCache, phaseIdx, phasePresence, targetFug);
//             return;
//         }
//         else
//             DUNE_THROW(NumericalProblem, "This constraint solver is not tested for non-ideal mixtures: Please refer computefromfugacities.hh for details" );
//     }
//
// protected:
//     // update the phase composition in case the phase is an ideal
//     // mixture, i.e. the component's fugacity coefficients are
//     // independent of the phase's composition.
//     template <class FluidState>
//     static void solveIdealMix_(FluidState &fluidState,
//                                ParameterCache &paramCache,
//                                int phaseIdx,
//                                int phasePresence,
//                                const ComponentVector &fugacities)
//     {
//         for (int i = 0; i < numComponents; ++ i)
//         {
//             Scalar phi = FluidSystem::fugacityCoefficient(fluidState,
//                                                           paramCache,
//                                                           phaseIdx,
//                                                           i);
//             Scalar gamma = phi * fluidState.pressure(phaseIdx);
//             fluidState.setFugacityCoefficient(phaseIdx, i, phi);
//             fluidState.setMoleFraction(phaseIdx, i, fugacities[i]/gamma);
//             // Special situation for drying PM and n-phase only situation.set the mole fraction of salt to  1e-10.
//             if (phaseIdx == 0 && i >= numMajorComponents && phasePresence == 2 /*nPhaseOnly*/)
//                 fluidState.setMoleFraction(phaseIdx, i, 1.0e-10);
//         }
//
//         paramCache.updatePhase(fluidState, phaseIdx);
//
//         Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
//         fluidState.setDensity(phaseIdx, rho);
//         return;
//     }
};

template <class Scalar, class FluidSystem>
class DUNE_DEPRECATED_MSG("compositionFromFugacities2pncmin is deprecated. Use CompositionFromFugacities2pncmin (capital C) instead.")
  compositionFromFugacities2pncmin
  : public CompositionFromFugacities2pncmin<Scalar, FluidSystem>
{ };
} // end namespace Dumux

#endif
