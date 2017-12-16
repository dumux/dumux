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
 * \brief Contains the mass conservation part of the volume variables
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_MASS_KINETIC_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_MASS_KINETIC_HH

#include <dumux/material/fluidstates/nonequilibrium.hh>
#include <dumux/porousmediumflow/mpnc/implicit/mass/volumevariables.hh>
#include <dumux/material/constraintsolvers/fluidsystemcomputefromreferencephase.hh>
#include <dumux/material/constraintsolvers/fluidsystemconstraintsolver.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>


namespace Dumux
{
/*!
 * \brief The compositional part of the volume variables if chemical
 *        equilibrium is _not_ assumed
 *
 *        The interface for basing mass transfer on chemical potentials was present in revision 12743
 */
template <class TypeTag>
class MPNCVolumeVariablesMass<TypeTag, /*bool enableKinetic=*/true>
{
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { moleFrac00Idx = Indices::moleFrac00Idx };

    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    //     here,  we need a constraint solver, which gives me the composition of a phase, if the other one is given
    using ConstraintReferencePhaseSolver = FluidSystemComputeFromReferencePhase<Scalar, FluidSystem>;

    //     here,  we need a constraint solver, which gives me the composition of all phases if p,T is given
    using ConstraintSolver = FluidSystemConstraintSolver<Scalar, FluidSystem>;

    //     this is the constraint solver that applies "fugacities"
//    using ConstraintSolver = MiscibleMultiPhaseComposition<Scalar, FluidSystem>;
public:

    /*!
     * \brief Update composition of all phases in the mutable
     *        parameters from the primary variables.
     *
     *        \param actualFluidState Container for all the secondary variables concerning the fluids
     *        \param paramCache Container for cache parameters
     *        \param priVars The primary Variables
     *        \param *hint the volume variables, usable for initial guess of composition
     *        \param problem The problem
     *        \param element The finite element
     *        \param fvGeometry The finite-volume geometry in the fully implicit scheme
     *        \param scvIdx The index of the sub-control volume
     */
    void update(FluidState & actualFluidState,
                ParameterCache & paramCache,
                const PrimaryVariables & priVars,
                const VolumeVariables * hint,
                const Problem & problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const unsigned int scvIdx)
    {
        // setting the mole fractions of the fluid state
        for(int smallLoopPhaseIdx=0; smallLoopPhaseIdx<numPhases; ++smallLoopPhaseIdx){
                // set the component mole fractions
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    actualFluidState.setMoleFraction(smallLoopPhaseIdx,
                           compIdx,
                           priVars[moleFrac00Idx +
                                   smallLoopPhaseIdx*numComponents +
                                   compIdx]);
                }
            }

//            // For using the ... other way of calculating equilibrium
//             THIS IS ONLY FOR silencing Valgrind but is not used in this model
            for(int smallLoopPhaseIdx=0; smallLoopPhaseIdx<numPhases; ++smallLoopPhaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar phi = FluidSystem::fugacityCoefficient(actualFluidState,
                                                                        paramCache,
                                                                        smallLoopPhaseIdx,
                                                                        compIdx);
                    actualFluidState.setFugacityCoefficient(smallLoopPhaseIdx,
                                                      compIdx,
                                                      phi);
            }

            FluidState equilFluidState; // the fluidState *on the interface* i.e. chemical equilibrium
            equilFluidState.assign(actualFluidState) ;
            ConstraintSolver::solve(equilFluidState,
                                    paramCache,
                                    /*setViscosity=*/false,
                                    /*setEnthalpy=*/false) ;

            // Setting the equilibrium composition (in a kinetic model not necessarily the same as the actual mole fraction)
            for(int smallLoopPhaseIdx=0; smallLoopPhaseIdx<numPhases; ++smallLoopPhaseIdx){
                for (int compIdx=0; compIdx< numComponents; ++ compIdx){
                    xEquil_[smallLoopPhaseIdx][compIdx] = equilFluidState.moleFraction(smallLoopPhaseIdx, compIdx);
                }
            }

            // compute densities of all phases
            for(int smallLoopPhaseIdx=0; smallLoopPhaseIdx<numPhases; ++smallLoopPhaseIdx){
                const Scalar rho = FluidSystem::density(actualFluidState, paramCache, smallLoopPhaseIdx);
                actualFluidState.setDensity(smallLoopPhaseIdx, rho);
            }

            // let Valgrind check whether everything is properly defined.
            checkDefined();
        }

    /*!
     * \brief The mole fraction we would have in the case of chemical equilibrium /
     *        on the interface.
     *
     *     \param phaseIdx The index of the fluid phase
     *     \param compIdx The local index of the component
     */
    const Scalar xEquil(const unsigned int phaseIdx, const unsigned int compIdx) const
    {
        return xEquil_[phaseIdx][compIdx] ;
    }


    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if HAVE_VALGRIND && !defined NDEBUG
        Valgrind::CheckDefined(xEquil_);
#endif
    }

private:
    Scalar xEquil_[numPhases][numComponents];
};

} // end namespace

#endif
