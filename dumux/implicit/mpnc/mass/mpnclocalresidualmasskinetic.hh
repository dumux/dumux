/*****************************************************************************
 *   Copyright (C) 2011 by Philipp Nuske
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \brief The local residual for the kinetic mass transfer module of
 *        the compositional multi-phase model.
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_MASS_KINETIC_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_MASS_KINETIC_HH

#include <dumux/implicit/mpnc/mass/mpnclocalresidualmass.hh>

// The idea here is to calculate the mass transfer into both phases (i.e. capture the resistivities of both sides of the interface),
// under the supplementary constrain, that mass fluxes add up to zero. From this the composition of both phases is to be calculated.
// This way the equality of fluxes into both phases does not have to be assumed, but can be calculated.
// NOT YET WORKING!
#if FUNKYMASSTRANSFER
#include <dumux/material/constraintsolvers/compositionfrommasstransfer.hh>
#endif


namespace Dumux
{
/*!
 * \brief The mass conservation part of the Mp-Nc model.
 *
 * This is the specialization for the case where kinetic mass transfer
 * *is* considered.
 */
template<class TypeTag>
class MPNCLocalResidualMass<TypeTag, /*enableKinetic=*/true>
{
    typedef MPNCLocalResidualMassCommon<TypeTag> MassCommon;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

#if FUNKYMASSTRANSFER
    //     here,  we need a constraint solver, which gives me the composition of the phases as determined by the mass transfer constraints
    typedef CompositionFromMassTransfer<Scalar, FluidSystem> NonEquilConstraintSolver;
#endif

    enum { numPhases =  GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { nCompIdx = FluidSystem::nCompIdx } ;
    enum { wCompIdx = FluidSystem::wCompIdx } ;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx} ;
    enum { nPhaseIdx = FluidSystem::nPhaseIdx} ;
    enum { sPhaseIdx = FluidSystem::sPhaseIdx} ;
    enum { enableKineticEnergy = GET_PROP_VALUE(TypeTag, EnableKineticEnergy) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef MPNCLocalResidualEnergy<TypeTag, enableEnergy, enableKineticEnergy> EnergyResid;
    typedef typename Dune::FieldVector<Scalar, numComponents> ComponentVector;

public:

    /*!
     * \brief Calculate the storage for all mass balance equations
     *
     *    \param storage The mass of the component within the sub-control volume
     *    \param volVars The volume variables
     */
    static void computeStorage(PrimaryVariables & storage,
                               const VolumeVariables & volVars)
    {
        for (int phaseIdx =0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                storage[conti0EqIdx + phaseIdx*numComponents + compIdx] = 0.0;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            addPhaseStorage(storage, volVars, phaseIdx);
        }
    }

    /*!
     * \brief Calculate the storage for all mass balance equations
     *        within a single fluid phase
     *
     *     \param storage The mass of the component within the sub-control volume
     *    \param volVars The volume variables
     *    \param phaseIdx phaseIdx The index of the fluid phase
     */
    static void addPhaseStorage(PrimaryVariables & storage,
                                const VolumeVariables & volVars,
                                const unsigned int phaseIdx)
    {
        if (phaseIdx == sPhaseIdx)
            return;

        // calculate the component-wise mass storage
        ComponentVector phaseComponentValues;
        MassCommon::computePhaseStorage(phaseComponentValues,
                                        volVars,
                                        phaseIdx);
        Valgrind::CheckDefined(phaseComponentValues);
        Valgrind::CheckDefined(storage);

        // copy to the primary variables
        for (int compIdx = 0; compIdx < numComponents; ++compIdx){
        	storage[conti0EqIdx + phaseIdx*numComponents + compIdx]
        	        += phaseComponentValues[compIdx];

            #ifndef NDEBUG
            if (!std::isfinite(storage[conti0EqIdx + phaseIdx*numComponents + compIdx]))
            	DUNE_THROW(NumericalProblem, "Calculated non-finite storage");
            #endif
        }

       Valgrind::CheckDefined(storage);
    }

    /*!
     * \brief Calculate the storage for all mass balance equations
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param fluxVars The flux Variables
     *        \param elemVolVars The volume variables of the current element
     */
    static void computeFlux(PrimaryVariables & flux,
                            const FluxVariables & fluxVars,
                            const ElementVolumeVariables & elemVolVars)
    {
        ComponentVector phaseComponentValuesMassTransport[numPhases]; // what goes into the energy module

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            ComponentVector phaseComponentValuesAdvection(0.);
            ComponentVector phaseComponentValuesDiffusion(0.);
            MassCommon::computeAdvectivePhaseFlux(phaseComponentValuesAdvection, fluxVars, phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                flux[conti0EqIdx + phaseIdx*numComponents + compIdx] =
                    phaseComponentValuesAdvection[compIdx];
                Valgrind::CheckDefined(flux);
            }
            MassCommon::computeDiffusivePhaseFlux(phaseComponentValuesDiffusion, fluxVars, phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                flux[conti0EqIdx + phaseIdx*numComponents + compIdx] +=
                        phaseComponentValuesDiffusion[compIdx];
                Valgrind::CheckDefined(flux);

                #ifndef NDEBUG
                if (!std::isfinite(flux[conti0EqIdx + phaseIdx*numComponents + compIdx]))
                	DUNE_THROW(NumericalProblem, "Calculated non-finite flux");
                #endif

            }

            // Right now I think that adding the two contributions individually into the flux is best for debugging and understanding.
            // The Energy module needs both contributions.
            phaseComponentValuesMassTransport[phaseIdx] = phaseComponentValuesDiffusion + phaseComponentValuesAdvection ;
            Valgrind::CheckDefined(flux);

        }// phases

        // The computeflux() of the Energy module needs a component-wise flux (for the diffusive enthalpie transport)
        // It makes some sense calling energy from here, because energy is carried by mass
        // However, it is not really tidied up.
        // todo is there a better way to do this?
        EnergyResid::computeFlux(flux,
                                 fluxVars,
                                 elemVolVars,
                                 phaseComponentValuesMassTransport);
        Valgrind::CheckDefined(flux);
    }


    /*!
     * \brief Calculate the source terms for all mass balance equations
     *
     *         \param source The source/sink in the sub-control volume for each component
     *         \param volVars the volume variables
     */
    static void computeSource(PrimaryVariables & source,
                              const VolumeVariables & volVars)
    {

        // In the case of a kinetic consideration, mass transfer
        // between phases is realized via source terms there is a
        // balance equation for each component in each phase
        ComponentVector componentIntoPhaseMassTransfer[numPhases];

#if FUNKYMASSTRANSFER
        // calculate the mass transfer coefficient
        // Here mass transfer coefficient is everything, of the mass
        // flux from one phase to another except the difference in mole fraction
        Scalar massTransferCoefficient[numPhases][numComponents];
        for(int smallLoopPhaseIdx=0; smallLoopPhaseIdx<numPhases; ++smallLoopPhaseIdx){
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                static_assert(numComponents == 2,
                              "this works only for binary mixtures, where D_AB = D_BA in the case of isothermal, isobaric conditions for diffusion");
                const Scalar rhoN   = volVars.fluidState().molarDensity(smallLoopPhaseIdx);
                const Scalar D = FluidSystem::binaryDiffusionCoefficient(volVars.fluidState(),
                                                        smallLoopPhaseIdx,
                                                        wCompIdx,
                                                        nCompIdx);
                Scalar Sh(0.)     ;
                Sh     = volVars.sherwoodNumber(smallLoopPhaseIdx) ;
                massTransferCoefficient[smallLoopPhaseIdx][compIdx] = volVars.interfacialArea(wPhaseIdx, nPhaseIdx)
                                                                       * volVars.fluidState().molarDensity(smallLoopPhaseIdx)
                                                                       * FluidSystem::binaryDiffusionCoefficient(volVars.fluidState(),
                                                                                                        smallLoopPhaseIdx,
                                                                                                        wCompIdx,
                                                                                                        nCompIdx)// see static assert above
                                                                        * volVars.sherwoodNumber(smallLoopPhaseIdx)
                                                                        / volVars.characteristicLength();

                Valgrind::CheckDefined(massTransferCoefficient[smallLoopPhaseIdx][compIdx]);
            }
        }


        Scalar xTransferNonEquil[numPhases][numComponents];
        {
            FluidState nonEquilFluidState;
            nonEquilFluidState.assign(volVars.fluidState()) ;
            NonEquilConstraintSolver::solve(nonEquilFluidState, // for storing the non-equilibrium mole fractions determined by mass transfer
                                            massTransferCoefficient,
                                            volVars,
                                            /*setViscosity=*/false,
                                            /*setEnthalpy=*/false) ;

            for(int smallLoopPhaseIdx=0; smallLoopPhaseIdx<numPhases; ++smallLoopPhaseIdx){
                for (int compIdx=0; compIdx< numComponents; ++ compIdx){
                    xTransferNonEquil[smallLoopPhaseIdx][compIdx] = nonEquilFluidState.moleFraction(smallLoopPhaseIdx, compIdx);
                }
            }
        }





        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                const Scalar deltaX             = xTransferNonEquil[phaseIdx][compIdx] - volVars.xEquil(phaseIdx, compIdx);
                const Scalar massTransferCoeff  = massTransferCoefficient[phaseIdx][compIdx];
                componentIntoPhaseMassTransfer[phaseIdx][compIdx] =  (-deltaX*massTransferCoeff );
            }
        }

#else // FUNKYMASSTRANSFER
        // diffusion coefficients in wetting phase
        const Scalar diffCoeffNinW = volVars.diffCoeff(wPhaseIdx, wCompIdx, nCompIdx) ;
        Valgrind::CheckDefined(diffCoeffNinW);
        // diffusion coefficients in non-wetting phase
        const Scalar diffCoeffWinN = volVars.diffCoeff(nPhaseIdx, wCompIdx, nCompIdx) ;
        Valgrind::CheckDefined(diffCoeffWinN);
        Scalar phaseDensity[numPhases];
        Scalar x[numPhases][numComponents]; // mass fractions in wetting phase
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
            phaseDensity[phaseIdx] = volVars.fluidState().molarDensity(phaseIdx);
            for (int compIdx=0; compIdx< numComponents; ++ compIdx){
                x[phaseIdx][compIdx] = volVars.fluidState().moleFraction(phaseIdx, compIdx);
            }
        }
        Valgrind::CheckDefined(x);

//  //      "equilibrium" values: calculated in volume variables
        const Scalar x_wPhaseNCompEquil = volVars.xEquil(wPhaseIdx, nCompIdx) ;   // very 2p2c
        const Scalar x_nPhaseWCompEquil = volVars.xEquil(nPhaseIdx, wCompIdx);    // very 2p2c
        Valgrind::CheckDefined(x_wPhaseNCompEquil);
        Valgrind::CheckDefined(x_nPhaseWCompEquil);
        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        const Scalar factorMassTransfer     = volVars.factorMassTransfer()  ;
        const Scalar awn = volVars.interfacialArea(wPhaseIdx, nPhaseIdx);
        const Scalar gradNinWApprox  =  (x[wPhaseIdx][nCompIdx] - x_wPhaseNCompEquil) / characteristicLength;    // very 2p2c
        const Scalar gradWinNApprox  =  (x[nPhaseIdx][wCompIdx] - x_nPhaseWCompEquil) / characteristicLength;    // very 2p2c

        const Scalar sherwoodWPhase  = volVars.sherwoodNumber(wPhaseIdx);
        const Scalar sherwoodNPhase  = volVars.sherwoodNumber(nPhaseIdx);

        //      actual diffusion is always calculated for eq's 2,3
        //      Eq's 1,4 have to be the same with different sign, because no mass is accumulated in the interface
        //      i.e. automatically conserving mass that mvoes across the interface

        const Scalar nCompIntoWPhase  = - factorMassTransfer * gradNinWApprox * awn * phaseDensity[wPhaseIdx] * diffCoeffNinW * sherwoodWPhase;
        const Scalar nCompIntoNPhase  = - nCompIntoWPhase ;
        const Scalar wCompIntoNPhase  = - factorMassTransfer * gradWinNApprox * awn * phaseDensity[nPhaseIdx] * diffCoeffWinN * sherwoodNPhase;
        const Scalar wCompIntoWPhase  = - wCompIntoNPhase ;

        componentIntoPhaseMassTransfer[wPhaseIdx][nCompIdx] = nCompIntoWPhase;
        componentIntoPhaseMassTransfer[nPhaseIdx][nCompIdx] = nCompIntoNPhase;
        componentIntoPhaseMassTransfer[nPhaseIdx][wCompIdx] = wCompIntoNPhase;
        componentIntoPhaseMassTransfer[wPhaseIdx][wCompIdx] = wCompIntoWPhase;

#endif

#if MASS_TRANSFER_OFF
#warning MASS TRANSFER OFF
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
            for (int compIdx=0; compIdx< numComponents; ++ compIdx){
                componentIntoPhaseMassTransfer[phaseIdx][compIdx] =  0.;
            }
        }
#endif

        Valgrind::CheckDefined(componentIntoPhaseMassTransfer);

        // Call the (kinetic) Energy module, for the source term.
        // it has to be called from here, because the mass transfered has to be known.
        EnergyResid::computeSource(source,
                                   volVars,
                                   componentIntoPhaseMassTransfer);

        // Actually add the mass transfer to the sources which might
        // exist externally
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                const unsigned int eqIdx = conti0EqIdx + compIdx + phaseIdx*numComponents;
                        source[eqIdx] += componentIntoPhaseMassTransfer[phaseIdx][compIdx] ;
                #ifndef NDEBUG
                if (!std::isfinite(source[eqIdx]))
                	DUNE_THROW(NumericalProblem, "Calculated non-finite source");
                #endif
            }
        }
        Valgrind::CheckDefined(source);
    }
};

} // end namespace Dumux

#endif
