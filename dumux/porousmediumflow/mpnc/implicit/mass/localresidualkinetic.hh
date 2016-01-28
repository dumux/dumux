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
 * \brief The local residual for the kinetic mass transfer module of
 *        the compositional multi-phase model.
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_MASS_KINETIC_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_MASS_KINETIC_HH

#include <dumux/porousmediumflow/mpnc/implicit/mass/localresidual.hh>

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

    enum { numPhases =  GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { nCompIdx = FluidSystem::nCompIdx } ;
    enum { wCompIdx = FluidSystem::wCompIdx } ;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx} ;
    enum { nPhaseIdx = FluidSystem::nPhaseIdx} ;
    enum { sPhaseIdx = FluidSystem::sPhaseIdx} ;
    enum { numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef MPNCLocalResidualEnergy<TypeTag, enableEnergy, numEnergyEquations> EnergyResid;
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

            if (!std::isfinite(storage[conti0EqIdx + phaseIdx*numComponents + compIdx]))
                DUNE_THROW(NumericalProblem, "Calculated non-finite storage");
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

                if (!std::isfinite(flux[conti0EqIdx + phaseIdx*numComponents + compIdx]))
                    DUNE_THROW(NumericalProblem, "Calculated non-finite flux");
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
#define FUNKYMASSTRANSFER 0
#if FUNKYMASSTRANSFER
        const Scalar mu_nPhaseNCompEquil = volVars.chemicalPotentialEquil(nPhaseIdx, nCompIdx) ;   // very 2p2c
        const Scalar mu_wPhaseWCompEquil = volVars.chemicalPotentialEquil(wPhaseIdx, wCompIdx);    // very 2p2c

        const Scalar mu_wPhaseNComp = volVars.chemicalPotential(wPhaseIdx, nCompIdx) ;   // very 2p2c
        const Scalar mu_nPhaseWComp = volVars.chemicalPotential(nPhaseIdx, wCompIdx);    // very 2p2c

        Valgrind::CheckDefined(mu_nPhaseNCompEquil);
        Valgrind::CheckDefined(mu_wPhaseWCompEquil);
        Valgrind::CheckDefined(mu_wPhaseNComp);
        Valgrind::CheckDefined(mu_nPhaseWComp);

        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        const Scalar temperature            = volVars.temperature(wPhaseIdx);
        const Scalar pn                     = volVars.pressure(nPhaseIdx);
        const Scalar henry                  = FluidSystem::henry(temperature) ;
        const Scalar gradNinWApprox  = ( mu_wPhaseNComp - mu_nPhaseNCompEquil) / characteristicLength;    // very 2p2c // 1. / henry *
        const Scalar gradWinNApprox  = ( mu_nPhaseWComp - mu_wPhaseWCompEquil) / characteristicLength;    // very 2p2c // 1. / pn *

#else // FUNKYMASSTRANSFER
        Scalar x[numPhases][numComponents]; // mass fractions in wetting phase
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
            for (int compIdx=0; compIdx< numComponents; ++ compIdx){
                x[phaseIdx][compIdx] = volVars.moleFraction(phaseIdx, compIdx);
            }
        }
        Valgrind::CheckDefined(x);

//      "equilibrium" values: calculated in volume variables
        const Scalar x_wPhaseNCompEquil = volVars.xEquil(wPhaseIdx, nCompIdx) ;   // very 2p2c
        const Scalar x_nPhaseWCompEquil = volVars.xEquil(nPhaseIdx, wCompIdx);    // very 2p2c
        Valgrind::CheckDefined(x_wPhaseNCompEquil);
        Valgrind::CheckDefined(x_nPhaseWCompEquil);
        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        const Scalar gradNinWApprox  =  (x[wPhaseIdx][nCompIdx] - x_wPhaseNCompEquil) / characteristicLength;    // very 2p2c
        const Scalar gradWinNApprox  =  (x[nPhaseIdx][wCompIdx] - x_nPhaseWCompEquil) / characteristicLength;    // very 2p2c
#endif
        Scalar phaseDensity[numPhases];
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
            phaseDensity[phaseIdx] = volVars.molarDensity(phaseIdx);
        }

        // diffusion coefficients in wetting phase
        const Scalar diffCoeffNinW = volVars.diffCoeff(wPhaseIdx, wCompIdx, nCompIdx) ;
        Valgrind::CheckDefined(diffCoeffNinW);
        // diffusion coefficients in non-wetting phase
        const Scalar diffCoeffWinN = volVars.diffCoeff(nPhaseIdx, wCompIdx, nCompIdx) ;
        Valgrind::CheckDefined(diffCoeffWinN);

        const Scalar factorMassTransfer     = volVars.factorMassTransfer()  ;
        const Scalar awn = volVars.interfacialArea(wPhaseIdx, nPhaseIdx);

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

                        if (!std::isfinite(source[eqIdx]))
                            DUNE_THROW(NumericalProblem, "Calculated non-finite source");
            }
        }
        Valgrind::CheckDefined(source);
    }
};

} // end namespace Dumux

#endif
