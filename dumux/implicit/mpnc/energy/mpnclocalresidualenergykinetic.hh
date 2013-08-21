// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \brief This file contains the parts of the local residual to
 *        calculate the heat conservation in the thermal non-equilibrium M-phase
 *        N-component model.
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_KINETIC_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_KINETIC_HH

#include <dumux/implicit/mpnc/mpnclocalresidual.hh>

namespace Dumux {
template <class TypeTag>
class MPNCLocalResidualEnergy<TypeTag, /*enableEnergy=*/true, /*kineticEnergyTransfer=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim              = GridView::dimension };
    enum { dimWorld         = GridView::dimensionworld };
    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx      = Indices::conti0EqIdx };
    enum { energyEq0Idx     = Indices::energyEq0Idx };
    enum { numEnergyEqs     = Indices::NumPrimaryEnergyVars};
    enum { wPhaseIdx        = FluidSystem::wPhaseIdx};
    enum { nPhaseIdx        = FluidSystem::nPhaseIdx};
    enum { nCompIdx         = FluidSystem::nCompIdx};
    enum { wCompIdx         = FluidSystem::wCompIdx};
    enum { sPhaseIdx        = FluidSystem::sPhaseIdx};

    typedef Dune::FieldVector<Scalar, dim>                                  DimVector;
    typedef typename Dune::FieldVector<Scalar, numComponents>               ComponentVector;
    typedef typename Dune::FieldMatrix<Scalar, numPhases, numComponents>    PhaseComponentMatrix;

public:
    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub-control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param volVars the Volume Variables
     */
    static void computeStorage(PrimaryVariables & storage,
                                  const VolumeVariables & volVars)
    {
        for(int energyEqIdx=0; energyEqIdx< numEnergyEqs; energyEqIdx++)
            storage[energyEq0Idx+energyEqIdx] = 0;

        // energy of the fluids
        for (int phaseIdx = 0; phaseIdx < numEnergyEqs; ++phaseIdx) {
            addPhaseStorage(storage, volVars, phaseIdx);
        }
    }
    /*!
     * \brief BEWARE !!!
     * Here comes some problem with the phase-specific conservation of energy:
     * Imagine a displacement process, with the initial state being a fully saturated porous medium.
     * If the residing phase is now displaced by another phase, what happens at the front of the process?
     * At the very front there is one cell in which the invading phase enters and no invading (but residing) phase leaves.
     * With enthalpy in the flux term and internal energy in the storage term, the difference (pv) has to be converted into temperature in order to fulfill energy conservation.
     * -> A temperature peak at the front arises (if spatial discretization is sufficiently fine). This peak has a maximum value and does not increase with further refinement.
     * -> Further evidence for this explanation: in a simple setting (constant parameters, few cells) the temperature peak can be correctly predicted a priori.
     * -> -> For those situations with a distinct displacement process the same quantity has to be stored and transported
     * This is equivalent to neglecting volume changing work.
     *
     * \brief Calculate the storage for all mass balance equations
     *        within a single fluid phase
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param volVars the Volume Variables
     *  \param phaseIdx The local index of the phases
     */
    static void addPhaseStorage(PrimaryVariables & storage,
                                const VolumeVariables & volVars,
                                const unsigned int phaseIdx)
    {

        const FluidState & fs = volVars.fluidState();

        if (phaseIdx not_eq sPhaseIdx) {
            storage[energyEq0Idx + phaseIdx] +=
                    fs.density(phaseIdx) *
                    fs.internalEnergy(phaseIdx) *
                    fs.saturation(phaseIdx) *
                    volVars.porosity();
        }
        else if(phaseIdx == sPhaseIdx) {
            // heat stored in the rock matrix
            storage[energyEq0Idx+phaseIdx] += volVars.temperature(phaseIdx) *
                                               volVars.soilDensity() *
                                               (1.-volVars.porosity()) *
                                               volVars.heatCapacity();
        }
        else
            DUNE_THROW(Dune::NotImplemented,
            		"wrong index");
        #ifndef NDEBUG
        if (!std::isfinite(storage[energyEq0Idx+phaseIdx]))
        	DUNE_THROW(NumericalProblem, "Calculated non-finite storage");
        #endif

    }
    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param fluxVars The flux Variables
     * \param volVars The volume variables
     * \param molarPhaseComponentValuesMassTransport[numPhases]
     */
    static void computeFlux(PrimaryVariables & flux,
                            const FluxVariables & fluxVars,
                            const ElementVolumeVariables & elemVolVars,
                            const ComponentVector molarPhaseComponentValuesMassTransport[numPhases])
    {
        // reset all energy fluxes
        for(int energyEqIdx=0; energyEqIdx<numEnergyEqs; ++energyEqIdx)
            flux[energyEq0Idx + energyEqIdx] = 0.0;

        // only the fluid phases transport enthalpy
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
            computePhaseEnthalpyFlux(flux,
                                     fluxVars,
                                     elemVolVars,
                                     phaseIdx,
                                     molarPhaseComponentValuesMassTransport[phaseIdx]);

        // all phases take part in conduction
        for(int energyEqIdx=0; energyEqIdx<numEnergyEqs; ++energyEqIdx){
            computeHeatConduction(flux,
                                  fluxVars,
                                  elemVolVars,
                                  energyEqIdx);
            #ifndef NDEBUG
            if (!std::isfinite(flux[energyEq0Idx + energyEqIdx]))
            	DUNE_THROW(NumericalProblem, "Calculated non-finite flux");
            #endif
        }
    }
    /*!
          * \brief the advective Flux of the enthalpy
          *        \param flux The flux over the SCV (sub-control-volume) face for each component
          *        \param fluxVars The flux Variables
          *        \param elemVolVars The volume variables of the current element
          *        \param phaseIdx The local index of the phases
          */
    static void computePhaseEnthalpyFlux(PrimaryVariables & flux,
                                         const FluxVariables & fluxVars,
                                         const ElementVolumeVariables & elemVolVars,
                                         const unsigned int phaseIdx,
                                         const ComponentVector & molarComponentValuesMassTransport)
    {
        Scalar massFlux = 0; // mass flux is not the perfect term: sum_kappa (rho_alpha v_alpha x_alpha^kappa) = v_alpha rho_alpha

        // calculate the mass flux in the phase i.e. make mass flux out
        // of mole flux and add up the fluxes of a phase
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            massFlux += molarComponentValuesMassTransport[compIdx]
                          * FluidSystem::molarMass(compIdx)        ;

        unsigned int upIdx = fluxVars.face().i;
        if (massFlux < 0){
            upIdx = fluxVars.face().j;
        }

        // use the phase enthalpy of the upstream vertex to calculate
        // the enthalpy transport
        const VolumeVariables & up = elemVolVars[upIdx];

        /* todo
         * CAUTION: this is not exactly correct: does diffusion carry the upstream phase enthalpy? To be more precise this should be the components enthalpy. In the same vein: Counter current diffusion is not accounted for here.
         */
        const Scalar enthalpy =  up.fluidState().enthalpy(phaseIdx) ;

        flux[energyEq0Idx + phaseIdx] += enthalpy * massFlux  ;
    }
    /*!
        * \brief The heat conduction in the phase
        *
        *        \param flux The flux over the SCV (sub-control-volume) face for each component
        *        \param fluxVars The flux Variables
        *        \param elemVolVars The volume variables of the current element
        *        \param phaseIdx The local index of the phases
        */
    static void computeHeatConduction(PrimaryVariables & flux,
                                      const FluxVariables & fluxVars,
                                      const ElementVolumeVariables & elemVolVars,
                                      const unsigned int phaseIdx)
    {
        const unsigned int iIdx = fluxVars.face().i;
        const unsigned int kIdx = fluxVars.face().j; // k can be better distinguished from i

        const VolumeVariables & iVolVar = elemVolVars[iIdx];
        const VolumeVariables & kVolVar = elemVolVars[kIdx];

        const FluidState & iFluidState = iVolVar.fluidState();
        const FluidState & kFluidState = kVolVar.fluidState();

        const Scalar iPorosity      = iVolVar.porosity();
        const Scalar kPorosity      = kVolVar.porosity();
        const Scalar barPorosity    = Dumux::harmonicMean(iPorosity, kPorosity);

        const Scalar ilambda        = iVolVar.thermalConductivity(phaseIdx);
        const Scalar klambda        = kVolVar.thermalConductivity(phaseIdx);
        // todo which average?
        // Using a harmonic average is justified by its properties: if one phase does not conduct energy, there is no transfer
        const Scalar barLambda      = Dumux::harmonicMean(ilambda, klambda);

        const Scalar gradientNormal = fluxVars.fluxVarsEnergy().temperatureGradient(phaseIdx)
                                        * fluxVars.face().normal ;

        // heat conduction of the rock matrix and the fluid phases
        if (phaseIdx == sPhaseIdx){
            flux[energyEq0Idx + phaseIdx] -= barLambda * gradientNormal * (1.-barPorosity) ;
        }
        else if (phaseIdx == wPhaseIdx or phaseIdx == nPhaseIdx){
            const Scalar iSaturation    = iFluidState.saturation(phaseIdx);
            const Scalar kSaturation    = kFluidState.saturation(phaseIdx);
            const Scalar barSaturation = Dumux::harmonicMean(iSaturation, kSaturation);
            flux[energyEq0Idx + phaseIdx] -= barLambda * gradientNormal *  barPorosity * barSaturation  ;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                               "wrong index");
    }
    /*!
       * \brief Calculate the source term of the equation
       *
       * \param source The source/sink in the sub-control volume for each component
       * \param volVars The volume variables
       * \param componentIntoPhaseMassTransfer[numPhases]
       */

    static void computeSource(PrimaryVariables & source,
                              const VolumeVariables & volVars,
                              const ComponentVector componentIntoPhaseMassTransfer[numPhases])
    {
        const Scalar awn = volVars.interfacialArea(wPhaseIdx, nPhaseIdx);
        const Scalar aws = volVars.interfacialArea(wPhaseIdx, sPhaseIdx);
        const Scalar ans = volVars.interfacialArea(nPhaseIdx, sPhaseIdx);

        const Scalar Tw = volVars.temperature(wPhaseIdx);
        const Scalar Tn = volVars.temperature(nPhaseIdx);
        const Scalar Ts = volVars.temperature(sPhaseIdx);

        const  Scalar lambdaWetting     = volVars.thermalConductivity(wPhaseIdx);
        const  Scalar lambdaNonWetting  = volVars.thermalConductivity(nPhaseIdx);
        const  Scalar lambdaSolid       = volVars.thermalConductivity(sPhaseIdx);

        // todo which average?
        // Using a harmonic average is justified by its properties: if one phase does not conduct energy, there is no transfer
        const Scalar lambdaWN      = Dumux::harmonicMean(lambdaWetting, lambdaNonWetting);
        const Scalar lambdaWS      = Dumux::harmonicMean(lambdaWetting, lambdaSolid);
        const Scalar lambdaNS      = Dumux::harmonicMean(lambdaNonWetting, lambdaSolid);
//      |------------------------------------------------------|
//      |                          |                           |
//      |                          |                           |
//      |        alpha             |            i              |
//      |       T_alpha            |           T_i             |
//      |                          |                           |
//      |                          |                           |
//      |------------------------------------------------------|

//      T_i > T_\alpha i.e. heat going into \alpha phase

//        Q_{i \leadsto \alpha} = a_{\alpha i} lambda_{\alpha i} (T_i - T_\alpha) / d // i.e.: this is the r.h.s. of alpha
        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        const Scalar factorEnergyTransfer   = volVars.factorEnergyTransfer()  ;

        const Scalar nusseltWN      = Dumux::harmonicMean(volVars.nusseltNumber(wPhaseIdx), volVars.nusseltNumber(nPhaseIdx));
        const Scalar nusseltWS      = volVars.nusseltNumber(wPhaseIdx);
        const Scalar nusseltNS      = volVars.nusseltNumber(nPhaseIdx);

//#warning SET NUSSELT TO 1
//        const Scalar nusseltWN      = 1. ;
//        const Scalar nusseltWS      = 1. ;
//        const Scalar nusseltNS      = 1. ;

        const Scalar wettingToNonWettingEnergyExchange = factorEnergyTransfer * (Tw - Tn) / characteristicLength * awn * lambdaWN * nusseltWN  ;
        const Scalar wettingToSolidEnergyExchange      = factorEnergyTransfer * (Tw - Ts) / characteristicLength * aws * lambdaWS * nusseltWS  ;
        const Scalar nonWettingToSolidEnergyExchange   = factorEnergyTransfer * (Tn - Ts) / characteristicLength * ans * lambdaNS * nusseltNS  ;

//#warning HEAT TRANSFER OFF
//        const Scalar WettingToNonWettingEnergyExchange = 0. ;
//        const Scalar WettingToSolidEnergyExchange      = 0. ;
//        const Scalar NonWettingToSolidEnergyExchange   = 0. ;


        for(int phaseIdx =0; phaseIdx<numEnergyEqs; ++phaseIdx){
            switch (phaseIdx){
            case wPhaseIdx:
                source[energyEq0Idx + phaseIdx] =  ( - wettingToNonWettingEnergyExchange - wettingToSolidEnergyExchange);
                break;
            case nPhaseIdx:
                source[energyEq0Idx + phaseIdx] =  (+ wettingToNonWettingEnergyExchange - nonWettingToSolidEnergyExchange);
                break;
            case sPhaseIdx:
                source[energyEq0Idx + phaseIdx] =  (+ wettingToSolidEnergyExchange + nonWettingToSolidEnergyExchange);
                break;
            default:
                DUNE_THROW(Dune::NotImplemented,
                           "wrong index");
            } // end switch



            #ifndef NDEBUG
            if (!std::isfinite(source[energyEq0Idx + phaseIdx]))
            	DUNE_THROW(NumericalProblem, "Calculated non-finite source, " << "Tw="<< Tw << " Tn="<< Tn<< " Ts="<< Ts);
            #endif
        }// end phases

#define MASS_ENERGY_TRANSPORT 1
#if MASS_ENERGY_TRANSPORT
// Here comes the catch: We are not doing energy conservation for the whole
// system, but rather for each individual phase.
//        -> Therefore the energy fluxes over each phase boundary need be
// individually accounted for.
//        -> Each particle crossing a phase boundary does carry some mass and
//        thus energy!
//        -> Therefore, this contribution needs to be added.

//        -> the particle always brings the energy of the originating phase.
//        -> Energy advectivly transported into a phase = the moles of a component that go into a phase * molMass * enthalpy of the component in the *originating* phase

        // The fluidsystem likes to get a fluidstate. ...
        const FluidState & fluidState = volVars.fluidState();

        for(int phaseIdx =0; phaseIdx<numEnergyEqs; ++phaseIdx){
            switch (phaseIdx){
            case wPhaseIdx:
                source[energyEq0Idx + phaseIdx] += (componentIntoPhaseMassTransfer[wPhaseIdx][nCompIdx]
                                                                                              * FluidSystem::molarMass(nCompIdx)
                                                                                              * FluidSystem::componentEnthalpy(fluidState, nPhaseIdx, nCompIdx) );
                source[energyEq0Idx + phaseIdx] += (componentIntoPhaseMassTransfer[wPhaseIdx][wCompIdx]
                                                                                              * FluidSystem::molarMass(wCompIdx)
                                                                                              * FluidSystem::componentEnthalpy(fluidState, nPhaseIdx, wCompIdx));
                break;
            case nPhaseIdx:
                source[energyEq0Idx + phaseIdx] += (componentIntoPhaseMassTransfer[nPhaseIdx][nCompIdx]
                                                                                              * FluidSystem::molarMass(nCompIdx)
                                                                                              * FluidSystem::componentEnthalpy(fluidState, wPhaseIdx, nCompIdx));
                source[energyEq0Idx + phaseIdx] += (componentIntoPhaseMassTransfer[nPhaseIdx][wCompIdx]
                                                                                              * FluidSystem::molarMass(wCompIdx)
                                                                                              * FluidSystem::componentEnthalpy(fluidState, wPhaseIdx, wCompIdx));
                break;
            case sPhaseIdx:
                break; // no sorption
            default:
                DUNE_THROW(Dune::NotImplemented,
                           "wrong index");
            } // end switch
        }// end phases
#endif //MASS_ENERGY_TRANSPORT
        Valgrind::CheckDefined(source);
    }// end source
};

};

#endif // DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_KINETIC_HH
