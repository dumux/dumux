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

#include <dumux/porousmediumflow/mpnc/implicit/localresidual.hh>
#include <dumux/common/spline.hh>


namespace Dumux
{

/*!
 * \brief Specialization for the case of *3* energy balance equations.
 */
template <class TypeTag>
class MPNCLocalResidualEnergy<TypeTag, /*enableEnergy=*/true, /*numEnergyEquations=*/3>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { energyEq0Idx     = Indices::energyEq0Idx };
    enum { numEnergyEqs     = Indices::numPrimaryEnergyVars};
    enum { wPhaseIdx        = FluidSystem::wPhaseIdx};
    enum { nPhaseIdx        = FluidSystem::nPhaseIdx};
    enum { nCompIdx         = FluidSystem::nCompIdx};
    enum { wCompIdx         = FluidSystem::wCompIdx};
    enum { sPhaseIdx        = FluidSystem::sPhaseIdx};

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
     * With enthalpy in the flux term and internal energy in the storage term,
     * the difference (pv) has to be converted into temperature in order to fulfill energy conservation.
     * -> A temperature peak at the front arises (if spatial discretization is sufficiently fine).
     * This peak has a maximum value and does not increase with further refinement.
     * -> Further evidence for this explanation: in a simple setting (constant parameters,
     * few cells) the temperature peak can be correctly predicted a priori.
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
                                               volVars.solidDensity() *
                                               (1.-volVars.porosity()) *
                                               volVars.solidHeatCapacity();
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                    "wrong index");

        if (!std::isfinite(storage[energyEq0Idx+phaseIdx]))
            DUNE_THROW(NumericalProblem, "Calculated non-finite storage");

    }
    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param fluxVars The flux Variables
     * \param elemVolVars The volume variables of the current element
     * \param molarPhaseComponentValuesMassTransport
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
            if (!std::isfinite(flux[energyEq0Idx + energyEqIdx]))
                DUNE_THROW(NumericalProblem, "Calculated non-finite flux in phase " << energyEqIdx);
        }
    }
    /*!
          * \brief the advective Flux of the enthalpy
          *        \param flux The flux over the SCV (sub-control-volume) face for each component
          *        \param fluxVars The flux Variables
          *        \param elemVolVars The volume variables of the current element
          *        \param phaseIdx The local index of the phases
          *        \param molarComponentValuesMassTransport
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
        for (int compIdx = 0; compIdx < numComponents; ++compIdx){
            massFlux += molarComponentValuesMassTransport[compIdx]
                          * FluidSystem::molarMass(compIdx)        ;
        }

        unsigned int upIdx = fluxVars.face().i;
        unsigned int dnIdx = fluxVars.face().j;
        if (massFlux < 0){
            upIdx = fluxVars.face().j;
            dnIdx = fluxVars.face().i;
        }

        // use the phase enthalpy of the upstream vertex to calculate
        // the enthalpy transport
        const VolumeVariables & up = elemVolVars[upIdx];
        const VolumeVariables & dn = elemVolVars[dnIdx];


        // CAUTION: this is not exactly correct: does diffusion carry the upstream phase enthalpy?
        // To be more precise this should be the components enthalpy.
        // In the same vein: Counter current diffusion is not accounted for here.
        const Scalar transportedThingUp =  up.enthalpy(phaseIdx) ;
        const Scalar transportedThingDn =  dn.enthalpy(phaseIdx) ;


        const Scalar massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
        flux[energyEq0Idx + phaseIdx] += massFlux *
                                            (massUpwindWeight_ * transportedThingUp
                                                    +
                                            (1.-massUpwindWeight_) * transportedThingDn ) ;
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
        const Scalar barPorosity    = harmonicMean(iPorosity, kPorosity);

        const Scalar ilambda        = iVolVar.thermalConductivity(phaseIdx);
        const Scalar klambda        = kVolVar.thermalConductivity(phaseIdx);

        // Using a harmonic average is justified by its properties: if one phase does not conduct energy, there is no transfer
        const Scalar barLambda      = harmonicMean(ilambda, klambda) ;

        const Scalar gradientNormal = fluxVars.fluxVarsEnergy().temperatureGradient(phaseIdx)
                                        * fluxVars.face().normal ;

        // heat conduction of the rock matrix and the fluid phases
        if (phaseIdx == sPhaseIdx){
            flux[energyEq0Idx + phaseIdx] -= barLambda * gradientNormal * (1.-barPorosity) ;
        }
        else if (phaseIdx == wPhaseIdx or phaseIdx == nPhaseIdx){
            const Scalar iSaturation    = iFluidState.saturation(phaseIdx);
            const Scalar kSaturation    = kFluidState.saturation(phaseIdx);
            const Scalar barSaturation = harmonicMean(iSaturation, kSaturation);
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
       * \param componentIntoPhaseMassTransfer
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

        // Using a harmonic average is justified by its properties: if one phase does not conduct energy, there is no transfer
        const Scalar lambdaWN      = harmonicMean(lambdaWetting, lambdaNonWetting);
        const Scalar lambdaWS      = harmonicMean(lambdaWetting, lambdaSolid);
        const Scalar lambdaNS      = harmonicMean(lambdaNonWetting, lambdaSolid);
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

        const Scalar nusseltWN      = harmonicMean(volVars.nusseltNumber(wPhaseIdx), volVars.nusseltNumber(nPhaseIdx));
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
//        const Scalar wettingToNonWettingEnergyExchange = 0. ;
//        const Scalar wettingToSolidEnergyExchange      = 0. ;
//        const Scalar nonWettingToSolidEnergyExchange   = 0. ;


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



            if (!std::isfinite(source[energyEq0Idx + phaseIdx]))
                DUNE_THROW(NumericalProblem, "Calculated non-finite source, " << "Tw="<< Tw << " Tn="<< Tn<< " Ts="<< Ts);
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
//        -> Energy advectivly transported into a phase = the moles of a component that go into a
//           phase * molMass * enthalpy of the component in the *originating* phase

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


/*!
 * \brief Specialization for the case of *2* energy balance equations.
 */
template <class TypeTag>
class MPNCLocalResidualEnergy<TypeTag, /*enableEnergy=*/true, /*numEnergyEquations=*/2>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { energyEq0Idx     = Indices::energyEq0Idx };
    enum { numEnergyEqs     = Indices::numPrimaryEnergyVars};
    enum { wPhaseIdx        = FluidSystem::wPhaseIdx};
    enum { nPhaseIdx        = FluidSystem::nPhaseIdx};
    enum { nCompIdx         = FluidSystem::nCompIdx};
    enum { wCompIdx         = FluidSystem::wCompIdx};
    enum { sPhaseIdx        = FluidSystem::sPhaseIdx};
    enum { energyEqSolidIdx = Indices::energyEqSolidIdx};
    enum { temperatureFluidIdx = Indices::temperatureFluidIdx};
    enum { temperatureSolidIdx = Indices::temperatureSolidIdx};
    enum { dim = GridView::dimension}; // Grid and world dimension


    typedef typename Dune::FieldVector<Scalar, numComponents>               ComponentVector;
    typedef typename Dune::FieldMatrix<Scalar, numPhases, numComponents>    PhaseComponentMatrix;
    typedef typename FluidSystem::ParameterCache ParameterCache;

    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef Dumux::Spline<Scalar> Spline;


public:
    /*! \brief Evaluate the amount all conservation quantities
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
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            addPhaseStorage(storage, volVars, phaseIdx);
        }
        // heat stored in the rock matrix
        storage[energyEqSolidIdx] +=
            volVars.temperature(temperatureSolidIdx)
            * volVars.solidDensity()
            * (1.0 - volVars.porosity())
            * volVars.solidHeatCapacity();

        if (!std::isfinite(storage[energyEqSolidIdx]))
            DUNE_THROW(NumericalProblem, "Calculated non-finite storage");
    }

    /*! \brief Calculate the storage for all mass balance equations
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
            storage[energyEq0Idx ] +=
                    fs.density(phaseIdx) *
                    fs.internalEnergy(phaseIdx) *
                    fs.saturation(phaseIdx) *
                    volVars.porosity();
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                    "wrong index");

        if (!std::isfinite(storage[energyEq0Idx]))
            DUNE_THROW(NumericalProblem, "Calculated non-finite storage");
    }

    /*!\brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param fluxVars The flux Variables
     * \param elemVolVars The volume variables of the current element
     * \param molarPhaseComponentValuesMassTransport
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

            if (!std::isfinite(flux[energyEq0Idx + energyEqIdx]))
                DUNE_THROW(NumericalProblem, "Calculated non-finite flux in phase " << energyEqIdx);
        }
    }

    /*! \brief the advective Flux of the enthalpy
      *        \param flux The flux over the SCV (sub-control-volume) face for each component
      *        \param fluxVars The flux Variables
      *        \param elemVolVars The volume variables of the current element
      *        \param phaseIdx The local index of the phases
      *        \param molarComponentValuesMassTransport
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
        for (int compIdx = 0; compIdx < numComponents; ++compIdx){
            massFlux += molarComponentValuesMassTransport[compIdx]
                          * FluidSystem::molarMass(compIdx)        ;
        }

        unsigned int upIdx = fluxVars.face().i;
        unsigned int dnIdx = fluxVars.face().j;
        if (massFlux < 0){
            upIdx = fluxVars.face().j;
            dnIdx = fluxVars.face().i;
        }

        // use the phase enthalpy of the upstream vertex to calculate
        // the enthalpy transport
        const VolumeVariables & up = elemVolVars[upIdx];
        const VolumeVariables & dn = elemVolVars[dnIdx];

        // CAUTION: this is not exactly correct: does diffusion carry the upstream phase enthalpy?
        // To be more precise this should be the components enthalpy.
        // In the same vein: Counter current diffusion is not accounted for here.
        const Scalar transportedThingUp =  up.enthalpy(phaseIdx) ;
        const Scalar transportedThingDn =  dn.enthalpy(phaseIdx) ;

        const Scalar massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
        flux[energyEq0Idx] += massFlux *
                                            (massUpwindWeight_ * transportedThingUp
                                                    +
                                            (1.-massUpwindWeight_) * transportedThingDn ) ;
    }

    /*! \brief The heat conduction in the phase
     *
     *  \param flux The flux over the SCV (sub-control-volume) face for each component
     *  \param fluxVars The flux Variables
     *  \param elemVolVars The volume variables of the current element
     *  \param energyEqIdx The index of the phase energy equation
     */
    static void computeHeatConduction(PrimaryVariables & flux,
                                      const FluxVariables & fluxVars,
                                      const ElementVolumeVariables & elemVolVars,
                                      const unsigned int energyEqIdx)
    {
        const unsigned int iIdx = fluxVars.face().i;
        const unsigned int kIdx = fluxVars.face().j; // k can be better distinguished from i

        const VolumeVariables & iVolVar = elemVolVars[iIdx];
        const VolumeVariables & kVolVar = elemVolVars[kIdx];

        const Scalar iPorosity      = iVolVar.porosity();
        const Scalar kPorosity      = kVolVar.porosity();
        const Scalar barPorosity    = harmonicMean(iPorosity, kPorosity);

        const Scalar iSolidLambda        = iVolVar.thermalConductivity(sPhaseIdx);
        const Scalar kSolidLambda        = kVolVar.thermalConductivity(sPhaseIdx);

        // Using a harmonic average is justified by its properties: if one phase does not conduct energy, there is no transfer
        const Scalar barSolidLambda      = (iSolidLambda+kSolidLambda) / 2.0;

        const Scalar lumpedConductivity   = fluxVars.fluxVarsEnergy().lambdaEff() ;

        const Scalar gradientNormal = fluxVars.fluxVarsEnergy().temperatureGradient(energyEqIdx)
                                        * fluxVars.face().normal ;

        // heat conduction of the rock matrix and the fluid phases
        if (energyEqIdx == temperatureSolidIdx){
            flux[energyEqSolidIdx] -= barSolidLambda * (1.-barPorosity) * gradientNormal  ;
        }
        else if (energyEqIdx == temperatureFluidIdx){
            flux[energyEq0Idx ] -= lumpedConductivity /*already includes porosity*/ * gradientNormal  ;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                               "wrong index");
    }

    /*! \brief Calculate the source term of the equation
   *
   * \param source The source/sink in the sub-control volume for each component
   * \param volVars The volume variables
   * \param componentIntoPhaseMassTransfer
   */
    static void computeSource(PrimaryVariables & source,
                              const VolumeVariables & volVars,
                              const ComponentVector componentIntoPhaseMassTransfer[numPhases])
    {
        const Scalar solidToFluidEnergyExchange = qsf(volVars) ;

        for(int energyEqIdx =0; energyEqIdx<numEnergyEqs; ++energyEqIdx){
            switch (energyEqIdx){
            case 0 :
                source[energyEq0Idx + energyEqIdx] =  solidToFluidEnergyExchange;
                break;
            case 1 :
                source[energyEq0Idx + energyEqIdx] =  - solidToFluidEnergyExchange;
                break;
            default:
                DUNE_THROW(Dune::NotImplemented,
                           "wrong index");
            } // end switch
        }// end energyEqIdx
        Valgrind::CheckDefined(source);
    }// end source


    /*! \brief Calculate the whole energy transfer
   *
   * \param volVars The volume variables
   */
    static Scalar qsf(const VolumeVariables & volVars)
    {
        const FluidState & fs = volVars.fluidState() ;
        const Scalar characteristicLength   = volVars.characteristicLength()  ;

        // Shi & Wang, Transport in porous media (2011)
        const Scalar as = 6.0 * (1.0-volVars.porosity()) / characteristicLength ;

        const Scalar TFluid     = volVars.temperature(temperatureFluidIdx);
        const Scalar TSolid     = volVars.temperature(temperatureSolidIdx);

        const Scalar satW       = fs.saturation(wPhaseIdx) ;
        const Scalar satN       = fs.saturation(nPhaseIdx) ;

        const Scalar eps = 1e-6 ;
        Scalar solidToFluidEnergyExchange ;

        Scalar fluidConductivity ;
        if (satW < 1.0 - eps)
            fluidConductivity = volVars.thermalConductivity(nPhaseIdx) ;
        else if (satW >= 1.0 - eps)
            fluidConductivity = volVars.thermalConductivity(wPhaseIdx) ;
        else
            DUNE_THROW(Dune::NotImplemented,
                       "wrong range");

        const Scalar factorEnergyTransfer   = volVars.factorEnergyTransfer()  ;

        solidToFluidEnergyExchange = factorEnergyTransfer * (TSolid - TFluid) / characteristicLength * as * fluidConductivity ;

        const Scalar epsRegul = 1e-3 ;

        if (satW < (0 + eps) )
            solidToFluidEnergyExchange *=  volVars.nusseltNumber(nPhaseIdx) ;

        else if ( (satW >= 0 + eps) and (satW < 1.0-eps) ){
            solidToFluidEnergyExchange *=  (volVars.nusseltNumber(nPhaseIdx) * satN );
            Scalar qBoil ;

            if (satW<=epsRegul){// regularize
                Spline sp(0.0,                      epsRegul,                           // x1, x2
                          QBoilFunc(volVars, 0.0),  QBoilFunc(volVars, epsRegul),       // y1, y2
                          0.0,                      dQBoil_dSw(volVars, epsRegul) );    // m1, m2

                qBoil = sp.eval(satW) ;
            }

            else if (satW>= (1.0-epsRegul) ){// regularize
                Spline sp(1.0-epsRegul,                                     1.0,    // x1, x2
                          QBoilFunc(volVars, 1.0-epsRegul),                 0.0,    // y1, y2
                          dQBoil_dSw(volVars, 1.0-epsRegul),                    0.0 );      // m1, m2

                qBoil = sp.eval(satW) ;
            }
            else
                qBoil = QBoilFunc(volVars, satW)  ;

            solidToFluidEnergyExchange += qBoil;
        }
        else if (satW >= 1.0-eps)
            solidToFluidEnergyExchange *=  volVars.nusseltNumber(wPhaseIdx) ;
        else
            DUNE_THROW(Dune::NotImplemented,
                       "wrong range");

        if (!std::isfinite(solidToFluidEnergyExchange))
                        DUNE_THROW(NumericalProblem, "Calculated non-finite source, " << "TFluid="<< TFluid << " TSolid="<< TSolid  );

        return solidToFluidEnergyExchange ;
    }

    /*! \brief Calculate the energy transfer during boiling, i.e. latent heat
   *
   * \param volVars The volume variables
   * \param satW The wetting phase saturation. Not taken from volVars, because we regularize.
   */
    static Scalar QBoilFunc(const VolumeVariables & volVars,
                            const  Scalar satW)
    {
        // using saturation as input (instead of from volVars)
        // in order to make regularization (evaluation at different points) easyer
        const FluidState & fs = volVars.fluidState() ;
    const Scalar g( 9.81 ) ;
    const Scalar gamma(0.0589) ;
        const Scalar TSolid     = volVars.temperature(temperatureSolidIdx);
        const Scalar characteristicLength   = volVars.characteristicLength()  ;

    const Scalar as = 6.0 * (1.0-volVars.porosity()) / characteristicLength ;
    const Scalar mul = fs.viscosity(wPhaseIdx) ;
        const Scalar deltahv = fs.enthalpy(nPhaseIdx) - fs.enthalpy(wPhaseIdx);
        const Scalar deltaRho = fs.density(wPhaseIdx) - fs.density(nPhaseIdx) ;
    const Scalar firstBracket = std::pow(g * deltaRho / gamma, 0.5);
    const Scalar cp = FluidSystem::heatCapacity(fs, wPhaseIdx) ;
    // This use of Tsat is only justified if the fluid is always boiling (tsat equals boiling conditions)
    // If a different state is to be simulated, please use the actual fluid temperature instead.
    const Scalar Tsat = FluidSystem::vaporTemperature(fs, nPhaseIdx ) ;
    const Scalar deltaT = TSolid - Tsat ;
    const Scalar secondBracket = std::pow( (cp *deltaT / (0.006 * deltahv)  ) , 3.0 ) ;
    const Scalar Prl = volVars.prandtlNumber(wPhaseIdx) ;
    const Scalar thirdBracket = std::pow( 1/Prl , (1.7/0.33) );
    const Scalar QBoil = satW * as * mul * deltahv * firstBracket * secondBracket * thirdBracket ;
        return QBoil;
    }

    /*! \brief Calculate the derivative of the energy transfer function during boiling. Needed for regularization.
   *
   * \param volVars The volume variables
   * \param satW The wetting phase saturation. Not taken from volVars, because we regularize.
   */
    static Scalar dQBoil_dSw(const VolumeVariables & volVars,
                                const Scalar satW)
    {
        // on the fly derive w.r.t. Sw.
        // Only linearly depending on it (directly)
        return (QBoilFunc(volVars, satW) / satW ) ;
    }
};


} // end namespace Dumux

#endif // DUMUX_MPNC_LOCAL_RESIDUAL_ENERGY_KINETIC_HH
