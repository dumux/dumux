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
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief The local residual for the kinetic mass transfer module of
 *        the compositional multi-phase model.
 */
#ifndef DUMUX_NONEQUILIBRIUM_LOCAL_RESIDUAL_HH
#define DUMUX_NONEQUILIBRIUM_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/nonequilibrium/thermal/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, bool enableThermalNonEquilibrium, bool enableChemicalNonEquilibrium>
class NonEquilibriumLocalResidualImplementation;

template <class TypeTag>
using NonEquilibriumLocalResidual = NonEquilibriumLocalResidualImplementation<TypeTag, GET_PROP_TYPE(TypeTag, ModelTraits)::enableThermalNonEquilibrium(), GET_PROP_TYPE(TypeTag, ModelTraits)::enableChemicalNonEquilibrium()>;

/*!
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief The mass conservation part of the nonequilibrium model for a model without chemical non-equilibrium
 */
template<class TypeTag>
class NonEquilibriumLocalResidualImplementation<TypeTag, true, false>: public GET_PROP_TYPE(TypeTag, EquilibriumLocalResidual)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ParentType = typename GET_PROP_TYPE(TypeTag, EquilibriumLocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numPhases();
    static constexpr int numComponents = ModelTraits::numComponents();
    enum { conti0EqIdx = Indices::conti0EqIdx };
public:
    using ParentType::ParentType;

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param scv The sub-control volume over which we integrate the source term
     * \note This is the default implementation for all models as sources are computed
     *       in the user interface of the problem
     *
     */

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        // get upwind weights into local scope
        NumEqVector flux(0.0);

        const auto moleDensity = [](const auto& volVars, const int phaseIdx)
        { return volVars.molarDensity(phaseIdx); };

        const auto moleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return  volVars.moleFraction(phaseIdx, compIdx); };

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = conti0EqIdx + compIdx;

                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&moleDensity, &moleFraction, phaseIdx, compIdx] (const auto& volVars)
                { return moleDensity(volVars, phaseIdx)*moleFraction(volVars, phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                    flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                // diffusive fluxes (only for the component balances)
                    flux[eqIdx] += diffusiveFluxes[compIdx];
            }

            //! Add advective and diffusive phase energy fluxes
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, elemVolVars,scvf, phaseIdx);
        }
        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;


    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);
        // Call the (kinetic) Energy module, for the source term.
        // it has to be called from here, because the mass transfered has to be known.
        EnergyLocalResidual::computeSourceEnergy(source,
                                                 element,
                                                 fvGeometry,
                                                 elemVolVars,
                                                 scv);
        return source;
    }

};


/*!
 * \brief The mass conservation part of the nonequilibrium model for a model assuming chemical non-equilibrium and two phases */
template<class TypeTag>
class NonEquilibriumLocalResidualImplementation<TypeTag, true, true>: public GET_PROP_TYPE(TypeTag, EquilibriumLocalResidual)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, EquilibriumLocalResidual);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numPhases();
    static constexpr int numComponents = ModelTraits::numComponents();
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { comp1Idx = FluidSystem::comp1Idx } ;
    enum { comp0Idx = FluidSystem::comp0Idx } ;
    enum { phase0Idx = FluidSystem::phase0Idx} ;
    enum { phase1Idx = FluidSystem::phase1Idx} ;

public:
     using ParentType::ParentType;
    /*!
     * \brief Calculate the storage for all mass balance equations
     *
     *    \param storage The mass of the component within the sub-control volume
     *    \param volVars The volume variables
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
       NumEqVector storage(0.0);

     // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                 auto eqIdx = conti0EqIdx + phaseIdx*numComponents + compIdx;
                 storage[eqIdx] += volVars.porosity()
                                * volVars.saturation(phaseIdx)
                                * volVars.molarDensity(phaseIdx)
                                * volVars.moleFraction(phaseIdx, compIdx);
            }
            //! The energy storage in the fluid phase with index phaseIdx
            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }
         //! The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);
        return storage;
    }

    /*!
     * \brief Calculate the storage for all mass balance equations
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param fluxVars The flux Variables
     *        \param elemVolVars The volume variables of the current element
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        // get upwind weights into local scope
        NumEqVector flux(0.0);

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx =  conti0EqIdx + phaseIdx*numComponents + compIdx;
                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [phaseIdx, compIdx] (const auto& volVars)
                { return volVars.molarDensity(phaseIdx)*volVars.moleFraction(phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                    flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                // diffusive fluxes (only for the component balances)
                if (compIdx == phaseIdx) //do not add diffusive flux of main component, as that is not done in master as well
                        continue;
                    flux[eqIdx] += diffusiveFluxes[compIdx];
            }
            //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, elemVolVars,scvf, phaseIdx);
        }

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;
    }


    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
         NumEqVector source(0.0);
        // In the case of a kinetic consideration, mass transfer
        // between phases is realized via source terms there is a
        // balance equation for each component in each phase
        const auto& volVars = elemVolVars[scv];
        ComponentVector componentIntoPhaseMassTransfer[numPhases];
#define FUNKYMASSTRANSFER 0
#if FUNKYMASSTRANSFER
        const Scalar mu_nPhaseNCompEquil = volVars.chemicalPotentialEquil(phase1Idx, comp1Idx) ;   // very 2p2c
        const Scalar mu_wPhaseWCompEquil = volVars.chemicalPotentialEquil(phase0Idx, comp0Idx);    // very 2p2c

        const Scalar mu_wPhaseNComp = volVars.chemicalPotential(phase0Idx, comp1Idx) ;   // very 2p2c
        const Scalar mu_nPhaseWComp = volVars.chemicalPotential(phase1Idx, comp0Idx);    // very 2p2c

        Valgrind::CheckDefined(mu_nPhaseNCompEquil);
        Valgrind::CheckDefined(mu_wPhaseWCompEquil);
        Valgrind::CheckDefined(mu_wPhaseNComp);
        Valgrind::CheckDefined(mu_nPhaseWComp);

        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        const Scalar temperature            = volVars.temperature(phase0Idx);
        const Scalar pn                     = volVars.pressure(phase1Idx);
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
        const Scalar x_wPhaseNCompEquil = volVars.xEquil(phase0Idx, comp1Idx) ;   // very 2p2c
        const Scalar x_nPhaseWCompEquil = volVars.xEquil(phase1Idx, comp0Idx);    // very 2p2c
        Valgrind::CheckDefined(x_wPhaseNCompEquil);
        Valgrind::CheckDefined(x_nPhaseWCompEquil);
        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        const Scalar gradNinWApprox  =  (x[phase0Idx][comp1Idx] - x_wPhaseNCompEquil) / characteristicLength;    // very 2p2c
        const Scalar gradWinNApprox  =  (x[phase1Idx][comp0Idx] - x_nPhaseWCompEquil) / characteristicLength;    // very 2p2c
#endif
        Scalar phaseDensity[numPhases];
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
            phaseDensity[phaseIdx] = volVars.molarDensity(phaseIdx);
        }

        // diffusion coefficients in wetting phase
        const Scalar diffCoeffNinW = volVars.diffusionCoefficient(phase0Idx, comp1Idx) ;
        Valgrind::CheckDefined(diffCoeffNinW);
        // diffusion coefficients in non-wetting phase
        const Scalar diffCoeffWinN = volVars.diffusionCoefficient(phase1Idx, comp0Idx) ;
        Valgrind::CheckDefined(diffCoeffWinN);

        const Scalar factorMassTransfer     = volVars.factorMassTransfer()  ;
        const Scalar awn = volVars.interfacialArea(phase0Idx, phase1Idx);

        const Scalar sherwoodWPhase  = volVars.sherwoodNumber(phase0Idx);
        const Scalar sherwoodNPhase  = volVars.sherwoodNumber(phase1Idx);

        //      actual diffusion is always calculated for eq's 2,3
        //      Eq's 1,4 have to be the same with different sign, because no mass is accumulated in the interface
        //      i.e. automatically conserving mass that mvoes across the interface

        const Scalar nCompIntoWPhase  = - factorMassTransfer * gradNinWApprox * awn * phaseDensity[phase0Idx] * diffCoeffNinW * sherwoodWPhase;
        const Scalar nCompIntoNPhase  = - nCompIntoWPhase ;
        const Scalar wCompIntoNPhase  = - factorMassTransfer * gradWinNApprox * awn * phaseDensity[phase1Idx] * diffCoeffWinN * sherwoodNPhase;
        const Scalar wCompIntoWPhase  = - wCompIntoNPhase ;

        componentIntoPhaseMassTransfer[phase0Idx][comp1Idx] = nCompIntoWPhase;
        componentIntoPhaseMassTransfer[phase1Idx][comp1Idx] = nCompIntoNPhase;
        componentIntoPhaseMassTransfer[phase1Idx][comp0Idx] = wCompIntoNPhase;
        componentIntoPhaseMassTransfer[phase0Idx][comp0Idx] = wCompIntoWPhase;

        Valgrind::CheckDefined(componentIntoPhaseMassTransfer);

        // Actually add the mass transfer to the sources which might
        // exist externally
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                const unsigned int eqIdx = conti0EqIdx + compIdx + phaseIdx*numComponents;
                        source[eqIdx] += componentIntoPhaseMassTransfer[phaseIdx][compIdx] ;

                        using std::isfinite;
                        if (!isfinite(source[eqIdx]))
                            DUNE_THROW(NumericalProblem, "Calculated non-finite source");
            }
        }

        // Call the (kinetic) Energy module, for the source term.
        // it has to be called from here, because the mass transfered has to be known.
        EnergyLocalResidual::computeSourceEnergy(source,
                                                 element,
                                                 fvGeometry,
                                                 elemVolVars,
                                                 scv);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);
        Valgrind::CheckDefined(source);

        return source;
    }
 };

} // end namespace Dumux

#endif
