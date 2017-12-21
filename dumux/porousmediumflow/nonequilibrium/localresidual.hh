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
using NonEquilibriumLocalResidual = NonEquilibriumLocalResidualImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableThermalNonEquilibrium), GET_PROP_VALUE(TypeTag, EnableChemicalNonEquilibrium)>;

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
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
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

    ResidualVector computeFlux(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const SubControlVolumeFace& scvf,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        // get upwind weights into local scope
        ResidualVector flux(0.0);

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

    ResidualVector computeSource(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolume &scv) const
    {
        ResidualVector source(0.0);

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
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { nCompIdx = FluidSystem::nCompIdx } ;
    enum { wCompIdx = FluidSystem::wCompIdx } ;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx} ;
    enum { nPhaseIdx = FluidSystem::nPhaseIdx} ;
    enum { sPhaseIdx = FluidSystem::sPhaseIdx} ;

public:
     using ParentType::ParentType;
    /*!
     * \brief Calculate the storage for all mass balance equations
     *
     *    \param storage The mass of the component within the sub-control volume
     *    \param volVars The volume variables
     */
    ResidualVector computeStorage(const Problem& problem,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars) const
    {
       ResidualVector storage(0.0);

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
    ResidualVector computeFlux(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const SubControlVolumeFace& scvf,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        // get upwind weights into local scope
        ResidualVector flux(0.0);

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


    ResidualVector computeSource(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolume &scv) const
    {
         ResidualVector source(0.0);
        // In the case of a kinetic consideration, mass transfer
        // between phases is realized via source terms there is a
        // balance equation for each component in each phase
        const auto& localScvIdx = scv.indexInElement();
        const auto& volVars = elemVolVars[localScvIdx];
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
        const Scalar diffCoeffNinW = volVars.diffusionCoefficient(wPhaseIdx, nCompIdx) ;
        Valgrind::CheckDefined(diffCoeffNinW);
        // diffusion coefficients in non-wetting phase
        const Scalar diffCoeffWinN = volVars.diffusionCoefficient(nPhaseIdx, wCompIdx) ;
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
