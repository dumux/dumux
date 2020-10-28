// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \ingroup ThermalNonEquilibriumModel
 * \brief This file contains the parts of the local residual to
 *        calculate the heat conservation in the thermal non-equilibrium model.
 */

#ifndef DUMUX_ENERGY_NONEQUILIBRIUM_LOCAL_RESIDUAL_HH
#define DUMUX_ENERGY_NONEQUILIBRIUM_LOCAL_RESIDUAL_HH

#include <cmath>
#include <dumux/common/spline.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

/*!
 * \ingroup ThermalNonEquilibriumModel
 * \brief This file contains the parts of the local residual to
 *        calculate the heat conservation in the thermal non-equilibrium  model.
 */
// forward declaration
template <class TypeTag, int numEnergyEqFluid>
class EnergyLocalResidualNonEquilibrium;

template<class TypeTag>
class EnergyLocalResidualNonEquilibrium<TypeTag, 1/*numEnergyEqFluid*/>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr auto numEnergyEqFluid = ModelTraits::numEnergyEqFluid();
    static constexpr auto numEnergyEqSolid = ModelTraits::numEnergyEqSolid();
    static constexpr auto energyEq0Idx = Indices::energyEq0Idx;
    static constexpr auto energyEqSolidIdx = Indices::energyEqSolidIdx;

    static constexpr auto numPhases = ModelTraits::numFluidPhases();
    static constexpr auto numComponents = ModelTraits::numFluidComponents();

public:
    //! The energy storage in the fluid phase with index phaseIdx
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars,
                                  int phaseIdx)
    {
        //in case we have one energy equation for more than one fluid phase, add up  parts on the             one energy equation
        storage[energyEq0Idx] += volVars.porosity()
                                 * volVars.density(phaseIdx)
                                 * volVars.internalEnergy(phaseIdx)
                                 * volVars.saturation(phaseIdx);

    }


    //! The energy storage in the solid matrix
    static void solidPhaseStorage(NumEqVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars)
    {
        // heat conduction for the fluid phases
        for(int sPhaseIdx = 0; sPhaseIdx<numEnergyEqSolid; ++sPhaseIdx)
        {
            storage[energyEqSolidIdx+sPhaseIdx] += volVars.temperatureSolid()
                                                   * volVars.solidHeatCapacity()
                                                   * volVars.solidDensity()
                                                   * (1.0 - volVars.porosity());
        }
    }


    //! The advective phase energy fluxes
    static void heatConvectionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars,
                                   int phaseIdx)
    {
        auto upwindTerm = [phaseIdx](const auto& volVars)
        { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx)*volVars.enthalpy(phaseIdx); };

        //in case we have one energy equation for more than one fluid phase, add up advective parts on the one energy equation
        flux[energyEq0Idx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

        //now add the diffusive part
        const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
        const auto& elemVolVars = fluxVars.elemVolVars();
        const auto& scvf = fluxVars.scvFace();
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            //no diffusion of the main component, this is a hack to use normal fick's law which computes both diffusions (main and component). We only add the part from the component here
            if (phaseIdx == compIdx)
                continue;
            //we need the upwind enthalpy. Even better would be the componentEnthalpy
            auto enthalpy = 0.0;
            if (diffusiveFluxes[compIdx] > 0)
                enthalpy += insideVolVars.enthalpy(phaseIdx);
            else
                enthalpy += outsideVolVars.enthalpy(phaseIdx);

             //check for the reference system and adapt units of the diffusive flux accordingly.
            if (FluxVariables::MolecularDiffusionType::referenceSystemFormulation() == ReferenceSystemFormulation::massAveraged)
               flux[energyEq0Idx] += diffusiveFluxes[compIdx]*enthalpy;
            else
               flux[energyEq0Idx] += diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx)*enthalpy;
        }
    }

    //! The diffusive energy fluxes
    static void heatConductionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {
        //in case we have one energy equation for more than one fluid phase we use an effective law in the nonequilibrium fourierslaw
        flux[energyEq0Idx] += fluxVars.heatConductionFlux(0);
         //heat conduction for the solid phases
        for(int sPhaseIdx = 0; sPhaseIdx<numEnergyEqSolid; ++sPhaseIdx)
            flux[energyEqSolidIdx+sPhaseIdx] += fluxVars.heatConductionFlux(numPhases + sPhaseIdx);
    }

    /*!
     * \brief heat transfer between the phases for nonequilibrium models
     *
     * \param source The source which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scv The sub-control volume over which we integrate the source term
     */
    static void computeSourceEnergy(NumEqVector& source,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolume &scv)
    {
        // specialization for 2 fluid phases
        const auto& volVars = elemVolVars[scv];
        const Scalar characteristicLength = volVars.characteristicLength()  ;

        // interfacial area
        // Shi & Wang, Transport in porous media (2011)
        const Scalar as = volVars.fluidSolidInterfacialArea();

        // temperature fluid is the same for both fluids
        const Scalar TFluid = volVars.temperatureFluid(0);
        const Scalar TSolid = volVars.temperatureSolid();

        Scalar solidToFluidEnergyExchange ;

        const Scalar fluidConductivity = volVars.fluidThermalConductivity(0) ;

        const Scalar factorEnergyTransfer = volVars.factorEnergyTransfer()  ;

        solidToFluidEnergyExchange = factorEnergyTransfer * (TSolid - TFluid) / characteristicLength * as * fluidConductivity;

        solidToFluidEnergyExchange *= volVars.nusseltNumber(0);

        for(int energyEqIdx = 0; energyEqIdx < numEnergyEqFluid+numEnergyEqSolid; ++energyEqIdx)
        {
            switch (energyEqIdx)
            {
            case 0 :
                source[energyEq0Idx + energyEqIdx] += solidToFluidEnergyExchange;
                break;
            case 1 :
                source[energyEq0Idx + energyEqIdx] -= solidToFluidEnergyExchange;
                break;
            default:
                DUNE_THROW(Dune::NotImplemented,
                        "wrong index");
            } // end switch
        } // end energyEqIdx
    } // end source
};

template<class TypeTag>
class EnergyLocalResidualNonEquilibrium<TypeTag, 2/*numEnergyEqFluid*/>
: public EnergyLocalResidualNonEquilibrium<TypeTag, 1/*numEnergyEqFluid*/>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum { numPhases = ModelTraits::numFluidPhases() };
    enum { numEnergyEqFluid = ModelTraits::numEnergyEqFluid() };
    enum { numEnergyEqSolid = ModelTraits::numEnergyEqSolid() };
    enum { energyEq0Idx = Indices::energyEq0Idx };
    enum { energyEqSolidIdx = Indices::energyEqSolidIdx};
    enum { conti0EqIdx = Indices::conti0EqIdx };

    enum { numComponents = ModelTraits::numFluidComponents() };
    enum { phase0Idx = FluidSystem::phase0Idx};
    enum { phase1Idx = FluidSystem::phase1Idx};
    enum { sPhaseIdx = numPhases};

    static constexpr bool enableChemicalNonEquilibrium = ModelTraits::enableChemicalNonEquilibrium();

public:

    //! The energy storage in the fluid phase with index phaseIdx
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars,
                                  int phaseIdx)
    {
        storage[energyEq0Idx+phaseIdx] += volVars.porosity()
                                          * volVars.density(phaseIdx)
                                          * volVars.internalEnergy(phaseIdx)
                                          * volVars.saturation(phaseIdx);
    }

    //! The advective phase energy fluxes
    static void heatConvectionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars,
                                   int phaseIdx)
    {
        auto upwindTerm = [phaseIdx](const auto& volVars)
        { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx)*volVars.enthalpy(phaseIdx); };

        // in case we have one energy equation for more than one fluid phase, add up advective parts on the one energy equation
        flux[energyEq0Idx+phaseIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

        // add the diffusiv part
        const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
        const auto& elemVolVars = fluxVars.elemVolVars();
        const auto& scvf = fluxVars.scvFace();
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // no diffusion of the main component, this is a hack to use normal fick's law which computes both diffusions (main and component). We only add the part from the component here
            if (phaseIdx == compIdx)
                continue;
            // we need the upwind enthalpy. Even better would be the componentEnthalpy
            auto enthalpy = 0.0;
            if (diffusiveFluxes[compIdx] > 0)
                enthalpy += insideVolVars.enthalpy(phaseIdx);
            else
                enthalpy += outsideVolVars.enthalpy(phaseIdx);
            flux[energyEq0Idx+phaseIdx] += diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx)*enthalpy;
        }
    }

    //! The diffusive energy fluxes
    static void heatConductionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {
        for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            flux[energyEq0Idx+phaseIdx] += fluxVars.heatConductionFlux(phaseIdx);
        }
        for(int sPhaseIdx=0; sPhaseIdx<numEnergyEqSolid; ++sPhaseIdx)
        {
            flux[energyEqSolidIdx+sPhaseIdx] += fluxVars.heatConductionFlux(numPhases + sPhaseIdx);
        }
    }
    /*!
     * \brief Calculates the source term of the equation.
     *
     * \param source The source term which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scv The sub-control volume over which we integrate the source term
     */
    static void computeSourceEnergy(NumEqVector& source,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolume &scv)
    {
        // specialization for 2 fluid phases
        const auto &volVars = elemVolVars[scv];

        const Scalar awn = volVars.interfacialArea(phase0Idx, phase1Idx);
        const Scalar aws = volVars.interfacialArea(phase0Idx, sPhaseIdx);
        const Scalar ans = volVars.interfacialArea(phase1Idx, sPhaseIdx);

        const Scalar Tw = volVars.temperatureFluid(phase0Idx);
        const Scalar Tn = volVars.temperatureFluid(phase1Idx);
        const Scalar Ts = volVars.temperatureSolid();

        const  Scalar lambdaWetting = volVars.fluidThermalConductivity(phase0Idx);
        const  Scalar lambdaNonwetting = volVars.fluidThermalConductivity(phase1Idx);
        const  Scalar lambdaSolid = volVars.solidThermalConductivity();

        const Scalar lambdaWN = harmonicMean(lambdaWetting, lambdaNonwetting);
        const Scalar lambdaWS = harmonicMean(lambdaWetting, lambdaSolid);
        const Scalar lambdaNS = harmonicMean(lambdaNonwetting, lambdaSolid);

        const Scalar characteristicLength = volVars.characteristicLength()  ;
        const Scalar factorEnergyTransfer = volVars.factorEnergyTransfer()  ;

        const Scalar nusseltWN = harmonicMean(volVars.nusseltNumber(phase0Idx), volVars.nusseltNumber(phase1Idx));
        const Scalar nusseltWS = volVars.nusseltNumber(phase0Idx);
        const Scalar nusseltNS = volVars.nusseltNumber(phase1Idx);

        const Scalar wettingToNonwettingEnergyExchange = factorEnergyTransfer * (Tw - Tn) / characteristicLength * awn * lambdaWN * nusseltWN  ;
        const Scalar wettingToSolidEnergyExchange = factorEnergyTransfer * (Tw - Ts) / characteristicLength * aws * lambdaWS * nusseltWS  ;
        const Scalar nonwettingToSolidEnergyExchange = factorEnergyTransfer * (Tn - Ts) / characteristicLength * ans * lambdaNS * nusseltNS  ;

        for(int phaseIdx = 0; phaseIdx < numEnergyEqFluid+numEnergyEqSolid; ++phaseIdx)
        {
            switch (phaseIdx)
            {
            case phase0Idx:
                source[energyEq0Idx + phaseIdx] += ( - wettingToNonwettingEnergyExchange - wettingToSolidEnergyExchange);
                break;
            case phase1Idx:
                source[energyEq0Idx + phaseIdx] += (+ wettingToNonwettingEnergyExchange - nonwettingToSolidEnergyExchange);
                break;
            case sPhaseIdx:
                source[energyEq0Idx + phaseIdx] += (+ wettingToSolidEnergyExchange + nonwettingToSolidEnergyExchange);
                break;
            default:
                DUNE_THROW(Dune::NotImplemented,
                        "wrong index");
            } // end switch


            using std::isfinite;
            if (!isfinite(source[energyEq0Idx + phaseIdx]))
                DUNE_THROW(NumericalProblem, "Calculated non-finite source, " << "Tw="<< Tw << " Tn="<< Tn<< " Ts="<< Ts);
        }// end phases

        // we only need to do this for when there is more than 1 fluid phase
        if (enableChemicalNonEquilibrium)
        {
            // Here comes the catch: We are not doing energy conservation for the whole
            // system, but rather for each individual phase.
            //        -> Therefore the energy fluxes over each phase boundary need be
            //           individually accounted for.
            //        -> Each particle crossing a phase boundary does carry some mass and
            //           thus energy!
            //        -> Therefore, this contribution needs to be added.
            //        -> the particle always brings the energy of the originating phase.
            //        -> Energy advectivly transported into a phase = the moles of a component that go into a  phase
            //           * molMass * enthalpy of the component in the *originating* phase

            const auto& fluidState = volVars.fluidState();

            for(int phaseIdx = 0; phaseIdx < numEnergyEqFluid+numEnergyEqSolid; ++phaseIdx)
            {
                switch (phaseIdx)
                {
                case phase0Idx:
                //sum up the transfered energy by the components into the wetting phase
                    for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        const unsigned int eqIdx = conti0EqIdx + compIdx + phaseIdx*numComponents;
                        source[energyEq0Idx + phaseIdx] += (source[eqIdx]
                                                            * FluidSystem::molarMass(compIdx)
                                                            * FluidSystem::componentEnthalpy(fluidState, phase1Idx, compIdx) );
                    }
                break;
                case phase1Idx:
                //sum up the transfered energy by the components into the nonwetting phase
                    for(int compIdx =0; compIdx<numComponents; ++compIdx)
                    {
                        const unsigned int eqIdx = conti0EqIdx + compIdx + phaseIdx*numComponents;
                        source[energyEq0Idx + phaseIdx] += (source[eqIdx]
                                                            * FluidSystem::molarMass(compIdx)
                                                            *FluidSystem::componentEnthalpy(fluidState, phase0Idx, compIdx));
                    }
                    break;
                case sPhaseIdx:
                    break; // no sorption
                default:
                    DUNE_THROW(Dune::NotImplemented,
                                "wrong index");
                } // end switch
            } // end phases
        } // EnableChemicalNonEquilibrium
    } // end source
};
} // end namespace Dumux

#endif
