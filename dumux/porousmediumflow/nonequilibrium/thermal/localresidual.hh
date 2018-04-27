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
 * \ingroup PorousmediumThermalNonEquilibriumModel
 * \brief This file contains the parts of the local residual to
 *        calculate the heat conservation in the thermal non-equilibrium model.
 */
#ifndef DUMUX_ENERGY_NONEQUILIBRIUM_LOCAL_RESIDUAL_HH
#define DUMUX_ENERGY_NONEQUILIBRIUM_LOCAL_RESIDUAL_HH

#include <cmath>
#include <dumux/common/spline.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumThermalNonEquilibriumModel
 * \brief This file contains the parts of the local residual to
 *        calculate the heat conservation in the thermal non-equilibrium  model.
 */
// forward declaration
template <class TypeTag, int numEnergyEqFluid>
class EnergyLocalResidualNonEquilibrium;

template<class TypeTag>
class EnergyLocalResidualNonEquilibrium<TypeTag, 1/*numEnergyEqFluid*/>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;

    enum { numEnergyEqFluid = ModelTraits::numEnergyEqFluid() };
    enum { numEnergyEqSolid = ModelTraits::numEnergyEqSolid() };
    enum { energyEq0Idx = Indices::energyEq0Idx };
    enum { energyEqSolidIdx = Indices::energyEqSolidIdx};

    enum { numComponents    = ModelTraits::numComponents() };
    enum { phase0Idx        = FluidSystem::phase0Idx};
    enum { phase1Idx        = FluidSystem::phase1Idx};

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
         //heat conduction for the fluid phases
       for(int sPhaseIdx=0; sPhaseIdx<numEnergyEqSolid; ++sPhaseIdx)
       {
            storage[energyEqSolidIdx+sPhaseIdx] += volVars.temperatureSolid()
                                                * volVars.solidHeatCapacity()
                                                * volVars.solidDensity()
                                                * (1.0 - volVars.porosity());
       }
    }

     // this is to make nonequilibrium work with compositional local residual, compositional calls that for non-isothermal models
    static void heatConvectionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars,
                                   int phaseIdx)
    {}


   //! The advective phase energy fluxes
    static void heatConvectionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars,
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   int phaseIdx)
    {
        auto upwindTerm = [phaseIdx](const auto& volVars)
        { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx)*volVars.enthalpy(phaseIdx); };

        //in case we have one energy equation for more than one fluid phase, add up advective parts on the one energy equation
        flux[energyEq0Idx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

        //now add the diffusive part
        const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        auto insideEnthalpy = insideVolVars.enthalpy(phaseIdx);
        auto outsideEnthalpy = outsideVolVars.enthalpy(phaseIdx);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            //no diffusion of the main component, this is a hack to use normal fick's law which computes both diffusions (main and component). We only add the part from the component here
            if (phaseIdx == compIdx)
                continue;
            //we need the upwind enthapy. Even better would be the componentEnthalpy
            auto enthalpy = 0.0;
            if (diffusiveFluxes[compIdx] > 0)
                enthalpy += insideEnthalpy;
            else
                enthalpy += outsideEnthalpy;
            flux[energyEq0Idx] += diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx)*enthalpy;
        }
    }

    //! The diffusive energy fluxes
    static void heatConductionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {
        //in case we have one energy equation for more than one fluid phase we use an effective law in the nonequilibrium fourierslaw
         flux[energyEq0Idx] += fluxVars.heatConductionFlux(0);
         //heat conduction for the fluid phases
       for(int sPhaseIdx=0; sPhaseIdx<numEnergyEqSolid; ++sPhaseIdx)
       {
            flux[energyEqSolidIdx+sPhaseIdx] += fluxVars.heatConductionFlux(sPhaseIdx);
       }
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param scv The sub-control volume over which we integrate the source term
     */
    static void computeSourceEnergy(NumEqVector& source,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolume &scv)
    {
        //specialization for 2 fluid phases
        const auto& volVars = elemVolVars[scv];
        const auto& fs = volVars.fluidState() ;
        const Scalar characteristicLength = volVars.characteristicLength()  ;

        //interfacial area
        // Shi & Wang, Transport in porous media (2011)
        const Scalar as = 6.0 * (1.0-volVars.porosity()) / characteristicLength ;

        //temperature fluid is the same for both fluids
        const Scalar TFluid     = volVars.temperature(0);
        const Scalar TSolid     = volVars.temperatureSolid();

        const Scalar satW       = fs.saturation(phase0Idx) ;
        const Scalar satN       = fs.saturation(phase1Idx) ;

        const Scalar eps = 1e-6 ;
        Scalar solidToFluidEnergyExchange ;

        Scalar fluidConductivity ;
            if (satW < 1.0 - eps)
                fluidConductivity = volVars.fluidThermalConductivity(phase1Idx) ;
            else if (satW >= 1.0 - eps)
                fluidConductivity = volVars.fluidThermalConductivity(phase0Idx) ;
            else
                DUNE_THROW(Dune::NotImplemented,
                        "wrong range");

        const Scalar factorEnergyTransfer   = volVars.factorEnergyTransfer()  ;

        solidToFluidEnergyExchange = factorEnergyTransfer * (TSolid - TFluid) / characteristicLength * as * fluidConductivity ;
        const Scalar epsRegul = 1e-3 ;

        if (satW < (0 + eps) )
        {
                solidToFluidEnergyExchange *=  volVars.nusseltNumber(phase1Idx) ;
        }
        else if ( (satW >= 0 + eps) and (satW < 1.0-eps) )
        {
            solidToFluidEnergyExchange *=  (volVars.nusseltNumber(phase1Idx) * satN );
            Scalar qBoil ;
            if (satW<=epsRegul)
            {// regularize
                typedef Dumux::Spline<Scalar> Spline;
                    Spline sp(0.0, epsRegul,                           // x1, x2
                        QBoilFunc(volVars, 0.0), QBoilFunc(volVars, epsRegul),       // y1, y2
                        0.0,dQBoil_dSw(volVars, epsRegul));    // m1, m2

                qBoil = sp.eval(satW) ;
            }

            else if (satW>= (1.0-epsRegul) )
            {// regularize
                typedef Dumux::Spline<Scalar> Spline;
                Spline sp(1.0-epsRegul, 1.0,    // x1, x2
                        QBoilFunc(volVars, 1.0-epsRegul), 0.0,    // y1, y2
                        dQBoil_dSw(volVars, 1.0-epsRegul), 0.0 );      // m1, m2

                qBoil = sp.eval(satW) ;
            }
            else
            {
                qBoil = QBoilFunc(volVars, satW)  ;
             }

            solidToFluidEnergyExchange += qBoil;
        }
        else if (satW >= 1.0-eps)
        {
            solidToFluidEnergyExchange *=  volVars.nusseltNumber(phase0Idx) ;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                    "wrong range");

        using std::isfinite;
        if (!isfinite(solidToFluidEnergyExchange))
            DUNE_THROW(NumericalProblem, "Calculated non-finite source, " << "TFluid="<< TFluid << " TSolid="<< TSolid);

        for(int energyEqIdx =0; energyEqIdx<numEnergyEqFluid+numEnergyEqSolid; ++energyEqIdx)
        {
            switch (energyEqIdx)
            {
            case 0 :
                source[energyEq0Idx + energyEqIdx] += solidToFluidEnergyExchange;
                break;
            case 1 :
                source[energyEq0Idx + energyEqIdx] -=   solidToFluidEnergyExchange;
                break;
            default:
                DUNE_THROW(Dune::NotImplemented,
                        "wrong index");
            } // end switch
        }// end energyEqIdx
    }// end source


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
        const auto& fs = volVars.fluidState() ;
        const Scalar g( 9.81 ) ;
        const Scalar gamma(0.0589) ;
        const Scalar TSolid     = volVars.temperatureSolid();
        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        using std::pow;
        const Scalar as = 6.0 * (1.0-volVars.porosity()) / characteristicLength ;
        const Scalar mul = fs.viscosity(phase0Idx) ;
        const Scalar deltahv = fs.enthalpy(phase1Idx) - fs.enthalpy(phase0Idx);
        const Scalar deltaRho = fs.density(phase0Idx) - fs.density(phase1Idx) ;
        const Scalar firstBracket = pow(g * deltaRho / gamma, 0.5);
        const Scalar cp = FluidSystem::heatCapacity(fs, phase0Idx) ;
        // This use of Tsat is only justified if the fluid is always boiling (tsat equals boiling conditions)
        // If a different state is to be simulated, please use the actual fluid temperature instead.
        const Scalar Tsat = FluidSystem::vaporTemperature(fs, phase1Idx ) ;
        const Scalar deltaT = TSolid - Tsat ;
        const Scalar secondBracket = pow( (cp *deltaT / (0.006 * deltahv)  ) , 3.0 ) ;
        const Scalar Prl = volVars.prandtlNumber(phase0Idx) ;
        const Scalar thirdBracket = pow( 1/Prl , (1.7/0.33) );
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

template<class TypeTag>
class EnergyLocalResidualNonEquilibrium<TypeTag, 2 /*numEnergyEqFluid*/>
: public EnergyLocalResidualNonEquilibrium<TypeTag, 1 /*numEnergyEqFluid*/>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SolidSystem = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;

    enum { numPhases        = ModelTraits::numPhases() };
    enum { numEnergyEqFluid = ModelTraits::numEnergyEqFluid() };
    enum { numEnergyEqSolid = ModelTraits::numEnergyEqSolid() };
    enum { energyEq0Idx = Indices::energyEq0Idx };
    enum { energyEqSolidIdx = Indices::energyEqSolidIdx};
    enum { conti0EqIdx = Indices::conti0EqIdx };

    enum { numComponents    = ModelTraits::numComponents() };
    enum { phase0Idx        = FluidSystem::phase0Idx};
    enum { phase1Idx        = FluidSystem::phase1Idx};
    enum { sPhaseIdx        = numPhases};

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
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   int phaseIdx)
    {
        auto upwindTerm = [phaseIdx](const auto& volVars)
        { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx)*volVars.enthalpy(phaseIdx); };

        //in case we have one energy equation for more than one fluid phase, add up advective parts on the one energy equation
        flux[energyEq0Idx+phaseIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

        //add the diffusiv part
        const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        auto insideEnthalpy = insideVolVars.enthalpy(phaseIdx);
        auto outsideEnthalpy = outsideVolVars.enthalpy(phaseIdx);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            //no diffusion of the main component, this is a hack to use normal fick's law which computes both diffusions (main and component). We only add the part from the component here
            if (phaseIdx == compIdx)
                continue;
            //we need the upwind enthapy. Even better would be the componentEnthalpy
            auto enthalpy = 0.0;
            if (diffusiveFluxes[compIdx] > 0)
                enthalpy += insideEnthalpy;
            else
                enthalpy += outsideEnthalpy;
            flux[energyEq0Idx+phaseIdx] += diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx)*enthalpy;
        }
    }

    //! The diffusive energy fluxes
    static void heatConductionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {

        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
        {
            flux[energyEq0Idx+phaseIdx] += fluxVars.heatConductionFlux(phaseIdx);
        }
       for(int sPhaseIdx=0; sPhaseIdx<numEnergyEqSolid; ++sPhaseIdx)
       {
            flux[energyEqSolidIdx+sPhaseIdx] += fluxVars.heatConductionFlux(sPhaseIdx);
       }
    }
    /*!
     * \brief Calculate the source term of the equation
     *
     * \param scv The sub-control volume over which we integrate the source term
     */
    static void computeSourceEnergy(NumEqVector& source,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolume &scv)
    {
        //specialization for 2 fluid phases
        const auto &volVars = elemVolVars[scv];

        const Scalar awn = volVars.interfacialArea(phase0Idx, phase1Idx);
        const Scalar aws = volVars.interfacialArea(phase0Idx, sPhaseIdx);
        const Scalar ans = volVars.interfacialArea(phase1Idx, sPhaseIdx);

        const Scalar Tw = volVars.temperature(phase0Idx);
        const Scalar Tn = volVars.temperature(phase1Idx);
        const Scalar Ts = volVars.temperatureSolid();

        const  Scalar lambdaWetting     = volVars.fluidThermalConductivity(phase0Idx);
        const  Scalar lambdaNonWetting  = volVars.fluidThermalConductivity(phase1Idx);
        const  Scalar lambdaSolid       = volVars.solidThermalConductivity();

        const Scalar lambdaWN      = harmonicMean(lambdaWetting, lambdaNonWetting);
        const Scalar lambdaWS      = harmonicMean(lambdaWetting, lambdaSolid);
        const Scalar lambdaNS      = harmonicMean(lambdaNonWetting, lambdaSolid);

        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        const Scalar factorEnergyTransfer   = volVars.factorEnergyTransfer()  ;

        const Scalar nusseltWN      = harmonicMean(volVars.nusseltNumber(phase0Idx), volVars.nusseltNumber(phase1Idx));
        const Scalar nusseltWS      = volVars.nusseltNumber(phase0Idx);
        const Scalar nusseltNS      = volVars.nusseltNumber(phase1Idx);

        const Scalar wettingToNonWettingEnergyExchange = factorEnergyTransfer * (Tw - Tn) / characteristicLength * awn * lambdaWN * nusseltWN  ;
        const Scalar wettingToSolidEnergyExchange      = factorEnergyTransfer * (Tw - Ts) / characteristicLength * aws * lambdaWS * nusseltWS  ;
        const Scalar nonWettingToSolidEnergyExchange   = factorEnergyTransfer * (Tn - Ts) / characteristicLength * ans * lambdaNS * nusseltNS  ;

        for(int phaseIdx =0; phaseIdx<numEnergyEqFluid+numEnergyEqSolid; ++phaseIdx){
            switch (phaseIdx)
            {
            case phase0Idx:
                source[energyEq0Idx + phaseIdx] +=  ( - wettingToNonWettingEnergyExchange - wettingToSolidEnergyExchange);
                break;
            case phase1Idx:
                source[energyEq0Idx + phaseIdx] +=  (+ wettingToNonWettingEnergyExchange - nonWettingToSolidEnergyExchange);
                break;
            case sPhaseIdx:
                source[energyEq0Idx + phaseIdx] +=  (+ wettingToSolidEnergyExchange + nonWettingToSolidEnergyExchange);
                break;
            default:
                DUNE_THROW(Dune::NotImplemented,
                        "wrong index");
            } // end switch


            using std::isfinite;
            if (!isfinite(source[energyEq0Idx + phaseIdx]))
                DUNE_THROW(NumericalProblem, "Calculated non-finite source, " << "Tw="<< Tw << " Tn="<< Tn<< " Ts="<< Ts);

        }// end phases

        //we only need to do this for when there is more than 1 fluid phase
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

            for(int phaseIdx =0; phaseIdx<numEnergyEqFluid+numEnergyEqSolid; ++phaseIdx)
            {
                switch (phaseIdx)
                {
                case phase0Idx:
                //sum up the transfered energy by the components into the wetting phase
                    for(int compIdx =0; compIdx<numComponents; ++compIdx)
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
            }// end phases
        }// EnableChemicalNonEquilibrium
    }// end source
};
} // end namespace Dumux

#endif //
