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

#ifndef DUMUX_BIOFLUID_SIMPLE_CHEM_SYSTEM_HH
#define DUMUX_BIOFLUID_SIMPLE_CHEM_SYSTEM_HH

// ## The fluid system (`biominsimplechemistry.hh`)
//
// This file contains the __fluidsystem class__ which defines all functions needed to describe the fluids and their properties,
// which are needed to model biomineralization.
//
// [[content]]
//
// ### Include files
// [[details]]
// [[codeblock]]
#include <dumux/material/idealgas.hh>

// we include the base fluid system
#include <dumux/material/fluidsystems/base.hh>
// we include the brine adapter adapting component indices to use the brine fluidsystem
#include "icpcomplexsalinitybrine.hh"

// we include all necessary fluid components
#include <dumux/material/fluidstates/adapter.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/co2tablereader.hh>
#include <dumux/material/components/sodiumion.hh>
#include <dumux/material/components/chlorideion.hh>
#include <dumux/material/components/calciumion.hh>
#include <dumux/material/components/urea.hh>
#include <dumux/material/components/o2.hh>
#include <dumux/material/components/glucose.hh>
#include <examples/biomineralization/material/components/suspendedbiomass.hh>

// we include brine-co2 and h2o-co2 binary coefficients
#include <dumux/material/binarycoefficients/brine_co2.hh>
#include <dumux/material/binarycoefficients/h2o_o2.hh>
// [[/codeblock]]

#include <dumux/material/fluidsystems/nullparametercache.hh>
#include <dumux/common/exceptions.hh>

#include <assert.h>
// [[/details]]
//
// ### The fluidsystem class
// In the BioMinSimpleChemistryFluid fluid system, we define all functions needed to describe the fluids and their properties accounted for in our simulation.
// The simplified biogeochemistry biomineralization fluid system requires the CO2 tables and the H2OType as template parameters.
// We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don`t clash with symbols from other libraries you may want to use in conjunction with Dumux.
// [[codeblock]]
namespace Dumux::FluidSystems {

template <class Scalar,
          class CO2Table,
          class H2OType = Components::TabulatedComponent<Components::H2O<Scalar>> >
class BioMinSimpleChemistryFluid
: public Base<Scalar, BioMinSimpleChemistryFluid<Scalar, CO2Table, H2OType> >
{
    using ThisType = BioMinSimpleChemistryFluid<Scalar, CO2Table, H2OType>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;
    using IdealGas = Dumux::IdealGas<Scalar>;
// [[/codeblock]]
//
// #### Component and phase definitions
// With the following function we define what phases and components will be used by the fluid system and define the indices used to distinguish  phases and components in the course of the simulation.
// [[codeblock]]
public:
    // We use convenient declarations that we derive from the property system
    typedef Components::CO2<Scalar, CO2Table> CO2;
    using H2O = H2OType;
    // export the underlying brine fluid system for the liquid phase, as brine is used as a "pseudo component"
    using Brine = Dumux::FluidSystems::ICPComplexSalinityBrine<Scalar, H2OType>;
    using Na = Components::SodiumIon<Scalar>;
    using Cl = Components::ChlorideIon<Scalar>;
    using Ca = Components::CalciumIon<Scalar>;
    using Urea = Components::Urea<Scalar>;
    using O2 = Components::O2<Scalar>;
    using Glucose = Components::Glucose<Scalar>;
    using SuspendedBiomass = Components::SuspendedBiomass<Scalar>;

    // We define the binary coefficents file, which accounts for the interactions of the main fluids in our setup, water/brine and CO2
    using Brine_CO2 = BinaryCoeff::Brine_CO2<Scalar, CO2Table, true>;

    // the type of parameter cache objects. this fluid system does not
    // cache anything, so it uses Dumux::NullParameterCache
    typedef Dumux::NullParameterCache ParameterCache;

    // We define phase-related indices and properties
    static constexpr int numPhases = 2; // liquid and gas phases
    static constexpr int wPhaseIdx = 0; // index of the liquid phase
    static constexpr int nPhaseIdx = 1; // index of the gas phase
    static constexpr int phase0Idx = wPhaseIdx;
    static constexpr int phase1Idx = nPhaseIdx;

    // the phase names
    static std::string phaseName(int phaseIdx)
    {
        static std::string name[] = {
            "w",
            "n"
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    // We define component-related indices and properties
    static constexpr int numComponents = 9; // H2O, TotalC, Na, Cl, Ca,...
    static constexpr int numSecComponents = 6;

    static constexpr int H2OIdx = 0;
    static constexpr int BrineIdx = 0;
    static constexpr int TCIdx = 1;
    static constexpr int wCompIdx = BrineIdx;
    static constexpr int nCompIdx = TCIdx;
    static constexpr int comp0Idx = wCompIdx;
    static constexpr int comp1Idx = nCompIdx;

    static constexpr int NaIdx  = 2;
    static constexpr int ClIdx  = 3;
    static constexpr int CaIdx  = 4;
    static constexpr int UreaIdx  = 5;
    static constexpr int O2Idx  = 6;
    static constexpr int GlucoseIdx  = 7;
    static constexpr int SuspendedBiomassIdx  = 8;
// [[/codeblock]]
// [[exclude]]

    // We check whether a phase is liquid
    static constexpr bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx != nPhaseIdx;
    }

    // We return whether a phase is gaseous
    static constexpr bool isGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx == nPhaseIdx;
    }

    // We assume ideal mixtures
    static constexpr bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    // We assume compressible fluid phases
    static constexpr bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    // The component names
    static std::string componentName(int compIdx)
    {

        switch (compIdx) {
            case BrineIdx: return H2O::name();
            case TCIdx: return "TotalC";
            case CaIdx: return Ca::name();
            case NaIdx: return Na::name();
            case ClIdx: return Cl::name();
            case SuspendedBiomassIdx: return SuspendedBiomass::name();
            case GlucoseIdx: return Glucose::name();
            case O2Idx: return O2::name();
            case UreaIdx: return Urea::name();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx); break;
        };
    }

    // The component molar masses
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
            case H2OIdx: return H2O::molarMass();
            // actually, the molar mass of brine is only needed for diffusion
            // but since chloride and sodium are accounted for seperately
            // only the molar mass of water is returned.
            case TCIdx: return CO2::molarMass();
            case CaIdx: return Ca::molarMass();
            case NaIdx: return Na::molarMass();
            case ClIdx: return Cl::molarMass();
            case SuspendedBiomassIdx: return SuspendedBiomass::molarMass();
            case GlucoseIdx: return Glucose::molarMass();
            case O2Idx: return O2::molarMass();
            case UreaIdx: return Urea::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx); break;
        };
    }
// [[/exclude]]
//
// #### The Brine Adapter
// With the brine adapter, we link water and the components sodium, chloride, and calcium to form brine, to be able to use the "brine.hh" fluidsystem expecting only water and a single salt.
// Here, we define that the components water, sodium, chloride, and calcium contribute to brine.
// The real work is done by the adapter "icpcomplexsalinitybrine.hh".
// [[codeblock]]
private:
    struct BrineAdapterPolicy
    {
        using FluidSystem = Brine;

        static constexpr int phaseIdx(int brinePhaseIdx) { return wPhaseIdx; }
        static constexpr int compIdx(int brineCompIdx)
        {
            switch (brineCompIdx)
            {
                assert(brineCompIdx == Brine::H2OIdx
                        || brineCompIdx == Brine::NaIdx
                        || brineCompIdx == Brine::ClIdx
                        || brineCompIdx == Brine::CaIdx
                      );
                case Brine::H2OIdx: return H2OIdx;
                case Brine::NaIdx: return NaIdx;
                case Brine::ClIdx: return ClIdx;
                case Brine::CaIdx: return CaIdx;
                default: return 0; // this will never be reached, only needed to suppress compiler warning
            }
        }
    };

    template<class FluidState>
    using BrineAdapter = FluidStateAdapter<FluidState, BrineAdapterPolicy>;
// [[/codeblock]]

    // #### Initializing the fluid system
    // [[codeblock]]
public:
    static void init()
    {
        init(/*startTemp=*/275.15, /*endTemp=*/455.15, /*tempSteps=*/30,
             /*startPressure=*/1e4, /*endPressure=*/99e6, /*pressureSteps=*/500);
    }
    static void init(Scalar startTemp, Scalar endTemp, int tempSteps,
                     Scalar startPressure, Scalar endPressure, int pressureSteps)
    {
        std::cout << "Initializing tables for the pure-water properties.\n";
        H2O::init(startTemp, endTemp, tempSteps,
                            startPressure, endPressure, pressureSteps);
     }
// [[/codeblock]]

    // #### The Component-Phase Interactions
    // In the following, we define the component-phase interactions such as // each component's fugacity coefficient in each phase
    // and each component's and the phases' main components binary diffusion coefficient in the respective phase.
    // [[codeblock]]
    // The component fugacity coefficients
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        if (phaseIdx == nPhaseIdx)
            // use the fugacity coefficients of an ideal gas. the
            // actual value of the fugacity is not relevant, as long
            // as the relative fluid compositions are observed,
            return 1.0;

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (pressure<0)
        {
            typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
            ComponentVector moleFractionsw;
            ComponentVector massFractionsw;

            for (int compIdx = 0; compIdx<numComponents;++compIdx)
            {
                moleFractionsw[compIdx] = fluidState.moleFraction(wPhaseIdx,compIdx);
                massFractionsw[compIdx] = fluidState.massFraction(wPhaseIdx,compIdx);
            }
        }

        assert(temperature > 0);
        assert(pressure > 0);

        // calulate the equilibrium composition for the given
        // temperature and pressure.
        Scalar xgH2O, xlH2O;
        Scalar xlCO2, xgCO2;
        Scalar Xl_Sal =    fluidState.massFraction(wPhaseIdx, NaIdx)                    //Salinity= XNa+XCl+XCa
                        + fluidState.massFraction(wPhaseIdx, ClIdx)
                        + fluidState.massFraction(wPhaseIdx, CaIdx);
        Brine_CO2::calculateMoleFractions(temperature,
                                          pressure,
                                          Xl_Sal,
                                          /*knownPhaseIdx=*/ -1,
                                          xlCO2,
                                          xgH2O);

        xlCO2 = std::max(0.0, std::min(1.0, xlCO2));
        xgH2O = std::max(0.0, std::min(1.0, xgH2O));
        xlH2O = 1.0 - xlCO2;
        xgCO2 = 1.0 - xgH2O;

        // Only H2O, CO2, and O2 in the gas phase
        if (compIdx == BrineIdx) {
            Scalar phigH2O = 1.0;
            return phigH2O * xgH2O / xlH2O;
        }
        else if (compIdx == TCIdx)
        {
            Scalar phigCO2 = 1.0;
            return phigCO2 * xgCO2 / xlCO2;
        }
        else if (compIdx == O2Idx)
        {
            return Dumux::BinaryCoeff::H2O_O2::henry(temperature)/pressure;
        }
        // All other components are assumed not be present in the gas phase
        else
        {
            return 1e-20;
        }
    }

    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    // The binary diffusion coefficients of the components and the phases main component
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             const ParameterCache &paramCache,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIIdx && compIIdx < numComponents);
        assert(0 <= compJIdx && compJIdx < numComponents);

       const Scalar temperature = fluidState.temperature(phaseIdx);
       const Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
        {
            assert(compIIdx == H2OIdx);
            if (compJIdx == TCIdx)
                return Brine_CO2::liquidDiffCoeff(temperature, pressure);
            else if(compJIdx == O2Idx)
                return Dumux::BinaryCoeff::H2O_O2::liquidDiffCoeff(temperature, pressure);
            else //all other components
                return 1.587e-9;  //[m²/s]    //J. Phys. D: Appl. Phys. 40 (2007) 2769-2776 //old Value from Anozie 1e-9
        }
        else
        {
            assert(phaseIdx == nPhaseIdx);
            assert(compIIdx == TCIdx);
            if (compJIdx == H2OIdx)
                return Brine_CO2::gasDiffCoeff(temperature, pressure);
            else if (compJIdx == O2Idx)
                return Dumux::BinaryCoeff::H2O_O2::gasDiffCoeff(temperature, pressure);
            else //all other components
                return 0.0;
        }
    }
    // [[/codeblock]]

    // #### The Fluid Properties
    // In the following, all functions defining the phase properties are given:
    // the density, viscosity, enthalpy, thermal conductivities, and heat capacities
    // of each phase depending on temperature, pressure, and phase composition
    // [[codeblock]]
    // The phase density of the liquid phase is calculated accounting for the impact of solutes in the brine as well as the contribution of CO2.
    // For the gas phase, a mixture of water and CO2 is considered, as solutes do not partition into the gas phase.
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const Scalar T = fluidState.temperature(phaseIdx);

        if (phaseIdx == wPhaseIdx)
        {
            const Scalar pl = fluidState.pressure(phaseIdx);
            const Scalar xlCO2 = fluidState.moleFraction(wPhaseIdx, TCIdx);
            const Scalar xlH2O = fluidState.moleFraction(wPhaseIdx, H2OIdx);
            if(T < 273.15) {
                DUNE_THROW(NumericalProblem,
                        "Liquid density for Brine and CO2 is only "
                        "defined above 273.15K (is" << T << ")");
            }
            if(pl >= 2.5e8) {
            DUNE_THROW(NumericalProblem,
                       "Liquid density for Brine and CO2 is only "
                       "defined below 250MPa (is" << pl << ")");
            }

            // calling the brine adapter for brine density
            Scalar rho_brine = Brine::molarDensity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx)
                       *(H2O::molarMass()*fluidState.moleFraction(wPhaseIdx, H2OIdx)
                         + Na::molarMass()*fluidState.moleFraction(wPhaseIdx, NaIdx)
                         + Cl::molarMass()*fluidState.moleFraction(wPhaseIdx, ClIdx)
                         + Ca::molarMass()*fluidState.moleFraction(wPhaseIdx, CaIdx)
                        );
            // the density of pure water
            Scalar rho_pure = H2O::liquidDensity(T, pl);
            //we also use a private helper function to calculate the liquid density of the same amount of CO2 in pure water
            Scalar rho_lCO2 = liquidDensityWaterCO2_(T, pl, xlH2O, xlCO2);
            // calculating the density contribution of CO2 in pure water
            Scalar contribCO2 = rho_lCO2 - rho_pure;
            // assuming CO2 increases the density for brine by the same amount as for pure water
            return rho_brine + contribCO2;
        }

        else if (phaseIdx == nPhaseIdx)
            return H2O::gasDensity(T, fluidState.partialPressure(nPhaseIdx, H2OIdx))
                   + CO2::gasDensity(T, fluidState.partialPressure(nPhaseIdx, TCIdx));
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    // The phase molar density of the liquid phase is assumed to not change significantly from the molar density of the pure brine.
    // For the gas phase, a mixture of water and CO2 is considered, as solutes do not partition into the gas phase.
    using Base::molarDensity;
    template <class FluidState>
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx)
            return Brine::molarDensity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else if (phaseIdx == nPhaseIdx)
            return H2O::gasMolarDensity(fluidState.temperature(phaseIdx),
                                        fluidState.partialPressure(phaseIdx, H2OIdx))
                   + CO2::gasMolarDensity(fluidState.temperature(phaseIdx),
                                          fluidState.partialPressure(phaseIdx, TCIdx));
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    // The phase viscosities are assumed to not change significantly from those of the pure brine for the liquid and pure CO2 for the gas phase
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == wPhaseIdx)
            return Brine::viscosity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else if (phaseIdx == nPhaseIdx)
        {
            return CO2::gasViscosity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);

    }
    // [[/codeblock]]

    // [[exclude]]
    // non-isothermal phase properties
    // The phase enthalpies are assumed to not change significantly from those of the pure brine for the liquid and pure CO2 for the gas phase
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const Scalar temperature = fluidState.temperature(phaseIdx);
        const Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
            return Brine::enthalpy(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else
            return fluidState.massFraction(nPhaseIdx, H2OIdx)*Brine::gasEnthalpy(temperature, pressure)
                   + fluidState.massFraction(nPhaseIdx, TCIdx)*CO2::gasEnthalpy(temperature, pressure);
    }

    // The phase thermal conductivities are assumed to not change significantly from those of the pure brine for the liquid and pure CO2 for the gas phase
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx)
            return Brine::thermalConductivity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else if (phaseIdx ==nPhaseIdx)
            return CO2::gasThermalConductivity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    // The phase heat capacitiy of the liquid phase is assumed to not change significantly from the heat capacitiy of the pure brine.
    // For the gas phase, a mixture of water and CO2 is considered, as solutes do not partition into the gas phase.
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState, int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
            return Brine::heatCapacity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);

        // We assume NaCl not to be present in the gas phase here
        else if (phaseIdx ==nPhaseIdx)
            return CO2::gasHeatCapacity(T, p)*fluidState.moleFraction(nPhaseIdx, TCIdx)
                   + H2O::gasHeatCapacity(T, p)*fluidState.moleFraction(nPhaseIdx, H2OIdx);

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }


private:

    //a private helper function to calculate the liquid density of the same amount of CO2 as we have in our brine in pure water
    static Scalar liquidDensityWaterCO2_(Scalar temperature,
                                         Scalar pl,
                                         Scalar xlH2O,
                                         Scalar xlCO2)
    {
        const Scalar M_CO2 = CO2::molarMass();
        const Scalar M_H2O = H2O::molarMass();

        const Scalar tempC = temperature - 273.15;        /* tempC : temperature in °C */
        const Scalar rho_pure = H2O::liquidDensity(temperature, pl);
        xlH2O = 1.0 - xlCO2; // xlH2O is available, but in case of a pure gas phase
                             // the value of M_T for the virtual liquid phase can become very large
        const Scalar M_T = M_H2O * xlH2O + M_CO2 * xlCO2;
        const Scalar V_phi =
            (37.51 +
             tempC*(-9.585e-2 +
                    tempC*(8.74e-4 -
                           tempC*5.044e-7))) / 1.0e6;
        return 1/ (xlCO2 * V_phi/M_T + M_H2O * xlH2O / (rho_pure * M_T));
    }
};

} // end namespace Dumux::FluidSystems
// [[/exclude]]
// [[/content]]
#endif
