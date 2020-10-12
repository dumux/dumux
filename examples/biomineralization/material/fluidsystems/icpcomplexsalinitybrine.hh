// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

#ifndef DUMUX_ICP_COMPLEX_SALINITY_BRINE_FLUID_SYSTEM_HH
#define DUMUX_ICP_COMPLEX_SALINITY_BRINE_FLUID_SYSTEM_HH

// ## The brine adapter (`icpcomplexsalinitybrine.hh`)
//
// This file contains the brineadapter class__ which defines
// how to adapt values for communication between
// the "full" fluidsystem accounting for 4 components influencing the brine properties
// and the "brine" fluidsystem assuming water and a single salt determining the brine properties.
//
// [[content]]
//
// ### Include files
// [[details]] includes
// we include all necessary components
#include <dune/common/math.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/sodiumion.hh>
#include <dumux/material/components/chlorideion.hh>
#include <dumux/material/components/calciumion.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/common/exceptions.hh>

#include <dumux/io/name.hh>
// [[/details]]
//
// ### The brine adapter class
// In ICPComplexSalinityBrine we make the transition from 3 components additional to water contributing to salinity to the single component that is assumed in the brine fluid system.
// We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don`t clash with symbols from other libraries you may want to use in conjunction with Dumux.
// [[codeblock]]
namespace Dumux::FluidSystems {

template< class Scalar, class H2OType = Components::TabulatedComponent<Dumux::Components::H2O<Scalar>> >
class ICPComplexSalinityBrine : public Base< Scalar, ICPComplexSalinityBrine<Scalar, H2OType>>
{
    using ThisType = ICPComplexSalinityBrine<Scalar, H2OType>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

public:
    //! export the involved components
    using H2O = H2OType;
    using Na = Components::SodiumIon<Scalar>;
    using Cl = Components::ChlorideIon<Scalar>;
    using Ca = Components::CalciumIon<Scalar>;
    // [[/codeblock]]
    //
    // #### Component and phase definitions
    // With the following function we define what phases and components will be used by the fluid system and define the indices used to distinguish  phases and components in the course of the simulation.
    // [[codeblock]]
    // We use convenient declarations that we derive from the property system.
    static constexpr int numPhases = 1;     //!< Number of phases in the fluid system
    static constexpr int numComponents = 4; //!< Number of components in the fluid system (H2O, Na, Cl, Ca)

    static constexpr int phase0Idx = 0; //!< Index of the first (and only) phase
    static constexpr int liquidPhaseIdx = phase0Idx; //!< The one considered phase is liquid

    static constexpr int H2OIdx = 0;         //!< index of the water component
    static constexpr int NaIdx = 1;        //!< index of the Na component
    static constexpr int ClIdx = 2;        //!< index of the Cl component
    static constexpr int CaIdx = 3;        //!< index of the Ca component
    static constexpr int comp0Idx = H2OIdx;  //!< index of the first component
    static constexpr int comp1Idx = NaIdx; //!< index of the second component
    static constexpr int comp2Idx = ClIdx; //!< index of the third component
    static constexpr int comp3Idx = CaIdx; //!< index of the fourth component

    // the fluid phase index
    static std::string phaseName(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return IOName::liquidPhase();
    }

    // the miscibility
    static constexpr bool isMiscible()
    {
        return false;
    }

    // we do not have a gas phase
    static constexpr bool isGas(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return false;
    }

    // We need a ideal mixture
    static constexpr bool isIdealMixture(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return true;
    }

    // phase compressibility
    static constexpr bool isCompressible(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return H2O::liquidIsCompressible();
    }

    // no gas, thus no ideal gas
    static constexpr bool isIdealGas(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return false; /*we're a liquid!*/
    }
// [[/codeblock]]
// [[exclude]]
    // The component names
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case H2OIdx: return H2O::name();
            case CaIdx: return Ca::name();
            case NaIdx: return Na::name();
            case ClIdx: return Cl::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    // The component molar masses
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
            case H2OIdx: return H2O::molarMass();
            case CaIdx: return Ca::molarMass();
            case NaIdx: return Na::molarMass();
            case ClIdx: return Cl::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    //initializing
     static void init()
     {
         init(/*tempMin=*/273.15,
              /*tempMax=*/623.15,
              /*numTemp=*/100,
              /*pMin=*/-10.,
              /*pMax=*/20e6,
              /*numP=*/200);
     }
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {

        if (H2O::isTabulated)
        {
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            H2O::init(tempMin, tempMax, nTemp, pressMin, pressMax, nPress);
        }
    }
// [[/exclude]]
    //
    // #### The Component-Phase Interactions
    // In the following, we define the component-phase interactions such as // each component's fugacity coefficient in brine
    // and each component's binary diffusion coefficient with water in brine
    // [[codeblock]]
    // The component fugacity coefficients
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == compIdx)
            // We could calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    // The binary diffusion coefficients of the components and the phases main component
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
        {
            if (compIIdx > compJIdx)
            {
                using std::swap;
                swap(compIIdx, compJIdx);
            }
            if (compJIdx == NaIdx || compJIdx == ClIdx || compJIdx == CaIdx)
                return 1.587e-9;  //[m²/s]    //J. Phys. D: Appl. Phys. 40 (2007) 2769-2776
            else
                DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components "
                                                 << compIIdx << " and " << compJIdx
                                                 << " in phase " << phaseIdx);
         }

         DUNE_THROW(Dune::InvalidStateException, "Invalid phase index: " << phaseIdx);
    }
    // [[/codeblock]]
    //
    // #### The Fluid (Brine) Properties
    // In the following, all functions defining the phase properties are given:
    // the density, viscosity, enthalpy, thermal conductivities, and heat capacities
    // of each phase depending on temperature, pressure, and phase composition
    // [[codeblock]]
    // The phase density
    template <class FluidState>
    static Scalar density(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        const Scalar temperature = fluidState.temperature(phaseIdx);
        const Scalar pressure = fluidState.pressure(phaseIdx);
        const Scalar xNaCl = fluidState.massFraction(phaseIdx, NaIdx)
                           + fluidState.massFraction(phaseIdx, ClIdx)
                           + fluidState.massFraction(phaseIdx, CaIdx);

        using std::max;
        const Scalar TempC = temperature - 273.15;
        const Scalar pMPa = pressure/1.0E6;
        const Scalar salinity = max(0.0, xNaCl);

        const Scalar rhow = H2O::liquidDensity(temperature, pressure);
        const Scalar density = rhow + 1000*salinity*(0.668 +
                                                     0.44*salinity +
                                                     1.0E-6*(300*pMPa -
                                                             2400*pMPa*salinity +
                                                             TempC*(80.0 +
                                                                    3*TempC -
                                                                    3300*salinity -
                                                                    13*pMPa +
                                                                    47*pMPa*salinity)));
        assert(density > 0.0);
        return density;
    }



    // The phase viscosity
    template <class FluidState>
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        const Scalar temperature = fluidState.temperature(phaseIdx);
        const Scalar xNaCl = fluidState.massFraction(phaseIdx, NaIdx)
                           + fluidState.massFraction(phaseIdx, ClIdx)
                           + fluidState.massFraction(phaseIdx, CaIdx);

        using std::pow;
        using Dune::power;
        using std::exp;
        using std::max;
        const Scalar T = max(temperature, 275.0);
        const Scalar salinity = max(0.0, xNaCl);

        const Scalar T_C = T - 273.15;
        const Scalar A = ((0.42*power((pow(salinity, 0.8)-0.17), 2)) + 0.045)*pow(T_C, 0.8);
        const Scalar mu_brine = 0.1 + (0.333*salinity) + (1.65+(91.9*salinity*salinity*salinity))*exp(-A); // [cP]
        assert(mu_brine > 0.0);
        return mu_brine/1000.0; // [Pa·s]
    }
// [[/codeblock]]
// [[exclude]]
    // The vapor pressure
    template <class FluidState>
    static Scalar vaporPressure(const FluidState& fluidState, int compIdx)
    {
        if (compIdx == H2OIdx)
        {
            // simplified version of Eq 2.29 in Vishal Jambhekar's Promo
            const Scalar temperature = fluidState.temperature(H2OIdx);
            return H2O::vaporPressure(temperature)/fluidState.massFraction(phase0Idx, H2OIdx);
        }
        else if (compIdx == NaIdx)
            DUNE_THROW(Dune::NotImplemented, "Na::vaporPressure(t)");
        else if (compIdx == ClIdx)
            DUNE_THROW(Dune::NotImplemented, "Cl::vaporPressure(t)");
        else if (compIdx == CaIdx)
            DUNE_THROW(Dune::NotImplemented, "Ca::vaporPressure(t)");
        else
            DUNE_THROW(Dune::NotImplemented, "Invalid component index " << compIdx);
    }

    // The phase enthalpies
    template <class FluidState>
    static Scalar enthalpy(const FluidState& fluidState, int phaseIdx)
    {
        /*Numerical coefficients from PALLISER*/
        static constexpr Scalar f[] = { 2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10 };

        /*Numerical coefficients from MICHAELIDES for the enthalpy of brine*/
        static constexpr Scalar a[4][3] = { { +9633.6, -4080.0, +286.49 },
                                            { +166.58, +68.577, -4.6856 },
                                            { -0.90963, -0.36524, +0.249667E-1 },
                                            { +0.17965E-2, +0.71924E-3, -0.4900E-4 } };

        const Scalar p = fluidState.pressure(phaseIdx);
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar theta = T - 273.15;
        const Scalar salSat = f[0] + f[1]*theta + f[2]*theta*theta + f[3]*theta*theta*theta;

        /*Regularization*/
        using std::min;
        using std::max;
        const Scalar xNaCl = fluidState.massFraction(phaseIdx, NaIdx)
                           + fluidState.massFraction(phaseIdx, ClIdx)
                           + fluidState.massFraction(phaseIdx, CaIdx);
        const Scalar salinity = min(max(xNaCl,0.0), salSat);

        const Scalar hw = H2O::liquidEnthalpy(T, p)/1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        const Scalar h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
                              + ((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* U [kJ/kg] */

        const Scalar m = (1E3/58.44)*(salinity/(1-salinity));

        using Dune::power;
        Scalar d_h = 0;
        for (int i = 0; i<=3; i++)
            for (int j=0; j<=2; j++)
                d_h = d_h + a[i][j] * power(theta, i) * power(m, j);

        /* heat of dissolution for halite according to Michaelides 1971 */
        const Scalar delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without any dissolved gas */
        const Scalar h_ls1 =(1-salinity)*hw + salinity*h_NaCl + salinity*delta_h; /* kJ/kg */
        return h_ls1*1E3; /*J/kg*/
    }
    // The molar density
    template <class FluidState>
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    {
        return density(fluidState, phaseIdx)/fluidState.averageMolarMass(phaseIdx);
    }

    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState& fluidState, int phaseIdx, int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::ICPComplexSalinityBrine::diffusionCoefficient()");
    }

    // The thermal conductivity
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return H2O::liquidThermalConductivity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index: " << phaseIdx);
    }

    // The phase heat capacity
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return H2O::liquidHeatCapacity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
};

} // end namespace Dumux::FluidSystems
// [[/exclude]]
// [[/content]]
#endif
