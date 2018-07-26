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
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::FluidSystems::BrineCO2
 */
#ifndef DUMUX_BRINE_CO2_SYSTEM_HH
#define DUMUX_BRINE_CO2_SYSTEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/idealgas.hh>
#include <dumux/material/fluidsystems/base.hh>

#include <dumux/material/components/co2.hh>
#include <dumux/material/components/co2tablereader.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/material/binarycoefficients/brine_co2.hh>

namespace Dumux
{
// include the default tables for CO2
#include <dumux/material/components/co2tables.inc>

namespace FluidSystems
{
/*!
 * \ingroup Fluidsystems
 * \brief A compositional fluid with brine and carbon as
 *        components in both, the liquid and the gas (supercritical) phase.
 */
template<class Scalar,
         class CO2Table,
         class H2Otype = Components::TabulatedComponent<Components::H2O<Scalar> >,
         class BrineRawComponent = Components::Brine<Scalar, Components::H2O<Scalar> >,
         class Brinetype = Components::TabulatedComponent<BrineRawComponent> >
class BrineCO2
: public Base<Scalar, BrineCO2<Scalar, CO2Table, H2Otype, BrineRawComponent, Brinetype> >
{
    using ThisType = BrineCO2<Scalar, CO2Table, H2Otype, BrineRawComponent, Brinetype>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;


    using Brine_CO2 = BinaryCoeff::Brine_CO2<Scalar, CO2Table>;

public:
    using ParameterCache = NullParameterCache;
    using H2O = H2Otype;
    using Brine = Brinetype;
    using CO2 = Dumux::Components::CO2<Scalar, CO2Table>;

    static constexpr int numComponents = 2;
    static constexpr int numPhases = 2;

    static constexpr int liquidPhaseIdx = 0; //!< index of the liquid phase
    static constexpr int gasPhaseIdx = 1; //!< index of the gas phase
    static constexpr int phase0Idx = liquidPhaseIdx; //!< index of the first phase
    static constexpr int phase1Idx = gasPhaseIdx; //!< index of the second phase

    static constexpr int comp0Idx = 0;
    static constexpr int comp1Idx = 1;

    static constexpr int BrineIdx = comp0Idx;
    static constexpr int CO2Idx = comp1Idx;

    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        static std::string name[] = {
            std::string("l"),
            std::string("g")
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \brief Returns whether the fluids are miscible
     */
    static constexpr bool isMiscible()
    { return true; }

    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx != gasPhaseIdx;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if Henry's law and Raoult's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        static std::string name[] = {
            Brine::name(),
            CO2::name(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        static const Scalar M[] = {
            Brine::molarMass(),
            CO2::molarMass(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return M[compIdx];
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    // Initializing with salinity and default tables
    static void init(Scalar salinity)
    {
        init(/*startTemp=*/273.15, /*endTemp=*/623.15, /*tempSteps=*/100,
             /*startPressure=*/1e4, /*endPressure=*/40e6, /*pressureSteps=*/200, salinity);
    }

    // Initializing with salinity and custom tables
    static void init(Scalar startTemp, Scalar endTemp, int tempSteps,
                     Scalar startPressure, Scalar endPressure, int pressureSteps,
                     Scalar salinity)
    {
        if(H2O::isTabulated)
        {
            std::cout << "Initializing tables for the pure-water properties.\n";
            H2O::init(startTemp, endTemp, tempSteps,
                                startPressure, endPressure, pressureSteps);
        }
        // set the salinity of brine
        BrineRawComponent::constantSalinity = salinity;

        if(Brine::isTabulated)
        {
            std::cout << "Initializing tables for the brine fluid properties.\n";
            Brine::init(startTemp, endTemp, tempSteps,
                                  startPressure, endPressure, pressureSteps);
        }
    }

    // Initializing with input parameter salinity and default tables
    static void init()
    {
        const auto salinity = getParam<Scalar>("FluidSystem.Salinity", 0.3);
        init(/*startTemp=*/273.15, /*endTemp=*/623.15, /*tempSteps=*/100,
             /*startPressure=*/1e4, /*endPressure=*/40e6, /*pressureSteps=*/200, salinity);
    }

    // Initializing with input parameter salinity and custom tables
    static void init(Scalar startTemp, Scalar endTemp, int tempSteps,
                     Scalar startPressure, Scalar endPressure, int pressureSteps)
    {
        const auto salinity = getParam<Scalar>("FluidSystem.Salinity", 0.3);
        if(H2O::isTabulated)
        {
            std::cout << "Initializing tables for the pure-water properties.\n";
            H2O::init(startTemp, endTemp, tempSteps,
                                startPressure, endPressure, pressureSteps);
        }
        // set the salinity of brine
        BrineRawComponent::constantSalinity = salinity;

        if(Brine::isTabulated)
        {
            std::cout << "Initializing tables for the brine fluid properties.\n";
            Brine::init(startTemp, endTemp, tempSteps,
                                  startPressure, endPressure, pressureSteps);
        }
    }

    using Base::density;
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx) {
            // use normalized composition for to calculate the density
            // (the relations don't seem to take non-normalized
            // compositions too well...)
            using std::min;
            using std::max;
            Scalar xlBrine = min(1.0, max(0.0, fluidState.moleFraction(liquidPhaseIdx, BrineIdx)));
            Scalar xlCO2 = min(1.0, max(0.0, fluidState.moleFraction(liquidPhaseIdx, CO2Idx)));
            Scalar sumx = xlBrine + xlCO2;
            xlBrine /= sumx;
            xlCO2 /= sumx;

            Scalar result = liquidDensity_(temperature,
                                           pressure,
                                           xlBrine,
                                           xlCO2);

            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            assert(phaseIdx == gasPhaseIdx);

            //in the end it comes down to a simplification of just CO2
            return CO2::gasDensity(temperature, pressure);
        }
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the mean molar mass \f$\overline M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{\overline M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx)
        {
            return density(fluidState, phaseIdx)/fluidState.averageMolarMass(phaseIdx);
        }
        else
            return CO2::gasMolarDensity(temperature,pressure);
    }

    using Base::viscosity;
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * \note For the viscosity of the phases the contribution of the minor
     *       component is neglected. This contribution is probably not big, but somebody
     *       would have to find out its influence.
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        Scalar result = 0;

        if (phaseIdx == liquidPhaseIdx)
            result = Brine::liquidViscosity(temperature, pressure);
        else
            result = CO2::gasViscosity(temperature, pressure);

        Valgrind::CheckDefined(result);
        return result;
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     *
     * The fugacity coefficient \f$\mathrm{\phi^\kappa_\alpha}\f$ of
     * component \f$\mathrm{\kappa}\f$ in phase \f$\mathrm{\alpha}\f$ is connected to
     * the fugacity \f$\mathrm{f^\kappa_\alpha}\f$ and the component's mole
     * fraction \f$\mathrm{x^\kappa_\alpha}\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     \f]
     * where \f$\mathrm{p_\alpha}\f$ is the pressure of the fluid phase.
     *
     * The fugacity itself is just an other way to express the
     * chemical potential \f$\mathrm{\zeta^\kappa_\alpha}\f$ of the component:
     *
     * \f[
     f^\kappa_\alpha := \exp\left\{\frac{\zeta^\kappa_\alpha}{k_B T_\alpha} \right\}
     \f]
     * where \f$\mathrm{k_B}\f$ is Boltzmann's constant.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        if (phaseIdx == gasPhaseIdx)
            // use the fugacity coefficients of an ideal gas. the
            // actual value of the fugacity is not relevant, as long
            // as the relative fluid compositions are observed,
            return 1.0;

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        assert(temperature > 0);
        assert(pressure > 0);

        // calulate the equilibrium composition for the given
        // temperature and pressure.
        Scalar xlH2O, xgH2O;
        Scalar xlCO2, xgCO2;
        Brine_CO2::calculateMoleFractions(temperature,
                                          pressure,
                                          BrineRawComponent::salinity,
                                          /*knowgasPhaseIdx=*/-1,
                                          xlCO2,
                                          xgH2O);

        // normalize the phase compositions
        using std::min;
        using std::max;
        xlCO2 = max(0.0, min(1.0, xlCO2));
        xgH2O = max(0.0, min(1.0, xgH2O));

        xlH2O = 1.0 - xlCO2;
        xgCO2 = 1.0 - xgH2O;

        if (compIdx == BrineIdx) {
            Scalar phigH2O = 1.0;
            return phigH2O * xgH2O / xlH2O;
        }

        assert(compIdx == CO2Idx);

        Scalar phigCO2 = 1.0;
        return phigCO2 * xgCO2 / xlCO2;
    }

    /*!
     * \brief Returns the equilibrium concentration of the dissolved component
     *        in a phase.
     * \param fluidState An arbitrary fluid state
     * \param paramCache Parameter cache
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar equilibriumMoleFraction(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        assert(temperature > 0);
        assert(pressure > 0);

        Scalar xgH2O;
        Scalar xlCO2;

        // calulate the equilibrium composition for the given
        // temperature and pressure.
        Brine_CO2::calculateMoleFractions(temperature,
                                                  pressure,
                                                  BrineRawComponent::constantSalinity,
                                                  /*knowgasPhaseIdx=*/-1,
                                                  xlCO2,
                                                  xgH2O);

        if(phaseIdx == gasPhaseIdx)
        {
            return xgH2O;
        }
        else
        {
            return xlCO2;
        }
    }


    using Base::diffusionCoefficient;
    /*!
     * \brief Calculate the molecular diffusion coefficient for a
     *        component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     *
     * Molecular diffusion of a compoent \f$\mathrm{\kappa}\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \textbf{grad} mu_\kappa \f]
     *
     * where \f$\mathrm{\mu_\kappa}\f$ is the component's chemical potential,
     * \f$D\f$ is the diffusion coefficient and \f$\mathrm{J}\f$ is the
     * diffusive flux. \f$\mathrm{mu_\kappa}\f$ is connected to the component's
     * fugacity \f$\mathrm{f_\kappa}\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$\mathrm{p_\alpha}\f$ and \f$\mathrm{T_\alpha}\f$ are the fluid phase'
     * pressure and temperature.
     *
     * Maybe see http://www.ddbst.de/en/EED/PCP/DIF_C1050.php
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given the phase compositions, return the binary
     *        diffusion coefficient \f$\mathrm{[m^2/s]}\f$ of two components in a phase.
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx Index of the component i
     * \param compJIdx Index of the component j
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIIdx && compIIdx < numComponents);
        assert(0 <= compJIdx && compJIdx < numComponents);

        if (compIIdx > compJIdx)
        {
            using std::swap;
            swap(compIIdx, compJIdx);
        }

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == liquidPhaseIdx) {
            assert(compIIdx == BrineIdx);
            assert(compJIdx == CO2Idx);

            Scalar result = Brine_CO2::liquidDiffCoeff(temperature, pressure);
            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            assert(phaseIdx == gasPhaseIdx);
            assert(compIIdx == BrineIdx);
            assert(compJIdx == CO2Idx);

            Scalar result = Brine_CO2::gasDiffCoeff(temperature, pressure);
            Valgrind::CheckDefined(result);
            return result;
        }
    }

    using Base::enthalpy;
    /*!
     * \brief Given the phase composition, return the specific
     *        phase enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx) {
            Scalar XlCO2 = fluidState.massFraction(phaseIdx, CO2Idx);

            Scalar result = liquidEnthalpyBrineCO2_(temperature,
                                                    pressure,
                                                    BrineRawComponent::constantSalinity,
                                                    XlCO2);
            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            Scalar result = 0;
            result +=
                Brine::gasEnthalpy(temperature, pressure) *
                fluidState.massFraction(gasPhaseIdx, BrineIdx);
            result +=
                CO2::gasEnthalpy(temperature, pressure) *
                fluidState.massFraction(gasPhaseIdx, CO2Idx);
            Valgrind::CheckDefined(result);
            return result;
        }
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * \note For the thermal conductivity of the phases the contribution of the minor
     *       component is neglected. This contribution is probably not big, but somebody
     *       would have to find out its influence.
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
        {
            return H2O::liquidThermalConductivity(fluidState.temperature(phaseIdx),
                                                  fluidState.pressure(phaseIdx));
        }
        else // gas phase
            return CO2::gasThermalConductivity(fluidState.temperature(phaseIdx),
                                               fluidState.pressure(phaseIdx));
    }

    using Base::heatCapacity;
    /*!
     * \copybrief Base::heatCapacity
     *
     * \note We employ the heat capacity of the pure phases.
     *
     * \todo Implement heat capacity for gaseous CO2
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        if(phaseIdx == liquidPhaseIdx)
            return H2O::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
        else
            return CO2::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
    }

private:

    /***********************************************************************/
    /*                                                                     */
    /* Total brine density with dissolved CO2                              */
    /* rho_{b,CO2} = rho_w + contribution(salt) + contribution(CO2)        */
    /*                                                                     */
    /***********************************************************************/
    static Scalar liquidDensity_(Scalar T,
                                 Scalar pl,
                                 Scalar xlH2O,
                                 Scalar xlCO2)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pl);
        Valgrind::CheckDefined(xlH2O);
        Valgrind::CheckDefined(xlCO2);

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

        Scalar rho_brine = Brine::liquidDensity(T, pl);
        Scalar rho_pure = H2O::liquidDensity(T, pl);
        Scalar rho_lCO2 = liquidDensityWaterCO2_(T, pl, xlH2O, xlCO2);
        Scalar contribCO2 = rho_lCO2 - rho_pure;

        return rho_brine + contribCO2;
    }

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

    static Scalar liquidEnthalpyBrineCO2_(Scalar T,
                                          Scalar p,
                                          Scalar S,
                                          Scalar X_CO2_w)
    {
        /* X_CO2_w : mass fraction of CO2 in brine */

        const Scalar h_ls1 = BrineRawComponent::liquidEnthalpy(T, p, S)/1E3; /* J/kg */

        /* heat of dissolution for CO2 according to Fig. 6 in Duan and Sun 2003. (kJ/kg)
           In the relevant temperature ranges CO2 dissolution is
           exothermal */
        const Scalar delta_hCO2 = (-57.4375 + T * 0.1325) * 1000/44;

        const Scalar hw = H2O::liquidEnthalpy(T, p) /1E3; /* kJ/kg */

        /* enthalpy contribution of CO2 (kJ/kg) */
        const Scalar hg = CO2::liquidEnthalpy(T, p)/1E3 + delta_hCO2;

        /* Enthalpy of brine with dissolved CO2 */
        const Scalar h_ls = (h_ls1 - X_CO2_w*hw + hg*X_CO2_w)*1E3; /*J/kg*/

        return h_ls;
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
