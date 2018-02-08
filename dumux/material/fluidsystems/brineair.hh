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
 * \brief @copybrief Dumux::FluidSystems::BrineAir
 */
#ifndef DUMUX_BRINE_AIR_SYSTEM_HH
#define DUMUX_BRINE_AIR_SYSTEM_HH

#include <cassert>
#include <dumux/material/idealgas.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/components/brine.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/nacl.hh>
#include <dumux/material/binarycoefficients/brine_air.hh>
#include <dumux/material/binarycoefficients/h2o_air.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/common/exceptions.hh>

namespace Dumux {
namespace FluidSystems {
/*!
 * \ingroup Fluidsystems
 * \brief A compositional two-phase fluid system with a liquid and a gaseous phase
 *        and \f$H_2O\f$, \f$Air\f$ and \f$S\f$ (dissolved minerals) as components.
 *
 *  This fluidsystem is applied by default with the tabulated version of
 *  water of the IAPWS-formulation.
 */
template <class Scalar,
          class H2Otype = TabulatedComponent<Scalar, H2O<Scalar>>,
          bool useComplexRelations=true>
class BrineAir
: public BaseFluidSystem<Scalar, BrineAir<Scalar, H2Otype, useComplexRelations>>
{
    using ThisType = BrineAir<Scalar, H2Otype, useComplexRelations>;
    using Base = BaseFluidSystem<Scalar, ThisType>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:

    using H2O = H2Otype;
    using H2O_Air = BinaryCoeff::H2O_Air;
    using Air = Dumux::Air<Scalar>;
    using Brine_Air = BinaryCoeff::Brine_Air<Scalar, Air>;
    using Brine = Dumux::Brine<Scalar,H2Otype>;
    using NaCl = Dumux::NaCl<Scalar>;

    // the type of parameter cache objects. this fluid system does not
    using ParameterCache = NullParameterCache;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    static const int numSecComponents = 0; //needed to set property not used for 2pncMin model
    static const int numPhases = 2; // liquid and gas phases
    static const int numSPhases = 1;// precipitated solid phases
    static const int lPhaseIdx = 0; // index of the liquid phase
    static const int gPhaseIdx = 1; // index of the gas phase
    static const int sPhaseIdx = 2; // index of the precipitated salt
    static const int wPhaseIdx = lPhaseIdx; // index of the wetting phase
    static const int nPhaseIdx = gPhaseIdx; // index of the non-wetting phase

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static const int wCompIdx = wPhaseIdx;
    static const int nCompIdx = nPhaseIdx;

    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx index of the phase
     */
    static std::string phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case wPhaseIdx: return "liquid";
        case nPhaseIdx: return "gas";
        case sPhaseIdx: return "NaCl";
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx != nPhaseIdx;
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
        // we assume Henry's and Raoult's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!

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
        // ideal gases are always compressible
        if (phaseIdx == nPhaseIdx)
            return true;
        // the water component decides for the liquid phase...
        return H2O::liquidIsCompressible();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        if (phaseIdx == nPhaseIdx)
            return H2O::gasIsIdeal() && Air::gasIsIdeal();
        return false; // not a gas
    }

    /****************************************
     * Component related static parameters
     ****************************************/
    static const int numComponents = 3; // H2O, Air, NaCl
    static const int numMajorComponents = 2;// H2O, Air

    static const int H2OIdx = wCompIdx;//0
    static const int AirIdx = nCompIdx;//1
    static const int NaClIdx  = 2;

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
        case H2OIdx: return H2O::name();
        case AirIdx: return Air::name();
        case NaClIdx:return "NaCl";
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
        case H2OIdx: return H2O::molarMass();
        case AirIdx: return Air::molarMass();
        case NaClIdx:return NaCl::molarMass();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the mass density of the precipitate \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx The index of the precipitated phase to consider
     */
    static Scalar precipitateDensity(int phaseIdx)
    {
        if(phaseIdx != sPhaseIdx)
            DUNE_THROW(Dune::InvalidStateException, "Invalid solid phase index " << sPhaseIdx);
        return NaCl::density();
    }

    /*!
     * \brief Return the saturation vapor pressure of the liquid phase \f$\mathrm{[Pa]}\f$.
     *
     * \param Temperature temperature of the liquid phase
     * \param salinity salinity (salt mass fraction) of the liquid phase
     */
     static Scalar vaporPressure(Scalar Temperature, Scalar salinity) //Vapor pressure dependence on Osmosis
     {
       return vaporPressure_(Temperature,salinity);
     }

    /*!
     * \brief Return the salt specific heat capacity \f$\mathrm{[J/molK]}\f$.
     *
     * \param phaseIdx The index of the precipitated phase to consider
     */
    static Scalar precipitateHeatCapacity(int phaseIdx)
    {
        return NaCl::heatCapacity();
    }

    /*!
     * \brief Return the molar density of the precipitate \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param phaseIdx The index of the precipitated phase to consider
     */
    static Scalar precipitateMolarDensity(int phaseIdx)
     {
        return precipitateDensity(phaseIdx)/molarMass(phaseIdx); //TODO this only works for this specific case here with phaseIdx=compIdx!
     }

    /****************************************
     * thermodynamic relations
     ****************************************/
    /*!
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/800.0,
             /*numTemptempSteps=*/200,
             /*startPressure=*/-10,
             /*endPressure=*/20e6,
             /*pressureSteps=*/200);
    }

   /*!
    * \brief Initialize the fluid system's static parameters using
    *        problem specific temperature and pressure ranges
    *
    * \param tempMin The minimum temperature used for tabulation of water \f$\mathrm{[K]}\f$
    * \param tempMax The maximum temperature used for tabulation of water \f$\mathrm{[K]}\f$
    * \param nTemp The number of ticks on the temperature axis of the  table of water
    * \param pressMin The minimum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
    * \param pressMax The maximum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
    * \param nPress The number of ticks on the pressure axis of the  table of water
    */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                      Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (useComplexRelations)
            std::cout << "Using complex Brine-Air fluid system\n";
        else
            std::cout << "Using fast Brine-Air fluid system\n";

        if (H2O::isTabulated) {
            std::cout << "Initializing tables for the H2O fluid properties ("
                        << nTemp*nPress
                        << " entries).\n";

            H2O::init(tempMin, tempMax, nTemp,
                                pressMin, pressMax, nPress);
        }
    }

    using Base::density;
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param fluidState the fluid state
     * \param phaseIdx index of the phase
     *
     * Equation given in:
     * - Batzle & Wang (1992) \cite batzle1992
     * - cited by: Bachu & Adams (2002)
     *   "Equations of State for basin geofluids" \cite adams2002
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        switch (phaseIdx) {
            case lPhaseIdx:
                return Brine::liquidDensity(temperature,
                              pressure,
                              fluidState.massFraction(lPhaseIdx, NaClIdx));
            case gPhaseIdx:
                return gasDensity_(temperature,
                            pressure,
                            fluidState.moleFraction(gPhaseIdx, H2OIdx));

            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
            }
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

        if (phaseIdx == lPhaseIdx)
        {
            Scalar XNaCl = fluidState.massFraction(lPhaseIdx, NaClIdx);
            result = Brine::liquidViscosity(temperature, pressure, XNaCl);
        }
        else if (phaseIdx == gPhaseIdx)
        {
            result = Air::gasViscosity(temperature, pressure);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);

        Valgrind::CheckDefined(result);
        return result;
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     * \param  fluidState The fluid state
     * \param phaseIdx Index of the phase
     * \param compIdx Index of the component
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
     * For liquids with very low miscibility this boils down to the
     * Henry constant for the solutes and the saturated vapor pressure
     * both divided by phase pressure.
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        assert(T > 0);
        assert(p > 0);

        if (phaseIdx == gPhaseIdx)
            return 1.0;

        else if (phaseIdx == lPhaseIdx)
        {
        if (compIdx == H2OIdx)
            return Brine::vaporPressure(T)/p;
        else if (compIdx == AirIdx)
            return BinaryCoeff::H2O_Air::henry(T)/p;
        else
            return 1/p;
        }
        else
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    template <class FluidState>
    static Scalar kelvinVaporPressure(const FluidState &fluidState,
                                      const int phaseIdx,
                                      const int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::BrineAir::kelvinVaporPressure()");
    }


    using Base::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$\mathrm{i}\f$ and \f$\mathrm{j}\f$ in this phase.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
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

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx) {
            Scalar result = 0.0;
            if(compJIdx == AirIdx)
                result = Brine_Air::liquidDiffCoeff(temperature, pressure);
            else if (compJIdx == NaClIdx)
                result = 0.12e-9; //http://webserver.dmt.upm.es/~isidoro/dat1/Mass%20diffusivity%20data.htm
            else
                DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components "
                                                 << compIIdx << " and " << compJIdx
                                                 << " in phase " << phaseIdx);
            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            assert(phaseIdx == gPhaseIdx);

            if (compIIdx != AirIdx)
            {
                using std::swap;
                swap(compIIdx, compJIdx);
            }

            Scalar result = 0.0;
            if(compJIdx == H2OIdx)
                result = Brine_Air::gasDiffCoeff(temperature, pressure);
            else if (compJIdx == NaClIdx)
                result = 0.12e-9; //Just added to avoid numerical problem. does not have any physical significance
            else
                DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components "
                                                 << compIIdx << " and " << compJIdx
                                                 << " in phase " << phaseIdx);
            Valgrind::CheckDefined(result);
            return result;
        }
    }

    using Base::enthalpy;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     *
     * See:
     * Class 2000
     * Theorie und numerische Modellierung nichtisothermer Mehrphasenprozesse in NAPL-kontaminierten porösen Medien
     * Chapter 2.1.13 Innere Energie, Wäremekapazität, Enthalpie \cite A3:class:2001
     *
     * Formula (2.42):
     * the specific enthalpy of a gas phase result from the sum of (enthalpies*mass fraction) of the components
     * For the calculation of enthalpy of brine we refer to (Michaelides 1981)
     *
     * \note For the phase enthalpy the contribution of gas-molecules in the liquid phase
     *       is neglected. This contribution is probably not big. Somebody would have to
     *       find out the enthalpy of solution for this system. ...
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx)
        {
            Scalar XlNaCl = fluidState.massFraction(phaseIdx, NaClIdx);
            Scalar result = Brine::liquidEnthalpy(T, p, XlNaCl);
                Valgrind::CheckDefined(result);
                return result;
        }
        else
        {
            Scalar XAir = fluidState.massFraction(gPhaseIdx, AirIdx);
            Scalar XH2O = fluidState.massFraction(gPhaseIdx, H2OIdx);

            Scalar result = 0;
            result += XH2O * H2O::gasEnthalpy(T, p);
            result += XAir * Air::gasEnthalpy(T, p);
            Valgrind::CheckDefined(result);
            return result;
        }
    }

   /*!
    * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in a specific phase
    * \param fluidState The fluid state
    * \param phaseIdx The index of the phase
    * \param componentIdx The index of the component
    */
    template <class FluidState>
    static Scalar componentEnthalpy(const FluidState &fluidState,
                                    int phaseIdx,
                                    int componentIdx)
    {
        Scalar T = fluidState.temperature(nPhaseIdx);
        Scalar p = fluidState.pressure(nPhaseIdx);
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);

        if (phaseIdx == wPhaseIdx)
        {
            DUNE_THROW(Dune::NotImplemented, "The component enthalpies in the liquid phase are not implemented.");
        }
        else if (phaseIdx == nPhaseIdx)
        {
            if (componentIdx ==  H2OIdx)
            {
                return H2O::gasEnthalpy(T, p);
            }
            else if (componentIdx == AirIdx)
            {
                return Air::gasEnthalpy(T, p);
            }
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }


    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState An abitrary fluid state
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
        if (phaseIdx == lPhaseIdx)
            return H2O::liquidThermalConductivity(fluidState.temperature(phaseIdx),
                                                  fluidState.pressure(phaseIdx));
        else // gas phase
            return Air::gasThermalConductivity(fluidState.temperature(phaseIdx),
                                               fluidState.pressure(phaseIdx));
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/(kg*K)}\f$.
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * \note The calculation of the isobaric heat capacity is preliminary. A better
     *       description of the influence of the composition on the phase property
     *       has to be found.
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        const Scalar temperature  = fluidState.temperature(phaseIdx);
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
        {
            return H2O::liquidHeatCapacity(temperature, pressure);
        }
        else if (phaseIdx == nPhaseIdx)
        {
            return Air::gasHeatCapacity(temperature, pressure) * fluidState.moleFraction(nPhaseIdx, AirIdx)
                   + H2O::gasHeatCapacity(temperature, pressure) * fluidState.moleFraction(nPhaseIdx, H2OIdx);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the molality of NaCl \f$\mathrm{[mol/m^3]}\f$.
     * \param fluidState An abitrary fluid state
     * \param paramCache parameter cache
     * \param salinity Salinity of brine
     */
    template <class FluidState>
    static Scalar molalityNaCl(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               Scalar salinity)
      {
        return Brine_Air::molalityNaCl(salinity);
      }

private:
    static Scalar gasDensity_(Scalar T, Scalar pg, Scalar xgH2O)
    {
        //Dalton' Law
        const Scalar pH2O = xgH2O*pg;
        const Scalar pAir = pg - pH2O;
        const Scalar gasDensityAir = Air::gasDensity(T, pAir);
        const Scalar gasDensityH2O = H2O::gasDensity(T, pH2O);
        const Scalar gasDensity = gasDensityAir + gasDensityH2O;
        return gasDensity;
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
