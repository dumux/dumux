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
 *
 * \brief A fluid system with a liquid and a gaseous phase and \f$H_2O\f$, \f$Air\f$ and \f$S\f$ (dissolved minerals) as components.
 *
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

#ifdef DUMUX_PROPERTIES_HH
#include <dumux/common/propertysystem.hh>
#include <dumux/common/basicproperties.hh>
#endif

namespace Dumux
{
namespace FluidSystems
{
/*!
 * \ingroup Fluidsystems
 *
 * \brief A compositional two-phase fluid system with water, air and dissolved salt as components in both, the liquid and the gas phase.
 *
 *  This fluidsystem is applied by default with the tabulated version of
 *  water of the IAPWS-formulation.
 *
 *  To change the component formulation (i.e. to use nontabulated or
 *  incompressible water), or to switch on verbosity of tabulation,
 *  specify the water formulation via template arguments or via the property
 *  system, as described in the TypeTag Adapter at the end of the file.
 *
            // Select fluid system
            SET_PROP(TestDecTwoPTwoCProblem, FluidSystem)
            {
                typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

                typedef Dumux::FluidSystems::BrineAir<Scalar, Dumux::SimpleH2O<Scalar>> type;
            };

 *   Also remember to initialize tabulated components (FluidSystem::init()), while this
 *   is not necessary for non-tabularized ones.
 *
 * This FluidSystem can be used without the PropertySystem that is applied in Dumux,
 * as all Parameters are defined via template parameters. Hence it is in an
 * additional namespace Dumux::FluidSystem::.
 * An adapter class using Dumux::FluidSystem<TypeTag> is also provided
 * at the end of this file.
 */
template <class Scalar,
          class H2Otype = Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar>>,
          bool useComplexRelations=true>
class BrineAir
: public BaseFluidSystem<Scalar, BrineAir<Scalar, H2Otype, useComplexRelations>>
{
    typedef BrineAir<Scalar, H2Otype, useComplexRelations> ThisType;
    typedef BaseFluidSystem <Scalar, ThisType> Base;

    typedef Dumux::IdealGas<Scalar> IdealGas;

public:

    typedef H2Otype H2O;
    typedef Dumux::BinaryCoeff::H2O_Air H2O_Air;
    typedef Dumux::Air<Scalar> Air;
    typedef Dumux::BinaryCoeff::Brine_Air<Scalar, Air> Brine_Air;
    typedef Dumux::Brine<Scalar, H2Otype> Brine;
    typedef Dumux::NaCl<Scalar> NaCl;

    // the type of parameter cache objects. this fluid system does not
    typedef Dumux::NullParameterCache ParameterCache;

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
    static const char *phaseName(int phaseIdx)
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
     * are indepent on the fluid composition. This assumtion is true
     * if Henry's law and Rault's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
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
    static bool isCompressible(int phaseIdx)
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
    static const char *componentName(int compIdx)
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
     * \brief Return the salt specific heat capacity \f$\mathrm{[J/(kg K)]}\f$.
     *
     * \param phaseIdx The index of the precipitated phase to consider
     */
     DUNE_DEPRECATED_MSG("saltSpecificHeatCapacity(int phaseIdx) is deprecated. Use precipitateSpecificHeatCapacity(int sPhaseIdx) instead.")
     static Scalar saltSpecificHeatCapacity(int phaseIdx)//Specific heat capacity per unit mole of solid salt phase (J/Kkg)
    {
        return 36.79/molarMass(phaseIdx);
    }

     /*!
     * \brief Return the salt specific heat capacity \f$\mathrm{[J/molK]}\f$.
     *
     * \param sPhaseIdx The index of the precipitated phase to consider
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
        return precipitateDensity(phaseIdx)/molarMass(phaseIdx);
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

     /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx index of the phase
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     * \param fluidState the fluid state
     *
     * Equation given in:
     *                         - Batzle & Wang (1992) \cite batzle1992 <BR>
     *                         - cited by: Bachu & Adams (2002)
     *                           "Equations of State for basin geofluids" \cite adams2002
     */
    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        switch (phaseIdx) {
            case lPhaseIdx:
                return liquidDensity_(temperature,
                                pressure,
                                fluidState.moleFraction(lPhaseIdx, AirIdx),
                                fluidState.moleFraction(lPhaseIdx, H2OIdx),
                                fluidState.massFraction(lPhaseIdx, NaClIdx));
            case gPhaseIdx:
                return gasDensity_(temperature,
                            pressure,
                            fluidState.moleFraction(gPhaseIdx, H2OIdx));

            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
            }
      }

    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * Equation given in:
     *                         - Batzle & Wang (1992) \cite batzle1992 <BR>
     *                         - cited by: Bachu & Adams (2002)
     *                           "Equations of State for basin geofluids" \cite adams2002
     */
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {

        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx)
        {
            // assume pure brine for the liquid phase. TODO: viscosity
            // of mixture
        Scalar XNaCl = fluidState.massFraction(lPhaseIdx, NaClIdx);
            Scalar result = Brine::liquidViscosity(temperature, pressure, XNaCl);
            Valgrind::CheckDefined(result);
            return result;
        }
        else if (phaseIdx == gPhaseIdx)
        {
            Scalar result = Air::gasViscosity(temperature, pressure);
            Valgrind::CheckDefined(result);
            return result;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);

    }

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
     * inverse Henry constant for the solutes and the saturated vapor pressure
     * both divided by phase pressure.
     */
    using Base::fugacityCoefficient;
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
            return vaporPressure_(T,fluidState.moleFraction(lPhaseIdx,NaClIdx))/p;
        else if (compIdx == AirIdx)
            return Dumux::BinaryCoeff::H2O_Air::henry(T)/p;
        else
            return 1/p;
        }
        else
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }


    using Base::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

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
    using Base::binaryDiffusionCoefficient;
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
            assert(compIIdx == H2OIdx);
            assert(compJIdx == AirIdx || compJIdx == NaClIdx);
            Scalar result = 0.0;
            if(compJIdx == AirIdx)
                result = Brine_Air::liquidDiffCoeff(temperature, pressure);
            else if (compJIdx == NaClIdx)
                result = 0.12e-9; //http://webserver.dmt.upm.es/~isidoro/dat1/Mass%20diffusivity%20data.htm
            else
                DUNE_THROW(Dune::NotImplemented, "Binary difussion coefficient : Incorrect compIdx");
            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            assert(phaseIdx == gPhaseIdx);

            if (compIIdx != AirIdx)
            std::swap(compIIdx, compJIdx);

            assert(compIIdx == AirIdx);
            assert(compJIdx == H2OIdx || compJIdx == NaClIdx);
            Scalar result = 0.0;
            if(compJIdx == H2OIdx)
                result = Brine_Air::gasDiffCoeff(temperature, pressure);
            else if (compJIdx == NaClIdx)
                result = 0.12e-9; //Just added to avoid numerical problem. does not have any physical significance
            else
                DUNE_THROW(Dune::NotImplemented, "Binary difussion coefficient : Incorrect compIdx");
            Valgrind::CheckDefined(result);
            return result;
        }
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     *
     * See:
     * Class 2000
     * Theorie und numerische Modellierung nichtisothermer Mehrphasenprozesse in NAPL-kontaminierten porösen Medien
     * Chapter 2.1.13 Innere Energie, Wäremekapazität, Enthalpie \cite A3:class:2001 <BR>
     *
     * Formula (2.42):
     * the specific enthalpy of a gas phase result from the sum of (enthalpies*mass fraction) of the components
     * For the calculation of enthalpy of brine we refer to (Michaelides 1981)
     *
     *  \todo This system neglects the contribution of gas-molecules in the liquid phase.
     *        This contribution is probably not big. Somebody would have to find out the enthalpy of solution for this system. ...
     */
    using Base::enthalpy;
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
            Scalar result = liquidEnthalpyBrine_(T, p, XlNaCl);
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


    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using Base::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        // TODO way too simple!
        if (phaseIdx == lPhaseIdx)
            return  0.59848; // conductivity of water[W / (m K ) ]

        else// gas phase
        return 0.0255535; // conductivity of air [W / (m K ) ]
    }

    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/(kg*K)}\f$.
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using Base::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        const Scalar temperature  = fluidState.temperature(phaseIdx);
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
        {
            // influence of air and dissolved salt is neglected here
            return H2O::liquidHeatCapacity(temperature, pressure);
        }
        else if (phaseIdx == nPhaseIdx)
        {
            //! \todo PRELIMINARY, right way to deal with solutes?
            return Air::gasHeatCapacity(temperature, pressure) * fluidState.moleFraction(nPhaseIdx, AirIdx)
                   + H2O::gasHeatCapacity(temperature, pressure) * fluidState.moleFraction(nPhaseIdx, H2OIdx);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
    /*!
     * \brief Return the molality of NaCl salt \f$\mathrm{[mol/m^3]}\f$.
     * \param fluidState An abitrary fluid state
     * \param paramCache parameter cache
     * \param salinity Salinity of brine
     */
    template <class FluidState>
    static Scalar molalityNaCl(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               Scalar salinity)
      {
        return Brine_Air::molalityNaCl(salinity);// massfraction
      }

private:
    static Scalar gasDensity_(Scalar T,
                              Scalar pg,
                              Scalar xgH2O)
    {
        Scalar pH2O = xgH2O*pg; //Dalton' Law
        Scalar pAir = pg - pH2O;
        Scalar gasDensityAir = Air::gasDensity(T, pAir);
        Scalar gasDensityH2O = H2O::gasDensity(T, pH2O);
        Scalar gasDensity = gasDensityAir + gasDensityH2O;
        return gasDensity;
    }

    /*!
     * \brief The density of pure brine at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \warning The influence of dissolved air in Brine is neglected
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * Equations given in:
     *                        - Batzle & Wang (1992) \cite batzle1992 <BR>
     *                        - cited by: Adams & Bachu in Geofluids (2002) 2, 257-271 \cite adams2002
     */

    static Scalar liquidDensity_(Scalar T,
                                 Scalar pl,
                                 Scalar xlAir,
                                 Scalar xlH2O,
                                 Scalar XlNaCl)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pl);
        Valgrind::CheckDefined(XlNaCl);
        Valgrind::CheckDefined(xlAir);

        if(T < 273.15 || T > 623.15) {
            DUNE_THROW(NumericalProblem,
                       "Liquid density for Brine and Air is only "
                       "defined between 273.15K and 623.15K (is " << T << ")");
              }
        if(pl >= 1.0e8) {
            DUNE_THROW(NumericalProblem,
                       "Liquid density for Brine and Air is only "
                       "defined below 100MPa (is " << pl << ")");
        }
        return Brine::liquidDensity(T, pl, XlNaCl); // The influence of dissolved air in Brine is neglected
    }

    static Scalar liquidEnthalpyBrine_(Scalar T,
                                       Scalar p,
                                       Scalar XlNaCl)
    {
        /* XlAir : mass fraction of Air in brine */
        /* same function as enthalpy_brine, only extended by Air content */
        /*Numerical coefficents from PALLISER (The values for F [] are given in ADAMS and BACHO)*/

        static const Scalar f[] =
        {2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10};

        /*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
        static const Scalar a[4][3] =
        {{ 9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }};

        Scalar theta, h_NaCl;
        Scalar m, h_ls, d_h, hw;
        Scalar S_lSAT, delta_h;
        int i, j;
        XlNaCl = std::abs(XlNaCl); // Added by Vishal

        theta = T - 273.15;

        S_lSAT = f[0] + f[1]*theta + f[2]*theta*theta + f[3]*theta*theta*theta; // This is NaCl specific (S_lSAT = 0.265234)
       // std::cout<<"saturation limit : " << S_lSAT<< std::endl;

        /*Regularization*/
        if (XlNaCl > S_lSAT) {
            XlNaCl = S_lSAT;
        }

        hw = H2O::liquidEnthalpy(T, p) /1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T // enthalpy of Halite (Confirm with Dr. Melani)
                        +((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        m = (1E3/58.44)*(XlNaCl/(1-XlNaCl));
        i = 0;
        j = 0;
        d_h = 0;

        for (i = 0; i<=3; i++) {
            for (j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * pow(theta, i) * pow(m, j);
            }
        }
        /* heat of dissolution for halite according to Michaelides 1971 */
        delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without Air */
        h_ls = (1-XlNaCl)*hw + XlNaCl*h_NaCl + XlNaCl*delta_h; /* kJ/kg */
        return (h_ls);
    }

    static Scalar vaporPressure_(Scalar T, Scalar x)
    {
     Scalar p_0 = H2O::vaporPressure(T);//Saturation vapor pressure for pure water
     Scalar p_s = p_0; // modified saturation vapor pressure for saline water
// #if SALINIZATION
//     Scalar vw = 18.0e-6;//[m3/mol] volume per unit mole of water
//     Scalar R = 8.314;//[j/K mol] universal gas constant
//     Scalar pi = (R * T * std::log(1- x))/vw;
//     if (x > 0.26) // here we have hard coaded the solubility limit for NaCl
//      pi = (R * T * std::log(0.74))/vw;
//     p_s = p_0 * std::exp((pi*vw)/(R*T));// Kelvin's law for reduction in saturation vapor pressure due to osmotic potential
// #endif
      return p_s;
    }
};

} // end namespace
} // end namespace

#endif
