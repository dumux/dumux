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
/*!
 * \file
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::FluidSystems::BrineAir
 */
#ifndef DUMUX_BRINE_AIR_FLUID_SYSTEM_HH
#define DUMUX_BRINE_AIR_FLUID_SYSTEM_HH

#include <array>
#include <cassert>
#include <iomanip>

#include <dumux/material/idealgas.hh>
#include <dumux/material/fluidstates/adapter.hh>
#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/fluidsystems/brine.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/nacl.hh>
#include <dumux/material/binarycoefficients/h2o_air.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/common/exceptions.hh>

#include <dumux/io/name.hh>

#include "brine.hh"

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief Policy for the brine-air fluid system
 */
template<bool fastButSimplifiedRelations = false>
struct BrineAirDefaultPolicy
{
    static constexpr bool useBrineDensityAsLiquidMixtureDensity() { return fastButSimplifiedRelations;}
    static constexpr bool useIdealGasDensity() { return fastButSimplifiedRelations; }
};

/*!
 * \ingroup Fluidsystems
 * \brief A compositional two-phase fluid system with a liquid and a gaseous phase
 *        and \f$H_2O\f$, \f$Air\f$ and \f$S\f$ (dissolved minerals) as components.
 *
 * \note This fluidsystem is applied by default with the tabulated version of
 *       water of the IAPWS-formulation.
 */
template <class Scalar,
          class H2Otype = Components::TabulatedComponent<Components::H2O<Scalar>>,
          class Policy = BrineAirDefaultPolicy<>>
class BrineAir
: public Base<Scalar, BrineAir<Scalar, H2Otype, Policy>>
{
    using ThisType = BrineAir<Scalar, H2Otype, Policy>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    //! export the involved components
    using H2O = H2Otype;
    using Air = Components::Air<Scalar>;
    using NaCl = Components::NaCl<Scalar>;

    //! export the underlying brine fluid system for the liquid phase
    using Brine = Dumux::FluidSystems::Brine<Scalar, H2Otype>;

    //! export the binary coefficients between air and water
    using H2O_Air = BinaryCoeff::H2O_Air;

    //! the type of parameter cache objects
    using ParameterCache = NullParameterCache;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    static constexpr int numPhases = 2;     // one liquid and one gas phase
    static constexpr int numComponents = 3; // H2O, Air, NaCl

    static constexpr int liquidPhaseIdx = 0; // index of the liquid phase
    static constexpr int gasPhaseIdx = 1;    // index of the gas phase

    static constexpr int phase0Idx = liquidPhaseIdx; // index of the first phase
    static constexpr int phase1Idx = gasPhaseIdx;    // index of the second phase

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static constexpr int H2OIdx = 0;
    static constexpr int AirIdx = 1;
    static constexpr int NaClIdx = 2;
    static constexpr int comp0Idx = H2OIdx;
    static constexpr int comp1Idx = AirIdx;
    static constexpr int comp2Idx = NaClIdx;

private:
    struct BrineAdapterPolicy
    {
        using FluidSystem = Brine;

        static constexpr int phaseIdx(int brinePhaseIdx) { return liquidPhaseIdx; }
        static constexpr int compIdx(int brineCompIdx)
        {
            switch (brineCompIdx)
            {
                case Brine::H2OIdx: return H2OIdx;
                case Brine::NaClIdx: return NaClIdx;
                default: return 0; // this will never be reached, only needed to suppress compiler warning
            }
        }
    };

    template<class FluidState>
    using BrineAdapter = FluidStateAdapter<FluidState, BrineAdapterPolicy>;

public:

    /****************************************
     * phase related static parameters
     ****************************************/

    /*!
     * \brief Return the human readable name of a fluid phase
     * \param phaseIdx index of the phase
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        switch (phaseIdx)
        {
            case liquidPhaseIdx: return IOName::liquidPhase();
            case gasPhaseIdx: return IOName::gaseousPhase();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns whether the fluids are miscible
     */
    static constexpr bool isMiscible()
    { return true; }

    /*!
     * \brief Return whether a phase is gaseous
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx == gasPhaseIdx;
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
        if (phaseIdx == gasPhaseIdx)
            return true;
        // let brine decide for the liquid phase...
        return Brine::isCompressible(Brine::liquidPhaseIdx);
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // let the fluids decide
        if (phaseIdx == gasPhaseIdx)
            return H2O::gasIsIdeal() && Air::gasIsIdeal();
        return false; // not a gas
    }

    /*!
     * \brief Get the main component of a given phase if possible
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr int getMainComponent(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        if (phaseIdx == liquidPhaseIdx)
            return H2OIdx;
        else
            return AirIdx;
    }

    /****************************************
     * Component related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a component
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        switch (compIdx)
        {
            case H2OIdx: return H2O::name();
            case AirIdx: return Air::name();
            case NaClIdx: return NaCl::name();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        switch (compIdx)
        {
            case H2OIdx: return H2O::molarMass();
            case AirIdx: return Air::molarMass();
            case NaClIdx: return NaCl::molarMass();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Vapor pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param fluidState The fluid state
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar vaporPressure(const FluidState& fluidState, int compIdx)
    {
        // The vapor pressure of the water is affected by the
        // salinity, thus, we forward to the interface of Brine here
        if (compIdx == H2OIdx)
            return Brine::vaporPressure(BrineAdapter<FluidState>(fluidState), Brine::H2OIdx);
        else if (compIdx == NaClIdx)
            DUNE_THROW(Dune::NotImplemented, "NaCl::vaporPressure(t)");
        else
            DUNE_THROW(Dune::NotImplemented, "Invalid component index " << compIdx);
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
        std::cout << "The brine-air fluid system was configured with the following policy:\n";
        std::cout << " - use brine density as liquid mixture density: " << std::boolalpha << Policy::useBrineDensityAsLiquidMixtureDensity() << "\n";
        std::cout << " - use ideal gas density: " << std::boolalpha << Policy::useIdealGasDensity() << std::endl;

        if (H2O::isTabulated)
            H2O::init(tempMin, tempMax, nTemp, pressMin, pressMax, nPress);
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
    static Scalar density(const FluidState &fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto T = fluidState.temperature(phaseIdx);
        const auto p = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx)
        {
            // assume pure brine
            if (Policy::useBrineDensityAsLiquidMixtureDensity())
                return Brine::density(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);

            // assume one molecule of gas replaces one "brine" molecule
            else
                return Brine::molarDensity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx)
                       *(H2O::molarMass()*fluidState.moleFraction(liquidPhaseIdx, H2OIdx)
                         + NaCl::molarMass()*fluidState.moleFraction(liquidPhaseIdx, NaClIdx)
                         + Air::molarMass()*fluidState.moleFraction(liquidPhaseIdx, AirIdx));
        }
        else if (phaseIdx == phase1Idx)
        {
            // for the gas phase assume an ideal gas
            if (Policy::useIdealGasDensity())
                return IdealGas::density(fluidState.averageMolarMass(phase1Idx), T, p);

            // if useComplexRelations = true, compute density. NaCl is assumed
            // not to be present in gas phase, NaCl has only solid interfaces implemented
            return H2O::gasDensity(T, fluidState.partialPressure(phase1Idx, H2OIdx))
                   + Air::gasDensity(T, fluidState.partialPressure(phase1Idx, AirIdx));
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density for the complex relation is defined by the
     * mass density \f$\rho_\alpha\f$ and the mean molar mass \f$\overline M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{\overline M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return Brine::molarDensity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else if (phaseIdx == phase1Idx)
        {
            const Scalar T = fluidState.temperature(phaseIdx);

            // for the gas phase assume an ideal gas
            if (Policy::useIdealGasDensity())
                return IdealGas::molarDensity(T, fluidState.pressure(phaseIdx));

            // if useComplexRelations = true, compute density. NaCl is assumed
            // not to be present in gas phase, NaCl has only solid interfaces implemented
            return H2O::gasMolarDensity(T, fluidState.partialPressure(phase1Idx, H2OIdx))
                   + Air::gasMolarDensity(T, fluidState.partialPressure(phase1Idx, AirIdx));
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == liquidPhaseIdx)
            return Brine::viscosity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else
            return Air::gasViscosity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
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
    static Scalar fugacityCoefficient(const FluidState& fluidState, int phaseIdx, int compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == gasPhaseIdx)
            return 1.0;

        else if (phaseIdx == liquidPhaseIdx)
        {
            // TODO: Should we use the vapor pressure of the mixture (brine) here?
            //       The presence of NaCl lowers the vapor pressure, thus, we would
            //       expect the fugacity coefficient to be lower as well. However,
            //       with the fugacity coefficient being dependent on the salinity,
            //       the equation system for the phase equilibria becomes non-linear
            //       and our constraint solvers assume linear system of equations.
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;

            else if (compIdx == AirIdx)
                return BinaryCoeff::H2O_Air::henry(T)/p;

            // we assume nacl always stays in the liquid phase
            else
                return 0.0;
        }

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
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIIdx && compIIdx < numComponents);
        assert(0 <= compJIdx && compJIdx < numComponents);

        const auto T = fluidState.temperature(phaseIdx);
        const auto p = fluidState.pressure(phaseIdx);

        if (compIIdx > compJIdx)
            std::swap(compIIdx, compJIdx);

        if (phaseIdx == liquidPhaseIdx)
        {
            if(compIIdx == H2OIdx && compJIdx == AirIdx)
                return H2O_Air::liquidDiffCoeff(T, p);
            else if (compIIdx == H2OIdx && compJIdx == NaClIdx)
                return Brine::binaryDiffusionCoefficient(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx, Brine::H2OIdx, Brine::NaClIdx);
            else
                DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components "
                                                 << compIIdx << " and " << compJIdx
                                                 << " in phase " << phaseIdx);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            if (compIIdx == H2OIdx && compJIdx == AirIdx)
                return H2O_Air::gasDiffCoeff(T, p);

            // NaCl is expected to never be present in the gas phase. we need to
            // return a diffusion coefficient that does not case numerical problems.
            // We choose a very small value here.
            else if (compIIdx == AirIdx && compJIdx == NaClIdx)
                return 1e-12;

            else
                DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components "
                                                 << compIIdx << " and " << compJIdx
                                                 << " in phase " << phaseIdx);
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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
    static Scalar enthalpy(const FluidState& fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx)
            return Brine::enthalpy(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else
        {
            // This assumes NaCl not to be present in the gas phase
            const Scalar XAir = fluidState.massFraction(gasPhaseIdx, AirIdx);
            const Scalar XH2O = fluidState.massFraction(gasPhaseIdx, H2OIdx);
            return XH2O*H2O::gasEnthalpy(T, p) + XAir*Air::gasEnthalpy(T, p);
        }
    }

   /*!
    * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in a specific phase
    * \param fluidState The fluid state
    * \param phaseIdx The index of the phase
    * \param componentIdx The index of the component
    */
    template <class FluidState>
    static Scalar componentEnthalpy(const FluidState& fluidState, int phaseIdx, int componentIdx)
    {
        const Scalar T = fluidState.temperature(gasPhaseIdx);
        const Scalar p = fluidState.pressure(gasPhaseIdx);

        if (phaseIdx == liquidPhaseIdx)
            DUNE_THROW(Dune::NotImplemented, "The component enthalpies in the liquid phase are not implemented.");

        else if (phaseIdx == gasPhaseIdx)
        {
            if (componentIdx == H2OIdx)
                return H2O::gasEnthalpy(T, p);
            else if (componentIdx == AirIdx)
                return Air::gasEnthalpy(T, p);
            else if (componentIdx == NaClIdx)
                DUNE_THROW(Dune::InvalidStateException, "Implementation assumes NaCl not to be present in gas phase");
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
     * \todo TODO: For the thermal conductivity of the air phase the contribution of the minor
     *       components is neglected. This contribution is probably not big, but somebody
     *       would have to find out its influence.
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return Brine::thermalConductivity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else if (phaseIdx == gasPhaseIdx)
            return Air::gasThermalConductivity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/(kg*K)}\f$.
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * \todo TODO: The calculation of the isobaric heat capacity is preliminary. A better
     *       description of the influence of the composition on the air phase property
     *       has to be found.
     */
    using Base::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState, int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx)
            return Brine::heatCapacity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);

        // We assume NaCl not to be present in the gas phase here
        else if (phaseIdx == gasPhaseIdx)
            return Air::gasHeatCapacity(T, p)*fluidState.moleFraction(gasPhaseIdx, AirIdx)
                   + H2O::gasHeatCapacity(T, p)*fluidState.moleFraction(gasPhaseIdx, H2OIdx);

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
