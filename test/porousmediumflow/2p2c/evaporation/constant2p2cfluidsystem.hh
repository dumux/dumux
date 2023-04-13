// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief @copybrief Dumux::FluidSystems::ConstantTwoPTwoCFluidsystem
 */
#ifndef DUMUX_CONSTANT_2P2C_FLUID_SYSTEM_HH
#define DUMUX_CONSTANT_2P2C_FLUID_SYSTEM_HH

#include <cassert>

#include <dumux/material/idealgas.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/binarycoefficients/h2o_constant.hh>

#include <dumux/common/exceptions.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup TwoPTwoCTests
 *
 * \brief A two-phase fluid system with two components with constant properties that can be specified via the input file
 */
template <class Scalar, class Component0, class Component1>
class ConstantTwoPTwoCFluidsystem
    : public Base<Scalar, ConstantTwoPTwoCFluidsystem<Scalar, Component0, Component1> >
{
    using ThisType = ConstantTwoPTwoCFluidsystem<Scalar, Component0, Component1>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

public:

    using ComponentOne = Component0;
    using ComponentTwo = Component1;
    //! Number of phases in the fluid system
    static constexpr int numPhases = 2;

    static constexpr int numComponents = 2; //!< Number of components in the fluid system

    static constexpr int phase0Idx = 0; // index of the wetting phase
    static constexpr int phase1Idx = 1; // index of the nonwetting phase

    // export component indices to indicate the main component
    static constexpr int comp0Idx = 0; // index of the wetting phase
    static constexpr int comp1Idx = 1; // index of the nonwetting phase

    /*!
     * \brief Returns the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        switch (phaseIdx)
        {
            case phase0Idx: return "liq";
            case phase1Idx: return "gas";
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns whether a phase is gaseous
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx == phase1Idx;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
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
        // gases are always compressible
        if (phaseIdx == phase1Idx)
            return true;
        // the liquid component decides for the liquid phase...
        return ComponentOne::liquidIsCompressible();
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

        if (phaseIdx == phase1Idx)
            // let the components decide
            return ComponentOne::gasIsIdeal() && ComponentTwo::gasIsIdeal();
        return false; // not a gas
    }

    /*!
     * \brief Returns the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        static std::string name[] = {
            ComponentOne::name(),
            ComponentTwo::name()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \brief Returns the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        static const Scalar M[] = {
            ComponentOne::molarMass(),
            ComponentTwo::molarMass(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return M[compIdx];
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        static const Scalar Tcrit[] = {
            ComponentOne::criticalTemperature(), // ComponentOne
            ComponentTwo::criticalTemperature() // ComponentTwo
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return Tcrit[compIdx];
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        static const Scalar pcrit[] = {
            ComponentOne::criticalPressure(),
            ComponentTwo::criticalPressure()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    }


    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initializes the fluid system's static parameters generically
     *
     * If a tabulated ComponentOne component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initializes the fluid system's static parameters using
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
        std::cout << "Using very simple constant component fluid system\n";
    }

    using Base::density;
    /*!
     * \brief Calculates the density \f$\mathrm{[kg/m^3]}\f$ of a fluid phase
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == phase0Idx)
        {
            return ComponentOne::liquidDensity(T, p);
        }
        else if (phaseIdx == phase1Idx)// gas phase
        {
            return ComponentTwo::gasDensity(T, p);
        }
        else DUNE_THROW(Dune::NotImplemented, "Wrong phase index");
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * This is a specific molar density  for the combustion test
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        if (phaseIdx == phase0Idx)
        {
            return density(fluidState, phaseIdx)/fluidState.averageMolarMass(phaseIdx);
        }
        else if (phaseIdx == phase1Idx)
        {
            return density(fluidState, phaseIdx)/fluidState.averageMolarMass(phaseIdx);
        }
        else DUNE_THROW(Dune::NotImplemented, "Wrong phase index");
    }

    using Base::viscosity;
    /*!
     * \brief Calculates the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == phase0Idx)
        {
            return ComponentOne::liquidViscosity(T, p);
        }
        else if (phaseIdx == phase1Idx) // gas phase
        {
            return ComponentTwo::gasViscosity(T, p);
        }
        else DUNE_THROW(Dune::NotImplemented, "Wrong phase index");
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Calculates the fugacity coefficient \f$\mathrm{[-]}\f$ of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\phi^\kappa_\alpha\f$ of
     * component \f$\kappa\f$ in phase \f$\alpha\f$ is connected to
     * the fugacity \f$f^\kappa_\alpha\f$ and the component's mole
     * fraction \f$x^\kappa_\alpha\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     \f]
     * where \f$p_\alpha\f$ is the pressure of the fluid phase.
     *
     * The quantity "fugacity" itself is just an other way to express
     * the chemical potential \f$\zeta^\kappa_\alpha\f$ of the
     * component. It is defined via
     *
     * \f[
     f^\kappa_\alpha := \exp\left\{\frac{\zeta^\kappa_\alpha}{k_B T_\alpha} \right\}
     \f]
     * where \f$k_B = 1.380\cdot10^{-23}\;J/K\f$ is the Boltzmann constant.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        // liquid phase
        if (phaseIdx == phase0Idx)
        {
            if (compIdx == comp0Idx)
                return ComponentOne::vaporPressure(T)/p;
            return BinaryCoeff::H2O_Component<Scalar, ComponentTwo>::henryCompInWater(T)/p;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    using Base::diffusionCoefficient;
    /*!
     * \brief Calculates the molecular diffusion coefficient for a
     *        component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
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
     * \brief Given a phase's composition, temperature and pressure,
     *        returns the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$i\f$ and \f$j\f$ in this phase.
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
        using std::swap;
        if (compIIdx > compJIdx)
            swap(compIIdx, compJIdx);

        // we are in the liquid phase
        if (phaseIdx == phase0Idx)
        {
            if (compIIdx == comp0Idx && compJIdx == comp1Idx)
                return 1e-9;
            else
                DUNE_THROW(Dune::InvalidStateException,
                           "Binary diffusion coefficient of components "
                            << compIIdx << " and " << compJIdx
                            << " in phase " << phaseIdx << " is undefined!\n");
        }

        // we are in the gas phase
        else if (phaseIdx == phase1Idx)
        {
            if (compIIdx == comp0Idx && compJIdx == comp1Idx)
                return 1e-5;
            else
                DUNE_THROW(Dune::InvalidStateException,
                           "Binary diffusion coefficient of components "
                           << compIIdx << " and " << compJIdx
                           << " in phase " << phaseIdx << " is undefined!\n");
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::enthalpy;
    /*!
     * \brief Calculates specific enthalpy \f$\mathrm{[J/kg]}\f$.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == phase0Idx)
        {
             return Component0::liquidEnthalpy(T, p);
        }
        else if (phaseIdx == phase1Idx) // gas phase
        {
            // assume ideal mixture: which means
            // that the total specific enthalpy is the sum of the
            // "partial specific enthalpies" of the components.
            Scalar hComponentOne =
                fluidState.massFraction(phase1Idx, comp0Idx)
                * ComponentOne::gasEnthalpy(T, p);
            Scalar hComponentTwo =
                fluidState.massFraction(phase1Idx, comp1Idx)
                * ComponentTwo::gasEnthalpy(T, p);
            return hComponentOne + hComponentTwo;
        }
        else DUNE_THROW(Dune::NotImplemented, "Wrong phase index");
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     *
     * Use the conductivity of vapor and liquid water at 100°C
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        const Scalar temperature  = fluidState.temperature(phaseIdx) ;
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == phase0Idx)
            return ComponentOne::liquidThermalConductivity(temperature, pressure);
        else
            return ComponentTwo::gasThermalConductivity(temperature, pressure);
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/kg / K]}\f$.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        // liquid phase
        const Scalar temperature  = fluidState.temperature(phaseIdx);
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == phase0Idx)
        {
            return ComponentOne::liquidHeatCapacity(temperature, pressure);
        }
        else if (phaseIdx == phase1Idx) // gas phase
        {
            return ComponentOne::gasHeatCapacity(temperature, pressure) * fluidState.moleFraction(phase1Idx, comp0Idx)
                   + ComponentTwo::gasHeatCapacity(temperature, pressure) * fluidState.moleFraction(phase1Idx, comp1Idx);
        }
        else DUNE_THROW(Dune::NotImplemented, "Wrong phase index");
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
