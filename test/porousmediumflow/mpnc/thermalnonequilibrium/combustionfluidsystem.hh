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
 * \ingroup MPNCTests
 * \brief @copybrief Dumux::FluidSystems::CombustionFluidsystem
 */
#ifndef DUMUX_PURE_WATER_FLUID_SYSTEM_HH
#define DUMUX_PURE_WATER_FLUID_SYSTEM_HH

#include <cassert>

#include <dumux/material/idealgas.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>

#include <dumux/common/exceptions.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup MPNCTests
 *
 * \brief A two-phase fluid system with water as sole component.
 *
 * Values are taken from Shi and Wang, 2011 \cite Shi2011.
 */
template <class Scalar>
class CombustionFluidsystem
    : public Base<Scalar, CombustionFluidsystem<Scalar> >
{
    using ThisType = CombustionFluidsystem<Scalar>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

    // convenience using declarations
    using IdealGas = Dumux::IdealGas<Scalar>;
    using SimpleH2O = Dumux::Components::SimpleH2O<Scalar>;
    using SimpleN2 = Dumux::Components::N2<Scalar>;

public:
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 2;

    static constexpr int wPhaseIdx = 0; // index of the wetting phase
    static constexpr int nPhaseIdx = 1; // index of the nonwetting phase
    static constexpr int phase0Idx = 0; // index of the wetting phase
    static constexpr int phase1Idx = 1; // index of the nonwetting phase

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static constexpr int wCompIdx = wPhaseIdx;
    static constexpr int nCompIdx = nPhaseIdx;
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
            case wPhaseIdx: return "liq";
            case nPhaseIdx: return "gas";
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
        return phaseIdx == nPhaseIdx;
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
        // gases are always compressible
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

        if (phaseIdx == nPhaseIdx)
            // let the components decide
            return H2O::gasIsIdeal() && N2::gasIsIdeal();
        return false; // not a gas
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 2;

    static constexpr int H2OIdx = wCompIdx;
    static constexpr int N2Idx = nCompIdx;

    //! The components for pure water
    using H2O = SimpleH2O;

    //! The components for pure nitrogen
    using N2 = SimpleN2;

    /*!
     * \brief Returns the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        static std::string name[] = {
            H2O::name(),
            N2::name()
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
            H2O::molarMass(),
            N2::molarMass(),
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
            H2O::criticalTemperature(), // H2O
            N2::criticalTemperature() // N2
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
            H2O::criticalPressure(),
            N2::criticalPressure()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    }

    /*!
     * \brief Molar volume of a component at the critical point \f$\mathrm{[m^3/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "criticalMolarVolume()");
    }

    /*!
     * \brief The acentric factor of a component \f$\mathrm{[-]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            H2O::acentricFactor(), // H2O (from Reid, et al.)
            N2::acentricFactor()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return accFac[compIdx];
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initializes the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
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
        std::cout << "Using very simple pure water fluid system\n";
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
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        // liquid phase
        if (phaseIdx == wPhaseIdx)
        {
            return 1044.0;
        }
        else if (phaseIdx == nPhaseIdx)// gas phase
        {
            return 1.679;
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
        if (phaseIdx == wPhaseIdx)
        {
            return density(fluidState, phaseIdx)/fluidState.averageMolarMass(phaseIdx);
        }
        else if (phaseIdx == nPhaseIdx)
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
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        // liquid phase
        if (phaseIdx == wPhaseIdx)
        {
            return 2.694e-7 * density(fluidState, phaseIdx);
        }
        else if (phaseIdx == nPhaseIdx) // gas phase
        {
            return 7.16e-6 * density(fluidState, phaseIdx);
        }
        else DUNE_THROW(Dune::NotImplemented, "Wrong phase index");
    }

    /*!
     * \brief Calculates the temperature \f$\mathrm{[K]}\f$ of vapor at a given pressure on the vapor pressure curve.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar vaporTemperature(const FluidState &fluidState,
                                   const unsigned int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        Scalar pressure = fluidState.pressure(nPhaseIdx);

        return IAPWS::Region4<Scalar>::vaporTemperature( pressure );
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
        if (phaseIdx == wPhaseIdx)
        {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            return BinaryCoeff::H2O_N2::henry(T)/p;
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
        DUNE_THROW(Dune::NotImplemented, "Binary Diffusion coefficients");
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
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        Scalar temperature = fluidState.temperature(phaseIdx);

        const Scalar cp = heatCapacity(fluidState, phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx)
        {
            return cp * (temperature - 373.15);
        }
        else if (phaseIdx == nPhaseIdx) // gas phase
        {
            return cp * (temperature - 373.15) + 2.257e6;
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
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        // liquid phase
        if (phaseIdx == wPhaseIdx)
        {
            return 0.68;
        }
        else if (phaseIdx == nPhaseIdx) // gas phase
        {
            return 0.0248;
        }
        else DUNE_THROW(Dune::NotImplemented, "Wrong phase index");
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
        if (phaseIdx == wPhaseIdx)
        {
            return 4.217e3;
        }
        else if (phaseIdx == nPhaseIdx) // gas phase
        {
            return 2.029e3;
        }
        else DUNE_THROW(Dune::NotImplemented, "Wrong phase index");
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
