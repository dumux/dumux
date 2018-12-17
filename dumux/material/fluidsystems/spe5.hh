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
 * \brief @copybrief Dumux::FluidSystems::Spe5
 */
#ifndef DUMUX_SPE5_FLUID_SYSTEM_HH
#define DUMUX_SPE5_FLUID_SYSTEM_HH

#include "spe5parametercache.hh"

#include <dumux/common/spline.hh>
#include <dumux/material/eos/pengrobinsonmixture.hh>
#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief The fluid system for the SPE-5 benchmark problem.
 *
 * A three-phase fluid system featuring gas, oil and water as phases and the seven
 * components distilled water, Methane \f$(\mathrm{C_1})\f$, Propane \f$(\mathrm{C_3})\f$,
 * Pentane \f$(\mathrm{C_5})\f$, Heptane \f$(\mathrm{C_7})\f$, Decane
 * \f$(\mathrm{C_{10}})\f$, Pentadecane \f$(\mathrm{C_{15}})\f$ and Icosane
 * \f$(\mathrm{C_{20}})\f$. For the water phase the IAPWS-97 formulation is used as
 * equation of state, while for the gas and oil phases a Peng-Robinson
 * equation of state with slightly modified parameters is used. This fluid system is highly
 * non-linear, and the gas and oil phases also cannot be considered ideal.
 *
 * See:
 *
 * J.E. Killough, et al.: Fifth Comparative Solution Project:
 * Evaluation of Miscible Flood Simulators, Ninth SPE Symposium on
 * Reservoir Simulation, 1987 \cite SPE5
 */
template <class Scalar>
class Spe5
{
    using ThisType = FluidSystems::Spe5<Scalar>;

    using PengRobinsonMixture = Dumux::PengRobinsonMixture<Scalar, ThisType>;
    using PengRobinson = Dumux::PengRobinson<Scalar>;

public:
    using ParameterCache = Spe5ParameterCache<Scalar, ThisType>;

    /****************************************
     * Fluid phase parameters
     ****************************************/

    //! Number of phases in the fluid system
    static const int numPhases = 3;

    //! Index of the gas phase
    static const int gPhaseIdx = 0;
    //! Index of the water phase
    static const int wPhaseIdx = 1;
    //! Index of the oil phase
    static const int oPhaseIdx = 2;

    //! The component for pure water to be used
    using H2O = Dumux::Components::H2O<Scalar>;

    /*!
     * \brief Return the human readable name of a fluid phase
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        switch (phaseIdx)
        {
            case gPhaseIdx: return IOName::gaseousPhase();
            case wPhaseIdx: return IOName::aqueousPhase();
            case oPhaseIdx: return IOName::naplPhase();
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
        return phaseIdx == gPhaseIdx;
    }

    /*!
     * \brief Return whether a phase is compressible
     * \param phaseIdx The index of the fluid phase to consider
     *
     * In the SPE-5 problems all fluids are compressible...
     */
    static constexpr bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return true;
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
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx)
    {
        // always use the reference oil for the fugacity coefficients,
        // so they cannot be dependent on composition and they the
        // phases thus always an ideal mixture
        return phaseIdx == wPhaseIdx;
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
        return false; // gas is not ideal here!
    }

    /****************************************
     * Component related parameters
     ****************************************/

    //! Number of components in the fluid system
    static const int numComponents = 7;

    static const int H2OIdx = 0;
    static const int C1Idx = 1;
    static const int C3Idx = 2;
    static const int C6Idx = 3;
    static const int C10Idx = 4;
    static const int C15Idx = 5;
    static const int C20Idx = 6;

    /*!
     * \brief Return the human readable name of a component
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        static std::string name[] = {
            H2O::name(),
            std::string("C1"),
            std::string("C3"),
            std::string("C6"),
            std::string("C10"),
            std::string("C15"),
            std::string("C20")
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        static const Scalar M[] = {
            H2O::molarMass(),
            16.04e-3, // C1
            44.10e-3, // C3
            86.18e-3, // C6
            142.29e-3, // C10
            206.00e-3, // C15
            282.00e-3 // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return M[compIdx];
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        static const Scalar Tcrit[] = {
            H2O::criticalTemperature(), // H2O
            343.0*5/9, // C1
            665.7*5/9, // C3
            913.4*5/9, // C6
            1111.8*5/9, // C10
            1270.0*5/9, // C15
            1380.0*5/9 // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return Tcrit[compIdx];
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        static const Scalar pcrit[] = {
            H2O::criticalPressure(), // H2O
            667.8 * 6894.7573, // C1
            616.3 * 6894.7573, // C3
            436.9 * 6894.7573, // C6
            304.0 * 6894.7573, // C10
            200.0 * 6894.7573, // C15
            162.0 * 6894.7573  // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    }

    /*!
     * \brief Molar volume of a component at the critical point \f$\mathrm{[m^3/mol]}\f$.
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        static const Scalar R = 8.314472;
        static const Scalar vcrit[] = {
            H2O::criticalMolarVolume(), // H2O
            0.290*R*criticalTemperature(1)/criticalPressure(1), // C1
            0.277*R*criticalTemperature(2)/criticalPressure(2), // C3
            0.264*R*criticalTemperature(3)/criticalPressure(3), // C6
            0.257*R*criticalTemperature(4)/criticalPressure(4), // C10
            0.245*R*criticalTemperature(5)/criticalPressure(5), // C15
            0.235*R*criticalTemperature(6)/criticalPressure(6) // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return vcrit[compIdx];
    }

    /*!
     * \brief The acentric factor of a component \f$\mathrm{[-]}\f$.
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            0.344, // H2O (from Reid, et al.)
            0.0130, // C1
            0.1524, // C3
            0.3007, // C6
            0.4885, // C10
            0.6500, // C15
            0.8500  // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return accFac[compIdx];
    }

    /*!
     * \brief Returns the interaction coefficient for two components.
     * \param comp1Idx The index of component 1 to consider
     * \param comp2Idx The index of component 2 to consider
     *
     * The values are given by the SPE5 paper.
     */
    static Scalar interactionCoefficient(int comp1Idx, int comp2Idx)
    {
        using std::min;
        using std::max;
        int i = min(comp1Idx, comp2Idx);
        int j = max(comp1Idx, comp2Idx);
        if (i == C1Idx && (j == C15Idx || j == C20Idx))
            return 0.05;
        else if (i == C3Idx && (j == C15Idx || j == C20Idx))
            return 0.005;
        return 0;
    }

    /****************************************
     * Methods which compute stuff
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    {
        PengRobinsonParamsMixture<Scalar, ThisType, gPhaseIdx, /*useSpe5=*/true> prParams;

        // find envelopes of the 'a' and 'b' parameters for the range
        // 273.15K <= T <= 373.15K and 10e3 Pa <= p <= 100e6 Pa. For
        // this we take advantage of the fact that 'a' and 'b' for
        // mixtures is just a convex combination of the attractive and
        // repulsive parameters of the pure components
        Scalar minT = 273.15, maxT = 373.15;
        Scalar minP = 1e4, maxP = 100e6;

        Scalar minA = 1e100, maxA = -1e100;
        Scalar minB = 1e100, maxB = -1e100;

        prParams.updatePure(minT, minP);
        using std::min;
        using std::max;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = min(prParams.pureParams(compIdx).a(), minA);
            maxA = max(prParams.pureParams(compIdx).a(), maxA);
            minB = min(prParams.pureParams(compIdx).b(), minB);
            maxB = max(prParams.pureParams(compIdx).b(), maxB);
        }

        prParams.updatePure(maxT, minP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = min(prParams.pureParams(compIdx).a(), minA);
            maxA = max(prParams.pureParams(compIdx).a(), maxA);
            minB = min(prParams.pureParams(compIdx).b(), minB);
            maxB = max(prParams.pureParams(compIdx).b(), maxB);
        }

        prParams.updatePure(minT, maxP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = min(prParams.pureParams(compIdx).a(), minA);
            maxA = max(prParams.pureParams(compIdx).a(), maxA);
            minB = min(prParams.pureParams(compIdx).b(), minB);
            maxB = max(prParams.pureParams(compIdx).b(), maxB);
        }

        prParams.updatePure(maxT, maxP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = min(prParams.pureParams(compIdx).a(), minA);
            maxA = max(prParams.pureParams(compIdx).a(), maxA);
            minB = min(prParams.pureParams(compIdx).b(), minB);
            maxB = max(prParams.pureParams(compIdx).b(), maxB);
        }

        PengRobinson::init(/*aMin=*/minA, /*aMax=*/maxA, /*na=*/100,
                           /*bMin=*/minB, /*bMax=*/maxB, /*nb=*/200);
    }

    /*!
     * \brief Calculate the density \f$\mathrm{[kg/m^3]}\f$ of a fluid phase
     *
     *
     * \param fluidState An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        return fluidState.averageMolarMass(phaseIdx)/paramCache.molarVolume(phaseIdx);
    }

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
    static Scalar molarDensity(const FluidState &fluidState, const ParameterCache &paramCache, int phaseIdx)
    {
        return 1.0/paramCache.molarVolume(phaseIdx);
    }

    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     * \param fs An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fs,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        if (phaseIdx == gPhaseIdx) {
            // given by SPE-5 in table on page 64. we use a constant
            // viscosity, though...
            return 0.0170e-2 * 0.1;
        }
        else if (phaseIdx == wPhaseIdx)
            // given by SPE-5: 0.7 centi-Poise  = 0.0007 Pa s
            return 0.7e-2 * 0.1;
        else {
            assert(phaseIdx == oPhaseIdx);
            // given by SPE-5 in table on page 64. we use a constant
            // viscosity, though...
            return 0.208e-2 * 0.1;
        }
    }

    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[-]}\f$ of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\mathrm{\phi^\kappa_\alpha}\f$ is connected
     * to the fugacity \f$\mathrm{f^\kappa_\alpha}\f$ and the component's mole
     * fraction in a phase \f$\mathrm{x^\kappa_\alpha}\f$ by means of the
     * relation
     *
     * \f[ f^\kappa_\alpha = \phi^\kappa_\alpha \cdot x^\kappa_\alpha \cdot p_alpha \f]
     * \param fs An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fs,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        assert(0 <= compIdx  && compIdx <= numComponents);

        if (phaseIdx == oPhaseIdx || phaseIdx == gPhaseIdx)
            return PengRobinsonMixture::computeFugacityCoefficient(fs,
                                                                   paramCache,
                                                                   phaseIdx,
                                                                   compIdx);
        else {
            assert(phaseIdx == wPhaseIdx);
            return henryCoeffWater_(compIdx, fs.temperature(wPhaseIdx))
                   / fs.pressure(wPhaseIdx);
        }
    }


    /*!
     * \brief Calculate the binary molecular diffusion coefficient for
     *        a component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     *
     * Molecular diffusion of a compoent \f$\mathrm{\kappa}\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D grad \mu_\kappa \f]
     *
     * where \f$\mathrm{\mu_\kappa}\f$ is the component's chemical potential,
     * \f$\mathrm{D}\f$ is the diffusion coefficient and \f$\mathrm{J}\f$ is the
     * diffusive flux. \f$\mathrm{\mu_\kappa}\f$ is connected to the component's
     * fugacity \f$\mathrm{f_\kappa}\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$\mathrm{p_\alpha}\f$ and \f$\mathrm{T_\alpha}\f$ are the fluid phase'
     * pressure and temperature.
     * \param fs An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fs,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    { DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients"); }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     * \param fluidState An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             const ParameterCache &paramCache,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    { DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficients"); }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        calculate its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fs An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fs,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    { DUNE_THROW(Dune::NotImplemented, "Enthalpies"); }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        calculate its thermal conductivity \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    { DUNE_THROW(Dune::NotImplemented, "Thermal conductivities"); }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        calculate its heat capacity \f$\mathrm{[J/(kg K)]}\f$.
     * \param fluidState An arbitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               int phaseIdx)
    { DUNE_THROW(Dune::NotImplemented, "Heat capacities"); }


private:
    static Scalar henryCoeffWater_(int compIdx, Scalar temperature)
    {
        // use henry's law for the solutes and the vapor pressure for
        // the solvent.
        switch (compIdx) {
            case H2OIdx: return H2O::vaporPressure(temperature);

                // the values of the Henry constant for the solutes have
                // are faked so far...
            case C1Idx: return 5.57601e+09;
            case C3Idx: return 1e10;
            case C6Idx: return 1e10;
            case C10Idx: return 1e10;
            case C15Idx: return 1e10;
            case C20Idx: return 1e10;
            default: DUNE_THROW(Dune::InvalidStateException, "Unknown component index " << compIdx);
        }
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
