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
 * \brief @copybrief Dumux::FluidSystems::TwoPLiquidVaporFluidsystem
 */
#ifndef DUMUX_2P_LIQUID_VAPOR_FLUID_SYSTEM_HH
#define DUMUX_2P_LIQUID_VAPOR_FLUID_SYSTEM_HH

#include <limits>
#include <cassert>
#include <iostream>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include <dune/common/exceptions.hh>

#include "base.hh"

namespace Dumux
{
namespace FluidSystems
{

/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with only one component.
 *
 *
 * This FluidSystem can be used without the PropertySystem that is applied in Dumux,
 * as all Parameters are defined via template parameters. Hence it is in an
 * additional namespace Dumux::FluidSystem::.
 * An adapter class using Dumux::FluidSystem<TypeTag> is also provided
 * at the end of this file.
 */
template <class Scalar, class ComponentType>
class TwoPLiquidVaporFluidsystem
    : public BaseFluidSystem<Scalar, TwoPLiquidVaporFluidsystem<Scalar, ComponentType> >
{
    // do not try to instantiate this class, it has only static members!
    TwoPLiquidVaporFluidsystem()
    {}

    typedef TwoPLiquidVaporFluidsystem<Scalar, ComponentType> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;
    typedef ComponentType Component;

public:
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 2;

    static constexpr int wPhaseIdx = 0; // index of the wetting phase
    static constexpr int nPhaseIdx = 1; // index of the non-wetting phase


    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        static std::string name[] = {
            std::string("liquid"),
            std::string("vapor"),
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
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
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        if (phaseIdx == nPhaseIdx)
            return true;
        // the component decides for the liquid phase...
        return Component::liquidIsCompressible();
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
            return Component::gasIsIdeal();
        return false; // not a gas
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 1;

    static constexpr int ComponentIdx = 0;



    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        static std::string name[] = {
            Component::name()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        static const Scalar M[] = {
            Component::molarMass()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return M[compIdx];
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        static const Scalar Tcrit[] = {
            Component::criticalTemperature()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return Tcrit[compIdx];
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        static const Scalar pcrit[] = {
            Component::criticalPressure()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    }

    /*!
     * \brief Molar volume of a component at the critical point [m^3/mol].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "TwoPLiquidVaporFluidsystem::criticalMolarVolume()");
        return 0;
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            Component::acentricFactor(), // H2O (from Reid, et al.)
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return accFac[compIdx];
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
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (Component::isTabulated) {
            std::cout << "Initializing tables for the fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            Component::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    /*!
     * \brief Calculate the density [kg/m^3] of a fluid phase
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx) {
                return Component::liquidDensity(T, p);
        }
        else if (phaseIdx == nPhaseIdx)// gas phase
        {
            return Component::gasDensity(T, p);
        }
        else DUNE_THROW(Dune::NotImplemented,
                   "wrong index");
    }

    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s]
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx) {
                return Component::liquidViscosity(T, p);
        }
        else if (phaseIdx == nPhaseIdx) // gas phase
        {
            return Component::gasViscosity(T, p) ;
        }
        else DUNE_THROW(Dune::NotImplemented,
                   "wrong index");
    }

    /*!
     * \brief calculate the temperature of vapor at a given pressure on the vapor pressure curve.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar vaporTemperature(const FluidState &fluidState,
                                   const unsigned int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        Scalar pressure = fluidState.pressure(nPhaseIdx) ;

        return Component::vaporTemperature( pressure ) ;
    }

    /*!
     * \brief Calculate the fugacity coefficient [-] of an individual
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
    using Base::fugacityCoefficient;
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
           return Component::vaporPressure(T)/p;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        else
        {
            return 1.0;

        }
    }


    /*!
     * \brief Calculate the molecular diffusion coefficient for a
     *        component in a fluid phase [mol^2 * s / (kg*m^3)]
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    using Base::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        // TODO!
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient for components
     *        \f$i\f$ and \f$j\f$ in this phase.
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
        DUNE_THROW(Dune::NotImplemented, "Binary Diffusion coefficients");
    }

    /*!
     * \brief Calculate specific enthalpy [J/kg].
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using Base::enthalpy;
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        // liquid phase
        if (phaseIdx == wPhaseIdx) {
            return Component::liquidEnthalpy(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == nPhaseIdx) // gas phase
        {
            return Component::gasEnthalpy(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
        }
        else DUNE_THROW(Dune::NotImplemented,
                   "wrong index");
    }

    /*!
     * \brief Thermal conductivity of a fluid phase [W/(m K)].
     *
     * Use the conductivity of vapor and liquid water at 100°C
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using Base::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        // liquid phase
        if (phaseIdx == wPhaseIdx) {
                 return Component::liquidThermalConductivity(fluidState.temperature(phaseIdx),
                                                        fluidState.pressure(phaseIdx)); //0.68 ;
        }
        else if (phaseIdx == nPhaseIdx) // gas phase
        {
            return Component::gasThermalConductivity(fluidState.temperature(phaseIdx),
                                                        fluidState.pressure(phaseIdx)); //0.0248;
        }
        else DUNE_THROW(Dune::NotImplemented,
                   "wrong index");
    }

    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/kg / K]}\f$.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using Base::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        // liquid phase
        if (phaseIdx == wPhaseIdx) {
                 return Component::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));//4.217e3 ;
        }
        else if (phaseIdx == nPhaseIdx) // gas phase
        {
            return Component::gasHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));//2.029e3;
        }
        else DUNE_THROW(Dune::NotImplemented,
                   "wrong index");
    }
};

} // end namespace FluidSystems

#ifdef DUMUX_PROPERTIES_HH
// forward definitions of the property tags
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Components);
}
/*!
 * \brief A two-phase fluid system with only one component.
 *
 * This is an adapter to use Dumux::TwoPLiquidVaporFluidsystem<TypeTag>, as is
 * done with most other classes in Dumux.
 */
template<class TypeTag>
class TwoPLiquidVaporFluidsystem
: public FluidSystems::TwoPLiquidVaporFluidsystem<typename GET_PROP_TYPE(TypeTag, Scalar),
                             typename GET_PROP_VALUE(TypeTag, Components)::Component>
{};
#endif

} // end namespace

#endif
