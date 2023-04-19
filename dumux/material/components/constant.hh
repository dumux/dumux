// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Setting constant fluid properties via the input file.
 */
#ifndef DUMUX_COMPONENTS_CONSTANT_HH
#define DUMUX_COMPONENTS_CONSTANT_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/parameters.hh>

#include <dumux/material/idealgas.hh>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A component which returns run time specified values
 *        for all fluid properties.
 *
 * \tparam id  The id used to read from the input file / parametertree
 * \tparam Scalar  The type used for scalar values
 *
 * \note For the constant component with id=1 you would specify the parameters in the input file as follows
 *       \code{.ini}
 *       [1.Component]
 *       MolarMass = 0.018 # kg/mol
 *       \endcode
 * \note If you only have one component you can also omit the "1.".
 */
template<int id, class Scalar>
class Constant
: public Components::Base<Scalar, Constant<id, Scalar> >
, public Components::Liquid<Scalar, Constant<id, Scalar> >
, public Components::Gas<Scalar, Constant<id, Scalar> >
, public Components::Solid<Scalar, Constant<id, Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return false; }

    /*!
     * \brief Returns true if the gas phase viscosity is constant
     */
    static constexpr bool gasViscosityIsConstant()
    { return true; }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Returns true if the liquid phase viscosity is constant
     */
    static constexpr bool liquidViscosityIsConstant()
    { return true; }

    /*!
     * \brief A human readable name for the component.
     */
    static const std::string& name()
    {
        static const std::string name = getParamFromGroup<std::string>(std::to_string(id), "Component.Name", "component");
        return name;
    }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass()
    {
        static const Scalar molarMass = getParamFromGroup<Scalar>(std::to_string(id), "Component.MolarMass");
        return molarMass;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at the components's triple point.
     */
    static Scalar tripleTemperature()
    {
        static const Scalar tripleTemperature = getParamFromGroup<Scalar>(std::to_string(id), "Component.TripleTemperature");
        return tripleTemperature;
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure()
    {
        static const Scalar triplePressure = getParamFromGroup<Scalar>(std::to_string(id), "Component.TriplePressure");
        return triplePressure;
    }

    /*!
     * \brief The vaporization enthalpy in \f$\mathrm{[J/kg]}\f$ needed to vaporize one kilogram of the liquid component to the gaseous state
     */
    static Scalar vaporizationEnthalpy()
    {
        static const Scalar vaporizationEnthalpy = getParamFromGroup<Scalar>(std::to_string(id), "Component.EnthalpyOfVaporization");
        return vaporizationEnthalpy;
    }


    /*!
     * \brief Sets the liquid density in \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        static const Scalar density = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidDensity");
        return density;
    }

    /*!
     * \brief The molar density in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief Sets the liquid dynamic viscosity in \f$\mathrm{[Pa*s]}\f$.
     *
     * \note We look for Component.LiquidKinematicViscosity or Component.LiquidDynamicViscosity.
     *       If both parameters are specified, it is considered a configuration error because it can
     *       be ambiguous if defaults are specified for several constant components in the plain
     *       "Component" group (without ID-prefix).
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        static const Scalar dynamicViscosity = [&]
        {
            if (hasParamInGroup(std::to_string(id), "Component.LiquidKinematicViscosity"))
            {
                if (hasParamInGroup(std::to_string(id), "Component.LiquidDynamicViscosity"))
                    DUNE_THROW(Dune::InvalidStateException, "Found both Component.LiquidKinematicViscosity and Component.LiquidDynamicViscosity."
                        << " Please only specify either the kinematic or the dynamic viscosity for all constant components to avoid ambiguities.");

                return getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidKinematicViscosity") * liquidDensity(temperature, pressure);
            }
            else
                return getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidDynamicViscosity");
        }();

        return dynamicViscosity;
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a liquid.
     * \param temperature temperature of phase in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    {
        static const Scalar thermalConductivity = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidThermalConductivity");
        return thermalConductivity;
    }

    /*!
     * \brief Specific internal energy of the component \f$\mathrm{[J/kg]}\f$ as a liquid.
     * \param temperature temperature of phase in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    {
        // u = c * dT for incompressible fluids
        const Scalar heatCapacity = liquidHeatCapacity(temperature, pressure);
        static const Scalar tRef = getParamFromGroup<Scalar>(std::to_string(id), "Component.ReferenceTemperature", 293.15);
        return heatCapacity * (temperature - tRef);
    }

    /*!
     * \brief Specific enthalpy of the component \f$\mathrm{[J/kg]}\f$ as a liquid.
     *
     * \param temperature temperature of phase in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    {
        const Scalar u = liquidInternalEnergy(temperature, pressure);
        const Scalar rho = liquidDensity(temperature, pressure);
        return u + pressure / rho;
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a liquid.
     *
     * \param temperature temperature of phase in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        static const Scalar heatCapacity = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidHeatCapacity");
        return heatCapacity;
    }

    /*!
     * \brief Sets the gas density in \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        static const Scalar density = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasDensity");
        return density;
    }

    /*!
     * \brief The molar density in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return gasDensity(temperature, pressure)/molarMass(); }


    /*!
     * \brief Sets the gas dynamic viscosity in \f$\mathrm{[Pa*s]}\f$.
     *
     * \note We look for Component.GasKinematicViscosity or Component.GasDynamicViscosity.
     *       If both parameters are specified, it is considered a configuration error because it can
     *       be ambiguous if defaults are specified for several constant components in the plain
     *       "Component" group (without ID-prefix).
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        static const Scalar dynamicViscosity = [&]
        {
            if (hasParamInGroup(std::to_string(id), "Component.GasKinematicViscosity"))
            {
                if (hasParamInGroup(std::to_string(id), "Component.GasDynamicViscosity"))
                    DUNE_THROW(Dune::InvalidStateException, "Found both Component.GasKinematicViscosity and Component.GasDynamicViscosity."
                        << " Please only specify either the kinematic or the dynamic viscosity for all constant components to avoid ambiguities.");

                return getParamFromGroup<Scalar>(std::to_string(id), "Component.GasKinematicViscosity") * gasDensity(temperature, pressure);
            }
            else
                return getParamFromGroup<Scalar>(std::to_string(id), "Component.GasDynamicViscosity");
        }();

        return dynamicViscosity;
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a gas.
     * \param temperature temperature of phase in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        static const Scalar thermalConductivity = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasThermalConductivity");
        return thermalConductivity;
    }

    /*!
     * \brief Specific internal energy of the component \f$\mathrm{[J/kg]}\f$ as a gas.
     *
     *        Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     *
     *        Rearranging for internal energy yields: \f$u = h - pv\f$.
     *
     *        Exploiting the Ideal Gas assumption (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     *        The universal gas constant can only be used in the case of molar formulations.
     * \param temperature temperature of phase in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        // 1/molarMass: conversion from [J/(mol K)] to [J/(kg K)]
        // R*T/molarMass: pressure *spec. volume for an ideal gas
        return gasEnthalpy(temperature, pressure) - 1/molarMass()* IdealGas::R*temperature;

    }

    /*!
     * \brief Specific enthalpy of the component \f$\mathrm{[J/kg]}\f$ as a gas.
     *
     * \param temperature temperature of phase in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        static const Scalar tRef = getParamFromGroup<Scalar>(std::to_string(id), "Component.ReferenceTemperature", 293.15);
        return gasHeatCapacity(temperature, pressure)*(temperature - tRef) + vaporizationEnthalpy();
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a gas.
     *
     * \param temperature temperature of phase in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    {
        static const Scalar heatCapacity = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasHeatCapacity");
        return heatCapacity;
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of a the component
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * We use the Clausius-Clapeyron Equation to estimate the vapor pressure. Vapor pressure depends on the enthalpy of vaporization. We use the triple point pressure and temperature as a reference point.
     */
    static Scalar vaporPressure(Scalar T)
    {
        const Scalar p2 = triplePressure();
        const Scalar T2 = tripleTemperature();
        const Scalar exponent = -(vaporizationEnthalpy()*molarMass())/IdealGas::R*(1/T - 1/T2);

        using std::exp;
        const Scalar vaporPressure = p2*exp(exponent);
        return vaporPressure;
    }

        /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidDensity(Scalar temperature)
    {
        static const Scalar density = getParamFromGroup<Scalar>(std::to_string(id), "Component.SolidDensity");
        return density;
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    {
        static const Scalar solidThermalConductivity = getParamFromGroup<Scalar>(std::to_string(id), "Component.SolidThermalConductivity");
        return solidThermalConductivity;
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    {
        static const Scalar solidHeatCapacity = getParamFromGroup<Scalar>(std::to_string(id), "Component.SolidHeatCapacity");
        return solidHeatCapacity;
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
