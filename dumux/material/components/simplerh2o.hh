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
 * \ingroup Components
 * \brief A simpler version of pure water
 */
#ifndef DUMUX_SIMPLER_H2O_HH
#define DUMUX_SIMPLER_H2O_HH

#include <iostream>
#include <ostream>

#include <dumux/common/parameters.hh>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

#include <dumux/material/idealgas.hh>
#include <dumux/material/components/h2o.hh>

namespace Dumux::Components {

/*!
 * \ingroup Components
 * \brief A simple version of water
 * This water version uses constant or linear functions for most quantities
 * All parameters initialized once based on the H2OType and a user-provided
 * temperature and pressure. Using this H2O makes sense if temperature and
 * pressure only vary in a small range and most values can therefore be
 * assumed constants. In order for the values to best match your simulation
 * scenario, adjust the temperature and pressure by setting the parameters
 * SimplerH2O.Pressure and SimplerH2O.Temperature.
 * Water vapor is assumed to be an ideal gas.
 *
 * \tparam Scalar The type used for scalar values
 * \tparam H2OType The base water component type we extract values from
 */
template <class Scalar, class H2OType = H2O<Scalar>>
class SimplerH2O
: public Components::Base<Scalar, SimplerH2O<Scalar> >
, public Components::Liquid<Scalar, SimplerH2O<Scalar> >
, public Components::Gas<Scalar, SimplerH2O<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;
    using H2O = H2OType;

public:
    /*!
     * \brief The temperature paramater (see class description)
     */
    static Scalar temperature()
    {
        static const Scalar T = getParam<Scalar>("SimplerH2O.Temperature", 293.15);
        return T;
    }

    /*!
     * \brief The pressure paramater (see class description)
     */
    static Scalar pressure()
    {
        static const Scalar p = getParam<Scalar>("SimplerH2O.Pressure", 1.0e5);
        return p;
    }

    /*!
     * \brief A human readable name for the water.
     */
    static std::string name()
    { return "SimplerH2O"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of water.
     */
    static constexpr Scalar molarMass()
    { return H2O::molarMass(); }

    /*!
     * \brief The acentric factor \f$\mathrm{[-]}\f$ of water.
     */
    static constexpr Scalar acentricFactor()
    { return H2O::acentricFactor(); }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of water
     */
    static constexpr Scalar criticalTemperature()
    { return H2O::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of water.
     */
    static constexpr Scalar criticalPressure()
    { return H2O::criticalPressure(); }

    /*!
     * \brief Returns the molar volume \f$\mathrm{[m^3/mol]}\f$ of water at the critical point
     */
    static constexpr Scalar criticalMolarVolume()
    { return H2O::criticalMolarVolume(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at water's triple point.
     */
    static constexpr Scalar tripleTemperature()
    { return H2O::tripleTemperature(); }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at water's triple point.
     */
    static constexpr Scalar triplePressure()
    { return H2O::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    {
        static const dp = [&]{
            const auto epsT = temperature()*1e-10;
            const Scalar u0 = H20::vaporPressure(temperature());
            const Scalar u1 = H20::vaporPressure(temperature()+epsT);
            return (u1-u0)/epsT;
        }();

        return dp*(T - temperature());
    }

    /*!
     * \brief Specific enthalpy of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar T, Scalar p)
    { return gasInternalEnergy(T, p) + 1.0/molarMass()*IdealGas::R*T; }

    /*!
     * \brief Specific enthalpy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar T, Scalar p)
    { return H2O::liquidInternalEnergy(T, p) + p/liquidDensity(T, p); }

    /*!
     * \brief Specific internal energy of steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar T, Scalar p)
    {
        static const du = [&]{
            const auto epsT = temperature()*1e-10;
            const Scalar u0 = H20::gasInternalEnergy(temperature(), pressure());
            const Scalar u1 = H20::gasInternalEnergy(temperature()+epsT, pressure());
            return (u1-u0)/epsT;
        }();

        return du*(T - temperature());
    }

    /*!
     * \brief Specific internal energy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar T, Scalar p)
    {
        static const du = [&]{
            const auto epsT = temperature()*1e-10;
            const Scalar u0 = H20::liquidInternalEnergy(temperature(), pressure());
            const Scalar u1 = H20::liquidInternalEnergy(temperature()+epsT, pressure());
            return (u1-u0)/epsT;
        }();

        return du*(T - temperature());
    }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Returns true if the gas phase viscosity is constant
     */
    static constexpr bool gasViscosityIsConstant()
    { return true; }

    /*!
     * \brief Returns true if the liquid phase viscosity is constant
     */
    static constexpr bool liquidViscosityIsConstant()
    { return true; }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of steam at a given pressure and temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar T, Scalar p)
    { return IdealGas::density(molarMass(), T, p); }

    /*!
     *  \brief The molar density of steam in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasMolarDensity(Scalar T, Scalar p)
    { return IdealGas::molarDensity(T, p); }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar T, Scalar density)
    { return IdealGas::pressure(T, density/molarMass()); }

    /*!
     * \brief The density of pure water at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar T, Scalar p)
    { return liquidDensity(); }

    static Scalar liquidDensity()
    { return H2O::liquidDensity(temperature(), pressure()); }

    /*!
     * \brief The molar density of pure water in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar T, Scalar p)
    { return liquidDensity(T, p)/molarMass(); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar T, Scalar p)
    { return gasViscosity(); }

    static Scalar gasViscosity()
    { return H2O::gasViscosity(temperature(), pressure()); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar T, Scalar p)
    { return liquidViscosity(); }

    static Scalar liquidViscosity()
    { return H2O::liquidViscosity(temperature(), pressure()); }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a liquid
     * \param T absolute temperature in \f$\mathrm{[K]}\f$
     * \param p pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidHeatCapacity(Scalar T, Scalar p)
    { return liquidHeatCapacity(); }

    static Scalar liquidHeatCapacity()
    { return H2O::liquidHeatCapacity(temperature(), pressure()); }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of water
     * \param T absolute temperature in \f$\mathrm{[K]}\f$
     * \param p pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity(Scalar T, Scalar p)
    { return liquidThermalConductivity(); }

    static Scalar liquidThermalConductivity()
    { return H2O::liquidThermalConductivity(temperature(), pressure()); }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of steam
     * \param T absolute temperature in \f$\mathrm{[K]}\f$
     * \param p pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar T, Scalar p)
    { return gasThermalConductivity(); }

    static Scalar gasThermalConductivity()
    { return H2O::gasThermalConductivity(temperature(), pressure()); }

    /*!
     * \brief Specific isobaric heat capacity of water steam \f$\mathrm{[J/(kg*K)}\f$
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasHeatCapacity(Scalar T, Scalar p)
    { return gasHeatCapacity(); }

    static Scalar gasHeatCapacity()
    { return H2O::gasHeatCapacity(temperature(), pressure()); }

    /*!
     * \brief Create a report of the reference state
     */
    static void reportReference(std::ostream& stream = std::cout)
    {
        stream << "SimplerH2O reference state:\n"
               << "-- temperature: " << temperature() << " K"
               << "-- pressure: " << pressure() << " Pa\n"
               << std::endl;
    }

    /*!
     * \brief Create a report of the computed liquid quantities for the reference state
     */
    static void reportLiquid(std::ostream& stream = std::cout)
    {
        stream << "SimplerH2O configuration (liquid):\n"
               << "-- density: " << liquidDensity() << " kg/m^3\n"
               << "-- viscosity: " << liquidViscosity() << " Pa*s\n"
               << "-- heat capacity: " << liquidHeatCapacity() << " J/(kg*K)\n"
               << "-- thermal conductivity: " << liquidThermalConductivity() << " W(m*K)\n"
               << "-- enthalpy, internal energy and vapor pressure linearized around reference (T,p)\n"
               << std::endl;
    }

    /*!
     * \brief Create a report of the computed gas quantities for the reference state
     */
    static void reportGas(std::ostream& stream = std::cout)
    {
        stream << "SimplerH2O configuration (gas):\n"
               << "-- density: " << gasDensity() << " kg/m^3\n"
               << "-- viscosity: " << gasViscosity() << " Pa*s\n"
               << "-- heat capacity: " << gasHeatCapacity() << " J/(kg*K)\n"
               << "-- thermal conductivity: " << gasThermalConductivity() << " W(m*K)\n"
               << "-- enthalpy, internal energy and vapor pressure linearized around reference (T,p)\n"
               << std::endl;
    }

    /*!
     * \brief Report all quantities
     */
    static void report(std::ostream& stream = std::cout)
    {
        reportReference(stream);
        reportLiquid(stream);
        reportGas(stream);
    }
};

template <class Scalar>
struct IsAqueous<SimplerH2O<Scalar>> : public std::true_type {};

} // end namespace Dumux::Components

#endif
