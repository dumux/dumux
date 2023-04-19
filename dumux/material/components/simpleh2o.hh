// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief A simple implementation of pure water
 */
#ifndef DUMUX_SIMPLE_H2O_HH
#define DUMUX_SIMPLE_H2O_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/idealgas.hh>

#include <cmath>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux::Components {

/*!
 * \ingroup Components
 * \brief A simple implementation of pure water
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class SimpleH2O
: public Components::Base<Scalar, SimpleH2O<Scalar> >
, public Components::Liquid<Scalar, SimpleH2O<Scalar> >
, public Components::Gas<Scalar, SimpleH2O<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    /*!
     * \brief A human readable name for the water.
     */
    static std::string name()
    { return "SimpleH2O"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of water.
     */
    static constexpr Scalar molarMass()
    { return 18e-3; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of water.
     */
    static Scalar criticalTemperature()
    { return 647.096; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of water.
     */
    static Scalar criticalPressure()
    { return 22.064e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at water's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.16; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at water's triple point.
     */
    static Scalar triplePressure()
    { return 611.657; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static Scalar vaporPressure(Scalar T)
    {
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // water is solid: We don't take sublimation into account

        constexpr Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        const Scalar sigma = T + n[8]/(T - n[9]);

        const Scalar A = (sigma + n[0])*sigma + n[1];
        const Scalar B = (n[2]*sigma + n[3])*sigma + n[4];
        const Scalar C = (n[5]*sigma + n[6])*sigma + n[7];

        using std::sqrt;
        const Scalar term = 2.0*C/(sqrt(B*B - 4.0*A*C) - B);

        return 1e6*term*term*term*term;
    }

    /*!
     * \brief Specific enthalpy of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        static const Scalar tRef = getParam<Scalar>("SimpleH2O.ReferenceTemperature", 293.15);
        return gasHeatCapacity(temperature, pressure)*(temperature - tRef) + vaporizationEnthalpy();
    }

    /*!
     * \brief Specific enthalpy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        static const Scalar tRef = getParam<Scalar>("SimpleH2O.ReferenceTemperature", 293.15);
        return liquidHeatCapacity(temperature, pressure)*(temperature - tRef)
                + pressure/liquidDensity(temperature, pressure);
    }

   /*!
    * \brief The vaporization enthalpy in \f$\mathrm{[J/kg]}\f$ needed to vaporize one kilogram of the liquid water to the gaseous state depending on temperature as found in: C. O. Popiel & J. Wojtkowiak (1998) Simple Formulas for Thermophysical Properties of Liquid Water for Heat Transfer Calculations (from 0°C to 150°C), DOI:10.1080/01457639808939929
    */
    static Scalar vaporizationEnthalpy()
    {
        constexpr Scalar A = 2500.304;
        constexpr Scalar B = -2.2521025;
        constexpr Scalar C = -0.021465847;
        constexpr Scalar D = 3.1750136e-4 ;
        constexpr Scalar E = -2.8607959e-5;

        //tRef in °C
        static const Scalar tRef = getParam<Scalar>("SimpleH2O.ReferenceTemperature", 293.15) - 273.15;

        using std::pow;
        static const Scalar vaporizationEnthalpy = A + B*tRef + C*(pow(tRef, 1.5)) + D*(pow(tRef, 2.5)) + E*(pow(tRef, 3));
        return vaporizationEnthalpy;
    }


    /*!
     * \brief Specific internal energy of steam \f$\mathrm{[J/kg]}\f$.
     *
     *        Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     *        Rearranging for internal energy yields: \f$u = h - pv\f$.
     *        Exploiting the Ideal Gas assumption (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     *        The universal gas constant can only be used in the case of molar formulations.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        // 1/molarMass: conversion from [J/(mol K)] to [J/(kg K)]
        // R*T/molarMass: pressure *spec. volume for an ideal gas
        return gasEnthalpy(temperature, pressure)
                - 1.0/molarMass()*IdealGas::R*temperature;
    }

    /*!
     * \brief Specific internal energy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    {
        return liquidEnthalpy(temperature, pressure)
                - pressure/liquidDensity(temperature, pressure);
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
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     *  \brief The molar density of steam in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief The density of pure water at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1000.0;
    }

    /*!
     * \brief The molar density of pure water in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief The pressure of water in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The liquid pressure is undefined for incompressible fluids");
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        return 1e-05;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 1e-03;
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a liquid.
     * source: http://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&Action=Page
     * @ T= 281.15K (8°C) , p=0.1MPa)
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        return 4180.0;
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of water.
     *        source: http://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&Action=Page
     *        @ T= 372.76K (99.6°C) , p=0.1MPa)
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    {
       return 0.679;
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of steam.
     *        source: http://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&Action=Page
     *        @ T= 372.76K (99.6°C) , p=0.1MPa)
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
       return 0.025;
    }

    /*!
     * \brief Specific isobaric heat capacity of water steam \f$\mathrm{[J/(kg*K)]}\f$.
     *        source: http://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&Action=Page
     *        @ T= 372.76K (99.6°C) , p=0.1MPa)
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    {
        return 2.08e3;
    }

};

template <class Scalar>
struct IsAqueous<SimpleH2O<Scalar>> : public std::true_type {};

} // end namespace Dumux::Components

#endif
