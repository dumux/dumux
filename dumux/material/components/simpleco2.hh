// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief A much simpler (and thus potentially less buggy) version of
 *        pure CO2.
 */
#ifndef DUMUX_SIMPLE_CO2_HH
#define DUMUX_SIMPLE_CO2_HH

#include <dune/common/stdstreams.hh>

#include <dumux/common/parameters.hh>
#include <dumux/material/idealgas.hh>

#include <cmath>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux::Components {

/*!
 * \ingroup Components
 * \brief A simple version of pure CO2
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class SimpleCO2
: public Components::Base<Scalar, SimpleCO2<Scalar> >
, public Components::Gas<Scalar, SimpleCO2<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    /*!
     * \brief A human readable name for the CO2.
     */
    static std::string name()
    { return "SimpleCO2"; }

    /*!
     * \brief The mass in \f$\mathrm{[kg/mol]}\f$ of one mole of CO2.
     */
    static constexpr Scalar molarMass()
    { return 44.0e-3; /* [kg/mol] */ }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of CO2
     */
    static Scalar criticalTemperature()
    { return 273.15 + 30.95; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of CO2
     */
    static Scalar criticalPressure()
    { return 73.8e5; /* [Pa] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at CO2's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.15 - 56.35; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at CO2's triple point.
     */
    static Scalar triplePressure()
    { return 5.11e5; /* [N/m^2] */ }


    /*!
     * \brief Specific enthalpy of CO2 \f$\mathrm{[J/kg]}\f$.
     *        source: Shomate Equation for a temperature range of 298. to 1200K.
     *        with components published by NIST  \cite NIST
     *        https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1&Type=JANAFG&Table=on
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        const Scalar t = temperature/1000;
        constexpr double a = 24.99735;
        constexpr double b = 55.18696;
        constexpr double c = -33.69137;
        constexpr double d = 7.948387;
        constexpr double e = -0.136638;
        constexpr double f = -403.6075;
        constexpr double h = -393.5224;
        return (a*t + b*t*t/2 + c*t*t*t/3 + d*t*t*t*t/4 - e/t +f -h)*1000/molarMass(); //conversion from kJ/mol to J/kg
    }

    /*!
     * \brief Specific internal energy of CO2 \f$\mathrm{[J/kg]}\f$.
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
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true if the gas phase viscostiy is constant
     */
    static constexpr bool gasViscosityIsConstant()
    { return false; }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of CO2 at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::density(molarMass(), temperature, pressure); }

    /*!
     *  \brief The molar density of CO2 in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
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
     * \brief The pressure of CO2 in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    { return IdealGas::pressure(temperature, density/molarMass()); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of CO2.
     * Equations given in: - Vesovic et al., 1990
     *                     - Fenhour et al., 1998
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * TODO: this does not look like a really "simple" parameterization. Can this be simplified further?
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        constexpr double a0 = 0.235156;
        constexpr double a1 = -0.491266;
        constexpr double a2 = 5.211155E-2;
        constexpr double a3 = 5.347906E-2;
        constexpr double a4 = -1.537102E-2;

        constexpr double d11 = 0.4071119E-2;
        constexpr double d21 = 0.7198037E-4;
        constexpr double d64 = 0.2411697E-16;
        constexpr double d81 = 0.2971072E-22;
        constexpr double d82 = -0.1627888E-22;

        constexpr double ESP = 251.196;

        if(temperature < 275.0) // regularisation
        {
            temperature = 275.0;
            Dune::dgrave << "Temperature below 275K in viscosity function:"
                    << "Regularizing temperature to 275K. " << std::endl;
        }


        const double TStar = temperature/ESP;

        /* mu0: viscosity in zero-density limit */
        using std::exp;
        using std::log;
        using std::sqrt;
        const double logTStar = log(TStar);
        const double SigmaStar = exp(a0 + a1*logTStar
                        + a2*logTStar*logTStar
                        + a3*logTStar*logTStar*logTStar
                        + a4*logTStar*logTStar*logTStar*logTStar );
        const double mu0 = 1.00697*sqrt(temperature) / SigmaStar;

        /* dmu : excess viscosity at elevated density */
        const double rho = gasDensity(temperature, pressure); /* CO2 mass density [kg/m^3] */

        using Dune::power;
        const double dmu = d11*rho + d21*rho*rho + d64*power(rho, 6)/(TStar*TStar*TStar)
            + d81*power(rho, 8) + d82*power(rho, 8)/TStar;

        return (mu0 + dmu)/1.0E6;   /* conversion to [Pa s] */
    }


    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of CO2.
     *
     * Thermal conductivity of CO2 at T=20°C, see:
     * http://www.engineeringtoolbox.com/carbon-dioxide-d_1000.html
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        return 0.087;
    }

    /*!
     * \brief Specific isobaric heat capacity of CO2 \f$\mathrm{[J/(kg*K)]}\f$.
     *        source: Shomate Equation for a temperature range of 298. to 1200K.
     *        with components published by NIST  \cite NIST
     *        https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1&Type=JANAFG&Table=on
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    {
        const Scalar t = temperature/1000;
        constexpr double a = 24.99735;
        constexpr double b = 55.18696;
        constexpr double c = -33.69137;
        constexpr double d = 7.948387;
        constexpr double e = -0.136638;
        return (a + b*t + c*t*t + d*t*t*t + e/(t*t))/molarMass();
    }

};

} // end namespace Dumux::Components

#endif
