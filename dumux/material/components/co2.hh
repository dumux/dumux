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
 * \brief A class for the CO2 fluid properties
 */
#ifndef DUMUX_CO2_HH
#define DUMUX_CO2_HH

#include <cmath>
#include <iostream>

#include <dune/common/math.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the CO2 fluid properties
 *
 * Under reservoir conditions, CO2 is typically in supercritical state. These
 * properties can be provided in tabulated form, which is necessary for this
 * component implementation. The template is passed through the fluidsystem
 * brineco2fluidsystem.hh.
 * Depending on the used tabulation, the fluidsystem can also be used for gaseous CO2
 */

template <class Scalar, class CO2Tables>
class CO2
: public Components::Base<Scalar, CO2<Scalar, CO2Tables> >
, public Components::Liquid<Scalar, CO2<Scalar, CO2Tables> >
, public Components::Gas<Scalar, CO2<Scalar, CO2Tables> >
{
    static const Scalar R;

    static bool warningThrown;

public:
    /*!
     * \brief A human readable name for the CO2.
     */
    static std::string name()
    { return "CO2"; }

    /*!
     * \brief The mass in \f$\mathrm{[kg/mol]}\f$ of one mole of CO2.
     */
    static constexpr Scalar molarMass()
    { return 44e-3; /* [kg/mol] */ }

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
     * \brief Returns the minimal tabulated pressure \f$\mathrm{[Pa]}\f$ of the used table
     */
    static Scalar minTabulatedPressure()
    { return CO2Tables::tabulatedEnthalpy.minPress(); /* [Pa] */ }

    /*!
     * \brief Returns the maximal tabulated pressure \f$\mathrm{[Pa]}\f$ of the used table
     */
    static Scalar maxTabulatedPressure()
    { return CO2Tables::tabulatedEnthalpy.maxPress(); /* [Pa] */ }

    /*!
     * \brief Returns the minimal tabulated temperature \f$\mathrm{[K]}\f$ of the used table
     */
    static Scalar minTabulatedTemperature()
    { return CO2Tables::tabulatedEnthalpy.minTemp(); /* [K] */ }

    /*!
     * \brief Returns the maximal tabulated temperature \f$\mathrm{[K]}\f$ of the used table
     */
    static Scalar maxTabulatedTemperature()
    { return CO2Tables::tabulatedEnthalpy.maxTemp(); /* [K] */ }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return false; }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure CO2
     *        at a given temperature.
     * \param T the temperature \f$\mathrm{[K]}\f$
     * See:
     *
     * R. Span and W. Wagner (1996, pp. 1509-1596) \cite span1996
     */
    static Scalar vaporPressure(Scalar T)
    {
        static const Scalar a[4] =
            { -7.0602087, 1.9391218, -1.6463597, -3.2995634 };
        static const Scalar t[4] =
            { 1.0, 1.5, 2.0, 4.0 };

        // this is on page 1524 of the reference
        Scalar exponent = 0;
        Scalar Tred = T/criticalTemperature();

        using std::pow;
        for (int i = 0; i < 4; ++i)
            exponent += a[i]*pow(1 - Tred, t[i]);
        exponent *= 1.0/Tred;

        using std::exp;
        return exp(exponent)*criticalPressure();
    }

    /*!
     * \brief Specific enthalpy of gaseous CO2 \f$\mathrm{[J/kg]}\f$.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasEnthalpy(Scalar temperature,
                              Scalar pressure)
    {
        if ((temperature < criticalTemperature() || pressure < criticalPressure()) && !warningThrown)
        {
            Dune::dwarn << "Subcritical values: Be aware to use "
                        <<"Tables with sufficient resolution!"<< std::endl;
            warningThrown=true;
        }
        return
            CO2Tables::tabulatedEnthalpy.at(temperature, pressure);
    }

    /*!
     * \brief Specific enthalpy of liquid CO2 \f$\mathrm{[J/kg]}\f$.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidEnthalpy(Scalar temperature,
                                 Scalar pressure)
    {
        if ((temperature < criticalTemperature() || pressure < criticalPressure()) && !warningThrown)
        {
            Dune::dwarn << "Subcritical values: Be aware to use "
                        <<"Tables with sufficient resolution!"<< std::endl;
            warningThrown=true;
        }

        return gasEnthalpy(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of CO2 \f$\mathrm{[J/kg]}\f$.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasInternalEnergy(Scalar temperature,
                                    Scalar pressure)
    {
        Scalar h = gasEnthalpy(temperature, pressure);
        Scalar rho = gasDensity(temperature, pressure);

        return h - (pressure / rho);
    }

    /*!
     * \brief Specific internal energy of liquid CO2 \f$\mathrm{[J/kg]}\f$.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidInternalEnergy(Scalar temperature,
                                       Scalar pressure)
    {
        Scalar h = liquidEnthalpy(temperature, pressure);
        Scalar rho = liquidDensity(temperature, pressure);

        return h - (pressure / rho);
    }

    /*!
     * \brief The density of CO2 at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        if ((temperature < criticalTemperature() || pressure < criticalPressure()) && !warningThrown)
        {
            Dune::dwarn << "Subcritical values: Be aware to use "
                        <<"Tables with sufficient resolution!"<< std::endl;
            warningThrown=true;
        }
        return CO2Tables::tabulatedDensity.at(temperature, pressure);
    }

    /*!
     *  \brief The molar density of CO2 gas in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return gasDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief The density of pure CO2 at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        if ((temperature < criticalTemperature() || pressure < criticalPressure()) && !warningThrown)
        {
            Dune::dwarn << "Subcritical values: Be aware to use "
                        <<"Tables with sufficient resolution!"<< std::endl;
            warningThrown=true;
        }
        return CO2Tables::tabulatedDensity.at(temperature, pressure);
    }

    /*!
     * \brief The molar density of CO2 in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        DUNE_THROW(NumericalProblem, "CO2::gasPressure()");
    }

    /*!
     * \brief The pressure of liquid water in \f$\mathrm{[Pa]}\f$ at a given density and
     *        temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        DUNE_THROW(NumericalProblem, "CO2::liquidPressure()");
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a liquid.
     * USE WITH CAUTION! Exploits enthalpy function with artificial increment
     * of the temperature!
     * Equation with which the specific heat capacity is calculated : \f$ c_p = \frac{dh}{dT}\f$
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        //temperature difference :
        Scalar dT = 1.; // 1K temperature increment
        Scalar temperature2 = temperature+dT;

        // enthalpy difference
        Scalar hold = liquidEnthalpy(temperature, pressure);
        Scalar hnew = liquidEnthalpy(temperature2, pressure);
        Scalar dh = hold-hnew;

        //specific heat capacity
        return dh/dT ;
    }


    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of CO2.
     * Equations given in: - Vesovic et al., 1990
     *                     - Fenhour et al., 1998
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        static const double a0 = 0.235156;
        static const double a1 = -0.491266;
        static const double a2 = 5.211155E-2;
        static const double a3 = 5.347906E-2;
        static const double a4 = -1.537102E-2;

        static const double d11 = 0.4071119E-2;
        static const double d21 = 0.7198037E-4;
        static const double d64 = 0.2411697E-16;
        static const double d81 = 0.2971072E-22;
        static const double d82 = -0.1627888E-22;

        static const double ESP = 251.196;

        double mu0, SigmaStar, TStar;
        double dmu, rho;
        double visco_CO2;

        if(temperature < 275.) // regularisation
        {
            temperature = 275;
            Dune::dgrave << "Temperature below 275K in viscosity function:"
                    << "Regularizing tempereature to 275K. " << std::endl;
        }


        TStar = temperature/ESP;

        /* mu0: viscosity in zero-density limit */
        using std::exp;
        using std::log;
        using std::sqrt;
        SigmaStar = exp(a0 + a1*log(TStar)
                        + a2*log(TStar)*log(TStar)
                        + a3*log(TStar)*log(TStar)*log(TStar)
                        + a4*log(TStar)*log(TStar)*log(TStar)*log(TStar) );
        mu0 = 1.00697*sqrt(temperature) / SigmaStar;

        /* dmu : excess viscosity at elevated density */
        rho = gasDensity(temperature, pressure); /* CO2 mass density [kg/m^3] */

        using Dune::power;
        dmu = d11*rho + d21*rho*rho + d64*power(rho,6)/(TStar*TStar*TStar)
            + d81*power(rho,8) + d82*power(rho,8)/TStar;

        visco_CO2 = (mu0 + dmu)/1.0E6;   /* conversion to [Pa s] */

        return visco_CO2;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure CO2.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        // no difference for supercritical CO2
        return gasViscosity(temperature, pressure);
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
};

template <class Scalar, class CO2Tables>
const Scalar CO2<Scalar, CO2Tables>::R = Constants<Scalar>::R;

template <class Scalar, class CO2Tables>
bool CO2<Scalar, CO2Tables>::warningThrown = false;

} // end namespace Components

} // end namespace Dumux

#endif
