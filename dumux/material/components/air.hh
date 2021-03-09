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
 * \brief A simple class for the air fluid properties
 */
#ifndef DUMUX_AIR_HH
#define DUMUX_AIR_HH

#include <dune/common/math.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/base.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/idealgas.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the air fluid properties
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Air
: public Components::Base<Scalar, Air<Scalar> >
, public Components::Gas<Scalar, Air<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    /*!
     * \brief A human readable name for Air.
     */
    static std::string name()
    { return "Air"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Air.
     *
     * Taken from constrelair.hh.
     */
    static constexpr Scalar molarMass()
    { return 0.02896; /* [kg/mol] */ }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of Air.
     */
    static Scalar criticalTemperature()
    { return 132.6312; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of Air.
     */
    static Scalar criticalPressure()
    { return 37.86e5; /* [Pa] */ }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of Air at a given pressure and temperature.
     *
     * Ideal gas is assumed.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     * \brief The molar density of air in \f$\mathrm{[mol/m^3]}\f$,
     *   depending on pressure and temperature.
     * \param temperature The temperature of the gas
     * \param pressure The pressure of the gas
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }

    /*!
     * \brief Returns true, the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true, the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure \f$\mathrm{[Pa]}\f$ of gaseous Air at a given density and temperature.
     *
     * Ideal gas is assumed.
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
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of Air at a given pressure and temperature.
     *
     * Criticial specific volume calculated by \f$V_c = (R*T_c)/p_c\f$.
     *
     * Reid et al. (1987, pp 396-397, 667) \cite reid1987 <BR>
     * Poling et al. (2001, pp 9.7-9.8) \cite poling2001 <BR>
     *
     * Accentric factor taken from: <BR>
     * Adebiyi (2003) \cite adebiyi2003
     *
     * air is a non-polar substance,
     * thus dipole moment mu is zero, as well the dimensionless dipole moment mu_r
     * therefore not considered below
     * the same holds for the correction value kappa for highly polar substances
     *
     * This calculation was introduced into Dumux in 2012 although the method here
     * is designed for general polar substances. Air, however, is (a) non-polar,
     * and (b) there are more precise methods available
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar oldGasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 84.525138; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.078; // accentric factor
        const Scalar M = molarMass() * 1e3; // molar mas [g/mol]

        const Scalar Fc = 1.0 - 0.2756*omega;
        const Scalar Tstar = 1.2593*temperature/Tc;

        using std::exp;
        using std::pow;
        const Scalar Omega_v = 1.16145*pow(Tstar, -0.14874)
                               + 0.52487*exp(-0.77320*Tstar)
                               + 2.16178*exp(-2.43787*Tstar);

        using std::cbrt;
        using std::sqrt;
        const Scalar mu = 40.785 * Fc * sqrt(M * temperature)/(cbrt(Vc * Vc) * Omega_v);

        // convertion from micro poise to Pa s
        return mu/1.0e6/10.0;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of Air at a given pressure and temperature.
     *
     * Simple method, already implemented in MUFTE-UG, but pretty accurate.
     *
     * The pressure correction is even simpler and developed and tested by
     * Holger Class in 2016 against the results of the Lemmon and Jacobsen (2004)
     * approach \cite Lemmon2004a
     * It shows very reasonable results throughout realistic pressure and
     * temperature ranges up to several hundred Kelvin and up to 500 bar
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        // above 1200 K, the function becomes inaccurate
        // since this should realistically never happen, we can live with it
        const Scalar tempCelsius = temperature - 273.15;
        const Scalar pressureCorrectionFactor = 9.7115e-9*tempCelsius*tempCelsius - 5.5e-6*tempCelsius + 0.0010809;

        using std::sqrt;
        const Scalar mu = 1.496e-6 * sqrt(temperature * temperature * temperature) / (temperature + 120.0)
                          * (1.0 + (pressure/1.0e5 - 1.0)*pressureCorrectionFactor);
        return mu;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of Air at a given pressure and temperature.
     *
     * Simple method, already implemented in MUFTE-UG, but pretty accurate
     * at atmospheric pressures.
     * Gas viscosity is not very dependent on pressure. Thus, for
     * low pressures one might switch the pressure correction off
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar simpleGasViscosity(Scalar temperature, Scalar pressure)
    {
        // above 1200 K, the function becomes inaccurate
        // since this should realistically never happen, we can live with it
        using std::sqrt;
        return 1.496e-6 * sqrt(temperature * temperature * temperature) / (temperature + 120.0);
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of Air at a given pressure and temperature.
     *
     * This is a very exact approach by Lemmon and Jacobsen (2004) \cite Lemmon2004a
     * All the values and parameters used below are explained in their paper
     * Since they use ''eta'' for dyn. viscosity, we do it as well for easier
     * comparison with the paper
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar exactGasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar epsk = 103.3; // [K]

        using std::log;
        using std::exp;
        using std::sqrt;
        const Scalar logTstar = log(temperature/epsk);
        const Scalar Omega = exp(0.431
                                 - 0.4623*logTstar
                                 + 0.08406*logTstar*logTstar
                                 + 0.005341*logTstar*logTstar*logTstar
                                 - 0.00331*logTstar*logTstar*logTstar*logTstar);

        const Scalar sigma = 0.36; // [nm]
        const Scalar eta0 = 0.0266958*sqrt(1000.0*molarMass()*temperature)/(sigma*sigma*Omega);

        using std::pow;
        using Dune::power;
        const Scalar tau = criticalTemperature()/temperature;
        const Scalar rhoc = 10.4477; // [mol/m^3]
        const Scalar delta = 0.001*pressure/(temperature*8.3144598)/rhoc;
        const Scalar etaR = 10.72 * pow(tau, 0.2) * delta
                            + 1.122 * pow(tau, 0.05) * power(delta, 4)
                            + 0.002019 * pow(tau, 2.4) * power(delta, 9)
                            - 8.876 * pow(tau, 0.6) * delta * exp(-delta)
                            - 0.02916 * pow(tau, 3.6) * power(delta, 8) * exp(-delta);

        return (eta0 + etaR)*1e-6;
    }

    /*!
     * \brief Specific enthalpy of Air \f$\mathrm{[J/kg]}\f$
     *        with 273.15 \f$ K \f$ as basis.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * Kays et al. (2005, 431ff) \cite kays2005 <BR>
     */
    static Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        return gasHeatCapacity(temperature, pressure) * (temperature-273.15);
    }

    /*!
     * \brief Specific internal energy of Air \f$\mathrm{[J/kg]}\f$.
     *
     * Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     * Rearranging for internal energy yields: \f$u = h - pv\f$.
     * Exploiting the Ideal Gas assumption
     * (\f$pv = R_{\textnormal{specific}} T\f$) gives: \f$u = h - R / M T \f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return gasEnthalpy(temperature, pressure)
               - IdealGas::R * temperature // = pressure * molar volume for an ideal gas
                 / molarMass(); // conversion from [J/(mol K)] to [J/(kg K)]
    }

    /*!
     * \brief Specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$ of pure
     *        air.
     *
     *  This methods uses the formula for "zero-pressure" heat capacity that
     *  is only dependent on temperature, because the pressure dependence is rather small.
     *  This one should be accurate for a pressure of 1 atm.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     *  Values taken from Hollis (1996) \cite hollis1996 <BR>
     *  "Tables of Thermal Properties of Gases"
     */
    static const Scalar gasHeatCapacity(Scalar temperature,
                                        Scalar pressure)
    {
        // scale temperature with reference temp of 100K
        Scalar phi = temperature/100;

        using std::pow;
        using Dune::power;
        Scalar c_p = 0.661738E+01
                -0.105885E+01 * phi
                +0.201650E+00 * power(phi,2)
                -0.196930E-01 * power(phi,3)
                +0.106460E-02 * power(phi,4)
                -0.303284E-04 * power(phi,5)
                +0.355861E-06 * power(phi,6);
        c_p +=  -0.549169E+01 * power(phi,-1)
                +0.585171E+01 * power(phi,-2)
                -0.372865E+01 * power(phi,-3)
                +0.133981E+01 * power(phi,-4)
                -0.233758E+00 * power(phi,-5)
                +0.125718E-01 * power(phi,-6);
        c_p *= IdealGas::R / molarMass(); // in J/(mol*K) / (kg/mol)

        return  c_p;
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of air.
     *
     * Isobaric Properties for Nitrogen in: NIST Standard \cite NIST <BR>
     * evaluated at p=.1 MPa, T=20°C <BR>
     * Nitrogen: 0.025398 <BR>
     * Oxygen: 0.026105 <BR>
     * lambda_air is approximately 0.78*lambda_N2+0.22*lambda_O2
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        return 0.0255535;
    }
};

} // end namespace Components
} // end namespace Dumux

#endif
