// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Properties of xylene.
 */
#ifndef DUMUX_XYLENE_HH
#define DUMUX_XYLENE_HH

#include <cmath>
#include <dumux/material/idealgas.hh>
#include <dumux/material/constants.hh>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Properties of xylene.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Xylene
: public Components::Base<Scalar, Xylene<Scalar> >
, public Components::Liquid<Scalar, Xylene<Scalar> >
, public Components::Gas<Scalar, Xylene<Scalar> >
{
    using Consts = Constants<Scalar>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    /*!
     * \brief A human readable name for the xylene
     */
    static std::string name()
    { return "xylene"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of xylene
     */
    constexpr static Scalar molarMass()
    { return 0.106; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of xylene
     */
    constexpr static Scalar criticalTemperature()
    { return 617.1; }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of xylene
     */
    constexpr static Scalar criticalPressure()
    { return 35.4e5; }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at xylene's boiling point (1 atm).
     */
    constexpr static Scalar boilingTemperature()
    { return 412.3; }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at xylene's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for xylene");
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at xylene's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for xylene");
    }

    /*!
     * \brief The saturation vapor pressure in \f$\mathrm{[Pa]}\f$ of pure xylene
     *        at a given temperature according to Antoine after Betz 1997 ->  Gmehling et al 1980 \cite gmehling1980 <BR>
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        const Scalar A = 7.00909;
        const Scalar B = 1462.266;
        const Scalar C = 215.110;

        Scalar T = temperature - 273.15;

        using std::pow;
        Scalar psat = 1.334*pow(10.0, (A - (B/(T + C))));  // in [mbar]
        psat *= 100.0;  // in [Pa] (0.001*1.E5)

        return psat;
    }

    /*!
     * \brief Specific heat cap of liquid xylene \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Missenard group contrib. method (chap 5-7, Table 5-11, s. example 5-8) \cite reid1987 <BR>
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidHeatCapacity(Scalar temp, Scalar pressure)
    {
        Scalar CH3,C6H5,H;
        // after Reid et al. : Missenard group contrib. method (s. example 5-8) \cite reid1987 <BR>
        // Xylene: C9H12  : 3* CH3 ; 1* C6H5 (phenyl-ring) ; -2* H (this was too much!)
        // linear interpolation between table values [J/(mol K)]

        if(temp < 298.0){                              // take care: extrapolation for Temp<273
            H = 13.4 + 1.2*(temp - 273.0)/25.0;        // 13.4 + 1.2 = 14.6 = H(T=298K) i.e. interpolation of table values 273<T<298
            CH3 = 40.0 + 1.6*(temp - 273.0)/25.0;    // 40 + 1.6 = 41.6 = CH3(T=298K)
            C6H5 = 113.0 + 4.2*(temp - 273.0)/25.0; // 113 + 4.2 = 117.2 = C6H5(T=298K)
        }
        else if(temp < 323.0){
            H = 14.6 + 0.9*(temp - 298.0)/25.0;        // i.e. interpolation of table values 298<T<323
            CH3 = 41.6 + 1.9*(temp - 298.0)/25.0;
            C6H5 = 117.2 + 6.2*(temp - 298.0)/25.0;
        }
        else if(temp < 348.0){
            H = 15.5 + 1.2*(temp - 323.0)/25.0;        // i.e. interpolation of table values 323<T<348
            CH3 = 43.5 + 2.3*(temp - 323.0)/25.0;
            C6H5 = 123.4 + 6.3*(temp - 323.0)/25.0;
        }
        else {
            H = 16.7 + 2.1*(temp - 348.0)/25.0;         // i.e. interpolation of table values 348<T<373
            CH3 = 45.8 + 2.5*(temp - 348.0)/25.0;        // take care: extrapolation for Temp>373
            C6H5 = 129.7 + 6.3*(temp - 348.0)/25.0;        // most likely leads to underestimation
        }

        return (C6H5 + 2*CH3 - H)/molarMass();// J/(mol K) -> J/(kg K)
    }


    /*!
     * \brief Specific enthalpy of liquid xylene \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidEnthalpy(const Scalar temperature,
                                 const Scalar pressure)
    {
        // Gauss quadrature rule:
        // Interval: [0K; temperature (K)]
        // Gauss-Legendre-Integration with variable transformation:
        // \int_a^b f(T) dT  \approx (b-a)/2 \sum_i=1^n \alpha_i f( (b-a)/2 x_i + (a+b)/2 )
        // with: n=2, legendre -> x_i = +/- \sqrt(1/3), \apha_i=1
        // here: a=273.15K, b=actual temperature in Kelvin
        // \leadsto h(T) = \int_273.15^T c_p(T) dT
        //              \approx 0.5 (T-273.15) * (cp( 0.5(temperature-273.15)sqrt(1/3) ) + cp(0.5(temperature-273.15)(-1)sqrt(1/3))

        // Enthalpy may have arbitrary reference state, but the empirical/fitted heatCapacity function needs Kelvin as input and is
        // fit over a certain temperature range. This suggests choosing an interval of integration being in the actual fit range.
        // I.e. choosing T=273.15K  as reference point for liquid enthalpy.
        using std::sqrt;
        const Scalar sqrt1over3 = sqrt(1./3.);
        // evaluation points according to Gauss-Legendre integration
        const Scalar TEval1 = 0.5*(temperature-273.15)*        sqrt1over3 + 0.5*(273.15+temperature);
        // evaluation points according to Gauss-Legendre integration
        const Scalar TEval2 = 0.5*(temperature-273.15)* (-1)*  sqrt1over3 + 0.5*(273.15+temperature);

        const Scalar h_n = 0.5 * (temperature-273.15) * ( liquidHeatCapacity(TEval1, pressure) + liquidHeatCapacity(TEval2, pressure) );

        return h_n;
    }

    /*!
     * \brief Latent heat of vaporization for xylene \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Chen method (chap. 7-11, Delta H_v = Delta H_v (T) according to chap. 7-12) \cite reid1987
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar heatVap(Scalar temperature,
                          const Scalar pressure)
    {
        using std::min;
        using std::max;
        temperature = min(temperature, criticalTemperature()); // regularization
        temperature = max(temperature, 0.0); // regularization

        constexpr Scalar T_crit = criticalTemperature();
        constexpr Scalar Tr1 = boilingTemperature()/criticalTemperature();
        constexpr Scalar p_crit = criticalPressure();

        // Chen method, eq. 7-11.4 (at boiling)
        using std::log;
        const Scalar DH_v_boil = Consts::R * T_crit * Tr1
                                  * (3.978 * Tr1 - 3.958 + 1.555*log(p_crit * 1e-5 /*Pa->bar*/ ) )
                                  / (1.07 - Tr1); /* [J/mol] */

        /* Variation with temp according to Watson relation eq 7-12.1*/
        using std::pow;
        const Scalar Tr2 = temperature/criticalTemperature();
        const Scalar n = 0.375;
        const Scalar DH_vap = DH_v_boil * pow(((1.0 - Tr2)/(1.0 - Tr1)), n);

        return (DH_vap/molarMass());          // we need [J/kg]
    }

    /*!
     * \brief Specific enthalpy of xylene vapor \f$\mathrm{[J/kg]}\f$.
     *
     *          This relation is true on the vapor pressure curve, i.e. as long
     *          as there is a liquid phase present.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        return liquidEnthalpy(temperature, pressure) + heatVap(temperature, pressure);
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of xylene gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        return IdealGas::density(molarMass(),
                                 temperature,
                                 pressure);
    }

    /*!
     * \brief The molar gas density \f$\mathrm{[mol/m^3]}\f$ of xylene gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }

    /*!
     * \brief The molar liquid density of pure xylene at a given pressure and temperature
     * \f$\mathrm{[mol/m^3]}\f$.
     *
     * source : Reid et al. (fourth edition): Modified Racket technique (chap. 3-11, eq. 3-11.9) \cite reid1987 <BR>
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidMolarDensity(Scalar temp, Scalar pressure)
    {
        // saturated molar volume according to Lide, CRC Handbook of
        // Thermophysical and Thermochemical Data, CRC Press, 1994
        // valid for 245 < Temp < 600
        using std::min;
        using std::max;
        temp = min(temp, 500.0); // regularization
        temp = max(temp, 250.0); // regularization

        using std::pow;
        const Scalar A1 = 0.25919;           // from table
        const Scalar A2 = 0.0014569;         // from table
        const Scalar expo = 1.0 + pow((1.0 - temp/criticalTemperature()), (2.0/7.0));
        const Scalar V = A2*pow(A1, expo);   // liquid molar volume [m^3/mol]

        return 1.0/V;             // molar density [mol/m^3]
    }

    /*!
     * \brief The density of pure xylene at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return liquidMolarDensity(temperature, pressure)*molarMass();
    }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
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
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of xylene vapor
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temp, Scalar pressure)
    {
        using std::min;
        using std::max;
        temp = min(temp, 500.0); // regularization
        temp = max(temp, 250.0); // regularization

        using std::pow;
        using std::exp;
        const Scalar Tr = max(temp/criticalTemperature(), 1e-10);
        const Scalar Fp0 = 1.0;
        const Scalar xi = 0.004623;
        const Scalar eta_xi = Fp0*(0.807*pow(Tr, 0.618)
                                   - 0.357*exp(-0.449*Tr)
                                   + 0.34*exp(-4.058*Tr)
                                   + 0.018);
        Scalar r = eta_xi/xi; // [1e-6 P]
        r /= 1.0e7; // [Pa s]

        return r;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure xylene.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temp, Scalar pressure)
    {
        using std::min;
        using std::max;
        temp = min(temp, 500.0); // regularization
        temp = max(temp, 250.0); // regularization

        const Scalar A = -3.82;
        const Scalar B = 1027.0;
        const Scalar C = -6.38e-4;
        const Scalar D = 4.52e-7;

        using std::exp;
        Scalar r = exp(A + B/temp + C*temp + D*temp*temp); // in [cP]
        r *= 1.0e-3; // in [Pa s]

        return r; // [Pa s]
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of xylene
     *
     * Thermal Conductivity of p-Xylene taken from the Dortmund Data Bank, see:
     * http://www.ddbst.de/en/EED/PCP/TCN_C176.php
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity( Scalar temperature,  Scalar pressure)
    {
        return 0.13;
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
