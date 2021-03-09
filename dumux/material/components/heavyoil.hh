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
 * \brief Properties of the component heavyoil.
 */
#ifndef DUMUX_HEAVYOIL_HH
#define DUMUX_HEAVYOIL_HH

#include <algorithm>

#include <dune/common/math.hh>

#include <dumux/material/constants.hh>
#include <dumux/material/components/base.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/idealgas.hh>

namespace Dumux::Components {

/*!
 * \ingroup Components
 * \brief Properties of the component heavyoil
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class HeavyOil
: public Components::Base<Scalar, HeavyOil<Scalar> >
, public Components::Liquid<Scalar, HeavyOil<Scalar> >
, public Components::Gas<Scalar, HeavyOil<Scalar> >
{
    using Consts = Dumux::Constants<Scalar>;
    using IdealGas = Dumux::IdealGas<Scalar>;
public:
    /*!
     * \brief A human readable name for heavyoil
     */
    static std::string name()
    { return "heavyoil"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of heavyoil
     */
    constexpr static Scalar molarMass()
    { return .350; }

    /*!
     * \brief The MolecularWeight in \f$\mathrm{[kg/mol]}\f$ of refComponent
     */
    constexpr static Scalar refComponentMolecularWeight()
    { return .400; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of heavyoil
     */

    constexpr static Scalar molecularWeight()
    { return .350; }

    /*!
     * \brief The Specific Gravity \f$\mathrm{[  ]}\f$ of heavyoil
     */
    constexpr static Scalar specificGravity()
    { return 0.91; }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at heavyoil's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for heavyoil");
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at heavyoil's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for heavyoil");
    }

    static Scalar refComponentSpecificGravity()
    {
        constexpr Scalar A = 0.83;
        constexpr Scalar B = 89.9513;
        constexpr Scalar C = 139.6612;
        constexpr Scalar D = 3.2033;
        constexpr Scalar E = 1.0564;

        const Scalar mW = refComponentMolecularWeight() *1000. ;  // in [g/mol];

        using std::pow;
        return A+(B/mW)-(C/pow((mW+D),E));
    }

    static Scalar perbutationFactorBoilingTemperature()
    {
        constexpr Scalar A = -7.4120e-2;    //All factors for 1 atm / 101325 pascals [760 mmHg]
        constexpr Scalar B = -7.5041e-3;
        constexpr Scalar C = -2.6031;
        constexpr Scalar D = 9.0180e-2;
        constexpr Scalar E = -1.0482;

        using std::log;
        const Scalar deltaSpecificGravity = log(refComponentSpecificGravity()/specificGravity());
        const Scalar deltaMolecularWeight = log(refComponentMolecularWeight()/molecularWeight());

        using Dune::power;
        return A*power(deltaSpecificGravity,2) + B*deltaSpecificGravity + C*power(deltaMolecularWeight,2) + D*deltaMolecularWeight
                + E*deltaSpecificGravity*deltaMolecularWeight;
    }

    static Scalar perbutationFactorCriticalTemperature()
    {
        constexpr Scalar A = -6.1294e-2;
        constexpr Scalar B = -7.0862e-2;
        constexpr Scalar C = 6.1976e-1;
        constexpr Scalar D = -5.7090e-2;
        constexpr Scalar E = -8.4583e-2;

        using std::log;
        const Scalar deltaSpecificGravity = log(refComponentSpecificGravity()/specificGravity());
        const Scalar deltaMolecularWeight = log(refComponentMolecularWeight()/molecularWeight());

        using Dune::power;
        return A*power(deltaSpecificGravity,2) + B*deltaSpecificGravity + C*power(deltaMolecularWeight,2) + D*deltaMolecularWeight
                + E*deltaSpecificGravity*deltaMolecularWeight;
    }

    static Scalar perbutationFactorCriticalPressure()
    {
        constexpr Scalar A = 1.8270e-1;
        constexpr Scalar B = -2.4864e-1;
        constexpr Scalar C = 8.3611;
        constexpr Scalar D = -2.2389e-1;
        constexpr Scalar E = 2.6984;

        using std::log;
        const Scalar deltaSpecificGravity = log(refComponentSpecificGravity()/specificGravity());
        const Scalar deltaMolecularWeight = log(refComponentMolecularWeight()/molecularWeight());

        using Dune::power;
        return A*power(deltaSpecificGravity,2) + B*deltaSpecificGravity + C*power(deltaMolecularWeight,2) + D*deltaMolecularWeight
                + E*deltaSpecificGravity*deltaMolecularWeight;
    }

    static Scalar refComponentBoilingTemperature()
    {
        constexpr Scalar A = 477.63;    //All factors for 1 atm /  101325 pascals [760 mmHg]
        constexpr Scalar B = 88.51;
        constexpr Scalar C = 1007;
        constexpr Scalar D = 1214.40;

        using std::log;
        return A*log((1000.*refComponentMolecularWeight() + B)/(1000.*refComponentMolecularWeight()+C)) + D;
    }

    static Scalar refComponentCriticalTemperature()
    {
        constexpr Scalar A = 226.50;
        constexpr Scalar B = 6.78;
        constexpr Scalar C = 1.282e6;
        constexpr Scalar D = 2668;

        using std::log;
        return A*log((1000.*refComponentMolecularWeight() + B)/(1000.*refComponentMolecularWeight()+C)) + D ;
    }

    static Scalar refComponentCriticalPressure()
    {
        constexpr Scalar A = 141.20;
        constexpr Scalar B = 45.66e-2;
        constexpr Scalar C = 16.59e-3;
        constexpr Scalar D = 2.19;

        using std::pow;
        return (A*1000.*molecularWeight())/(pow(B + (C*1000.*molecularWeight()),D)) ;
    }

   /*!
    * \brief Returns the temperature \f$\mathrm{[K]}\f$ at heavyoil's boiling point (1 atm)
    */
    static Scalar boilingTemperature()
    {
        using Dune::power;
        return refComponentBoilingTemperature() * power((1 + 2*perbutationFactorBoilingTemperature())/(1 - 2*perbutationFactorBoilingTemperature()),2);
    }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of heavyoil
     */
    static Scalar criticalTemperature()
    {
        using Dune::power;
        return refComponentCriticalTemperature() * power((1 + 2*perbutationFactorCriticalTemperature())/(1 - 2*perbutationFactorCriticalTemperature()),2);
    }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of heavyoil
     */
    static Scalar criticalPressure()
    {
        using Dune::power;
        return refComponentCriticalPressure() * power((1 + 2*perbutationFactorCriticalPressure())/(1 - 2*perbutationFactorCriticalPressure()),2);
    }

    /*!
     * \brief The saturation vapor pressure in \f$\mathrm{[Pa]}\f$ of
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        constexpr Scalar A = 8.25990;
        constexpr Scalar B = 2830.065;
        constexpr Scalar C = 42.95101;

        const Scalar T = temperature - 273.15;

        using std::pow;
        return 100*1.334*pow(10.0, (A - (B/(T + C))));  // in [Pa]
    }

    static Scalar vaporTemperature(Scalar pressure)
    {
        constexpr Scalar A = 8.25990;
        constexpr Scalar B = 2830.065;
        constexpr Scalar C = 42.95101;

        const Scalar P = pressure;

        using std::log10;
        return  Scalar ((B/(A-log10(P/100*1.334)))-C);
    }

    /*!
     * \brief Specific enthalpy of liquid heavyoil \f$\mathrm{[J/kg]}\f$.
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
        const Scalar TEval1 = 0.5*(temperature-273.15)*        sqrt1over3 + 0.5*(273.15+temperature)  ; // evaluation points according to Gauss-Legendre integration
        const Scalar TEval2 = 0.5*(temperature-273.15)* (-1)*  sqrt1over3 + 0.5*(273.15+temperature)  ; // evaluation points according to Gauss-Legendre integration

        const Scalar h_n = 0.5 * (temperature-273.15) * ( liquidHeatCapacity(TEval1, pressure) + liquidHeatCapacity(TEval2, pressure) ) ;

        return h_n;
    }

    /*!
     * \brief Latent heat of vaporization for heavyoil \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Chen method (chap. 7-11, Delta H_v = Delta H_v (T) according to chap. 7-12)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar heatVap(Scalar temperature,
                   const  Scalar pressure)
    {
        using std::clamp;
        temperature = clamp(temperature, 0.0, criticalTemperature()); // regularization

        const Scalar T_crit = criticalTemperature();
        const Scalar Tr1 = boilingTemperature()/criticalTemperature();
        const Scalar p_crit = criticalPressure();

        //        Chen method, eq. 7-11.4 (at boiling)
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
     * \brief Specific enthalpy of heavyoil vapor \f$\mathrm{[J/kg]}\f$.
     *
     *      This relation is true on the vapor pressure curve, i.e. as long
     *      as there is a liquid phase present.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        return liquidEnthalpy(temperature,pressure) + heatVap(temperature, pressure);
    }

    /*!
     * \brief The (ideal) gas density of heavyoil vapor at a given temperature and pressure \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     * \brief The molar density of pure heavyoil in \f$\mathrm{[mol/m^3]}\f$,
     *   depending on pressure and temperature.
     * \param temperature The temperature of the gas
     * \param pressure The pressure of the gas
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }

    /*!
     * \brief The density of pure heavyoil at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {

        /* according to Lashanizadegan et al (2008) in Chemical Engineering Communications:  */
        /* Simultaneous Heat and Fluid Flow in Porous Media: Case Study: Steam Injection for Tertiary Oil Recovery */
        constexpr Scalar rhoReference = 906.; // [kg/m^3] at reference pressure and temperature
        constexpr Scalar compressCoeff = 1.e-8; // just a value without justification
        constexpr Scalar expansCoeff = 1.e-7; // also just a value
        return rhoReference * (1. + (pressure - 1.e5)*compressCoeff) * (1. - (temperature - 293.)*expansCoeff);
    }

    /*!
     * \brief The molar density of pure heavyoil in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

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
    { return true; }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of heavyoil vapor
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure, bool regularize=true)
    {
        using std::clamp;
        temperature = clamp(temperature, 250.0, 500.0); // regularization

        // reduced temperature
        const Scalar Tr = temperature/criticalTemperature();

        using std::pow;
        using std::exp;
        constexpr Scalar Fp0 = 1.0;
        constexpr Scalar xi = 0.00474;
        const Scalar eta_xi =
            Fp0*(0.807*pow(Tr,0.618)
                 - 0.357*exp(-0.449*Tr)
                 + 0.34*exp(-4.058*Tr)
                 + 0.018);

        return eta_xi/xi/1e7; // [Pa s]
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure heavyoil.
     *
     * Lashanizadegan et al. (2008) \cite lashanizadegan2008 <BR>
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {

        /* according to Lashanizadegan et al (2008) in Chemical Engineering Communications:  */
        /* Simultaneous Heat and Fluid Flow in Porous Media: Case Study: Steam Injection for Tertiary Oil Recovery */

        //using std::exp;
        //return 1027919.422*exp(-0.04862*temperature); // [Pa s]

        //according to http://www.ecltechnology.com/subsur/reports/pvt_tgb.pdf [Page 10]
        const Scalar temperatureFahrenheit = (9/5)*(temperature-273.15)+32;
        constexpr Scalar API = 9;

        using std::pow;
        using Dune::power;
        return ((pow(10,0.10231*power(API,2)-3.9464*API+46.5037))*(pow(temperatureFahrenheit,-0.04542*power(API,2)+1.70405*API-19.18)))*0.001;

    }
    /*!
     * \brief Specific heat cap of liquid heavyoil \f$\mathrm{[J/kg]}\f$.
     *
     * Lashanizadegan et al. (2008) \cite lashanizadegan2008 <BR>
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidHeatCapacity(const Scalar temperature,
                                     const Scalar pressure)
    {
        return 618.; // J/(kg K)
    }

     /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of heavy oil
     *
     * Lashanizadegan et al. (2008) \cite lashanizadegan2008 <BR>
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity( Scalar temperature,  Scalar pressure)
    {
        return 0.127;
    }
};

} // end namespace Dumux::Components

#endif
