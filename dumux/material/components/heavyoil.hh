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
 * \ingroup Components
 *
 * \brief Properties of heavyoil.
 *
 */
#ifndef DUMUX_HEAVYOIL_HH
#define DUMUX_HEAVYOIL_HH

#include <dumux/material/idealgas.hh>
#include <dumux/material/components/component.hh>
#include <dumux/material/constants.hh>

namespace Dumux
{
/*!
 * \ingroup Components
 * \brief heavyoil
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class HeavyOil : public Component<Scalar, HeavyOil<Scalar> >
{
    typedef Dumux::Constants<Scalar> Consts;

public:
    /*!
     * \brief A human readable name for heavyoil
     */
    static const char *name()
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
        const Scalar A = 0.83;
        const Scalar B = 89.9513;
        const Scalar C = 139.6612;
        const Scalar D = 3.2033;
        const Scalar E = 1.0564;

        const Scalar mW = refComponentMolecularWeight() *1000. ;  // in [g/mol];

        return A+(B/mW)-(C/std::pow((mW+D),E));

    }

    static Scalar perbutationFactorBoilingTemperature()
    {
        const Scalar A = -7.4120e-2;    //All factors for 1 atm / 101325 pascals [760 mmHg]
        const Scalar B = -7.5041e-3;
        const Scalar C = -2.6031;
        const Scalar D = 9.0180e-2;
        const Scalar E = -1.0482;

        Scalar deltaSpecificGravity = std::log(refComponentSpecificGravity()/specificGravity());
        Scalar deltaMolecularWeight = std::log(refComponentMolecularWeight()/molecularWeight());

        return A*std::pow(deltaSpecificGravity,2) + B*deltaSpecificGravity + C*std::pow(deltaMolecularWeight,2) + D*deltaMolecularWeight
                + E*deltaSpecificGravity*deltaMolecularWeight;

    }

    static Scalar perbutationFactorCriticalTemperature()
    {
        const Scalar A = -6.1294e-2;
        const Scalar B = -7.0862e-2;
        const Scalar C = 6.1976e-1;
        const Scalar D = -5.7090e-2;
        const Scalar E = -8.4583e-2;

        Scalar deltaSpecificGravity = std::log(refComponentSpecificGravity()/specificGravity());
        Scalar deltaMolecularWeight = std::log(refComponentMolecularWeight()/molecularWeight());

        return A*std::pow(deltaSpecificGravity,2) + B*deltaSpecificGravity + C*std::pow(deltaMolecularWeight,2) + D*deltaMolecularWeight
                + E*deltaSpecificGravity*deltaMolecularWeight;

    }

    static Scalar perbutationFactorCriticalPressure()
    {
        const Scalar A = 1.8270e-1;
        const Scalar B = -2.4864e-1;
        const Scalar C = 8.3611;
        const Scalar D = -2.2389e-1;
        const Scalar E = 2.6984;

        Scalar deltaSpecificGravity = std::log(refComponentSpecificGravity()/specificGravity());
        Scalar deltaMolecularWeight = std::log(refComponentMolecularWeight()/molecularWeight());

        return A*std::pow(deltaSpecificGravity,2) + B*deltaSpecificGravity + C*std::pow(deltaMolecularWeight,2) + D*deltaMolecularWeight
                + E*deltaSpecificGravity*deltaMolecularWeight;

    }

     static Scalar refComponentBoilingTemperature()
    {
        const Scalar A = 477.63;    //All factors for 1 atm /  101325 pascals [760 mmHg]
        const Scalar B = 88.51;
        const Scalar C = 1007;
        const Scalar D = 1214.40;

        return A*std::log((1000.*refComponentMolecularWeight() + B)/(1000.*refComponentMolecularWeight()+C)) + D;

    }

        static Scalar refComponentCriticalTemperature()
    {
        const Scalar A = 226.50;
        const Scalar B = 6.78;
        const Scalar C = 1.282e6;
        const Scalar D = 2668;

        return A*std::log((1000.*refComponentMolecularWeight() + B)/(1000.*refComponentMolecularWeight()+C)) + D ;

    }

        static Scalar refComponentCriticalPressure()
    {
        const Scalar A = 141.20;
        const Scalar B = 45.66e-2;
        const Scalar C = 16.59e-3;
        const Scalar D = 2.19;

        return (A*1000.*molecularWeight())/(std::pow(B + (C*1000.*molecularWeight()),D)) ;

    }

   /*!
    * \brief Returns the temperature \f$\mathrm{[K]}\f$ at heavyoil's boiling point (1 atm)
    */
       static Scalar boilingTemperature()
    {

           return refComponentBoilingTemperature() * std::pow((1 + 2*perbutationFactorBoilingTemperature())/(1 - 2*perbutationFactorBoilingTemperature()),2);

    }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of heavyoil
     */
       static Scalar criticalTemperature()
    {

           return refComponentCriticalTemperature() * std::pow((1 + 2*perbutationFactorCriticalTemperature())/(1 - 2*perbutationFactorCriticalTemperature()),2);

    }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of heavyoil
     */
        static Scalar criticalPressure()
    {

            return refComponentCriticalPressure() * std::pow((1 + 2*perbutationFactorCriticalPressure())/(1 - 2*perbutationFactorCriticalPressure()),2);

    }

  /*!
     * \brief The saturation vapor pressure in \f$\mathrm{[Pa]}\f$ of
     *
     *
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        const Scalar A = 8.25990;
        const Scalar B = 2830.065;
        const Scalar C = 42.95101;

        Scalar T = temperature - 273.15;
        return 100*1.334*std::pow(10.0, (A - (B/(T + C))));  // in [Pa]

    }

    static Scalar vaporTemperature(Scalar pressure)
    {
        const Scalar A = 8.25990;
        const Scalar B = 2830.065;
        const Scalar C = 42.95101;

        const Scalar P = pressure;

        return  Scalar ((B/(A-std::log10(P/100*1.334)))-C);

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

        const Scalar sqrt1over3 = std::sqrt(1./3.);
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
        temperature = std::min(temperature, criticalTemperature()); // regularization
        temperature = std::max(temperature, 0.0); // regularization

        Scalar T_crit = criticalTemperature();
        Scalar Tr1 = boilingTemperature()/criticalTemperature();
        Scalar p_crit = criticalPressure();

        //        Chen method, eq. 7-11.4 (at boiling)
        const Scalar DH_v_boil = Consts::R * T_crit * Tr1
                                        * (3.978 * Tr1 - 3.958 + 1.555*std::log(p_crit * 1e-5 /*Pa->bar*/ ) )
                                        / (1.07 - Tr1); /* [J/mol] */

        /* Variation with temp according to Watson relation eq 7-12.1*/
        const Scalar Tr2 = temperature/criticalTemperature();
        const Scalar n = 0.375;
        const Scalar DH_vap = DH_v_boil * std::pow(((1.0 - Tr2)/(1.0 - Tr1)), n);

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
     * \brief
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        return IdealGas<Scalar>::density(molarMass(),
                                         temperature,
                                         pressure);
    }

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
        Scalar rhoReference = 906.; // [kg/m^3] at reference pressure and temperature
        Scalar compressCoeff = 1.e-8; // just a value without justification
        Scalar expansCoeff = 1.e-7; // also just a value
        Scalar rho = rhoReference * (1. + (pressure - 1.e5)*compressCoeff) * (1. - (temperature - 293.)*expansCoeff);

        return rho; // [kg/m^3]
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
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
        temperature = std::min(temperature, 500.0); // regularization
        temperature = std::max(temperature, 250.0);

        // reduced temperature
        Scalar Tr = temperature/criticalTemperature();

        Scalar Fp0 = 1.0;
        Scalar xi = 0.00474;
        Scalar eta_xi =
            Fp0*(0.807*std::pow(Tr,0.618)
                 - 0.357*std::exp(-0.449*Tr)
                 + 0.34*std::exp(-4.058*Tr)
                 + 0.018);

        return eta_xi/xi/1e7; // [Pa s]
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure heavyoil.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {

        /* according to Lashanizadegan et al (2008) in Chemical Engineering Communications:  */
        /* Simultaneous Heat and Fluid Flow in Porous Media: Case Study: Steam Injection for Tertiary Oil Recovery */

        //return 1027919.422*std::exp(-0.04862*temperature); // [Pa s]

        //according to http://www.ecltechnology.com/subsur/reports/pvt_tgb.pdf [Page 10]
        Scalar temperatureFahrenheit = (9/5)*(temperature-273.15)+32;
        Scalar API = 9;
            return ((std::pow(10,0.10231*std::pow(API,2)-3.9464*API+46.5037))*(std::pow(temperatureFahrenheit,-0.04542*std::pow(API,2)+1.70405*API-19.18)))*0.001;

    }
    /*!
     * \brief Specific heat cap of liquid heavyoil \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Missenard group contrib. method (chap 5-7, Table 5-11, s. example 5-8)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidHeatCapacity(const Scalar temperature,
                                     const Scalar pressure)
    {
        /* according to Lashanizadegan et al (2008) in Chemical Engineering Communications:  */
        /* Simultaneous Heat and Fluid Flow in Porous Media: Case Study: Steam Injection for Tertiary Oil Recovery */

        return 618.; // J/(kg K)
    }

protected:
    /*!
     * \brief The molar density of pure heavyoil at a given pressure and temperature
     * \f$\mathrm{[mol/m^3]}\f$.
     *
     * source : Reid et al. (fourth edition): Modified Racket technique (chap. 3-11, eq. 3-11.9)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar molarLiquidDensity_(Scalar temperature)
    {
        temperature = std::min(temperature, 500.0); // regularization
        temperature = std::max(temperature, 250.0);

        const Scalar Z_RA = 0.2556; // from equation
        const Scalar expo = 1.0 + std::pow(1.0 - temperature/criticalTemperature(), 2.0/7.0);
        Scalar V = Consts::R*criticalTemperature()/criticalPressure()*std::pow(Z_RA, expo); // liquid molar volume [cm^3/mol]

        return 1.0/V; // molar density [mol/m^3]
    }

};

} // end namespace

#endif
