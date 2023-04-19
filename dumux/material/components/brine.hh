// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief A class for the brine fluid properties,.
 */
#ifndef DUMUX_BRINE_HH
#define DUMUX_BRINE_HH

#include <cmath>

#include <dune/common/math.hh>

#include <dumux/common/parameters.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/nacl.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux::Components {

/*!
 * \ingroup Components
 * \brief A class for the brine fluid properties.
 *
 * \tparam Scalar The type used for scalar values
 * \tparam H2O Static polymorphism: the Brine class can access all properties of the H2O class
 * \note This is an implementation of brine as a pseudo-component with a constant salinity.
 * \note the salinity is read from the input file and is a mandatory parameter
 */
template <class Scalar,
          class H2O_Tabulated = Components::TabulatedComponent<Components::H2O<Scalar>>>
class Brine
: public Components::Base<Scalar, Brine<Scalar, H2O_Tabulated> >
, public Components::Liquid<Scalar, Brine<Scalar, H2O_Tabulated> >
, public Components::Gas<Scalar, Brine<Scalar, H2O_Tabulated> >
{
    using ThisType = Brine<Scalar, H2O_Tabulated>;
public:
    using H2O = Components::TabulatedComponent<Dumux::Components::H2O<Scalar>>;

    //! The ideal gas constant \f$\mathrm{[J/mol/K]}\f$
    static constexpr Scalar R = Constants<Scalar>::R;

    /*!
     * \brief A human readable name for the brine.
     */
    static std::string name()
    { return "Brine"; }

    /*!
     * \brief Return the constant salinity
     */
    static Scalar salinity()
    {
        static const Scalar salinity = getParam<Scalar>("Brine.Salinity");
        return salinity;
    }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of brine.
     * This assumes that the salt is pure NaCl.
     */
   static Scalar molarMass()
   {
       const Scalar M1 = H2O::molarMass();
       const Scalar M2 = Components::NaCl<Scalar>::molarMass(); // molar mass of NaCl [kg/mol]
       return M1*M2/(M2 + ThisType::salinity()*(M1 - M2));
   };

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of brine. Here, it is assumed to be equal to that of pure water.
     * The critical temperature of brine is dependent on the salt concentration.
     * With increasing NaCl mass fraction the critical point shifts to higher temperatures, see Bakker, 2018 (DOI: 10.1016/j.cageo.2018.03.003).
     */
    static Scalar criticalTemperature()
    { return H2O::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of brine. Here, it is assumed to be equal to that of pure water.
     * The critical pressure of brine is dependent on the salt concentration.
     * With increasing NaCl mass fraction the critical point shifts to higher pressures, see Bakker, 2018 (DOI: 10.1016/j.cageo.2018.03.003).
     */
    static Scalar criticalPressure()
    { return H2O::criticalPressure(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at brine's triple point. Here, it is assumed to be equal to that of pure water.
     * The triple temperature of brine is dependent on the salt concentration, see Bakker, 2018 (DOI: 10.1016/j.cageo.2018.03.003).
     */
    static Scalar tripleTemperature()
    { return H2O::tripleTemperature(); }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at brine's triple point. Here, it is assumed to be equal to that of pure water.
     * The triple pressure of brine is dependent on the salt concentration, see Bakker, 2018 (DOI: 10.1016/j.cageo.2018.03.003).
     */
    static Scalar triplePressure()
    { return H2O::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure brine
     *        at a given temperature.
     * \note The vapor pressure of brine decreases with the mole fraction of water in the liquid phase.
     * This is described by Raoult's law, see Thomas Fetzer's Dissertation Eq. 2.11.
     * It is also the simplified version of the Kelvin equation, neglecting the influence of the capillary pressure here.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        //calculate mole fraction
        const Scalar M1 = H2O::molarMass();
        const Scalar M2 = Components::NaCl<Scalar>::molarMass(); // molar mass of NaCl [kg/mol]
        const Scalar xNaClLiquid = - M1 * ThisType::salinity() / ((M2-M1) * ThisType::salinity() - M2);

        // Raoult's law, see Thomas Fetzer's Dissertation Eq. 2.11.
        return H2O::vaporPressure(temperature) * (1 - xNaClLiquid);
    }

    /*!
     * \brief Specific enthalpy of gaseous brine \f$\mathrm{[J/kg]}\f$.
     * Only water volatile and salt is suppose to stay in the liquid phase.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    { return H2O::gasEnthalpy(temperature, pressure); }

    /*!
     * \brief Specific enthalpy of liquid brine \f$\mathrm{[J/kg]}\f$.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param p pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * Equations given in:
     * - Palliser & McKibbin (1998) \cite palliser1998 <BR>
     * - Michaelides (1981) \cite michaelides1981 <BR>
     * - Daubert & Danner (1989) \cite daubert1989
     *
     */
    static const Scalar liquidEnthalpy(Scalar T, Scalar p)
    {
        /*Numerical coefficients from PALLISER*/
        static const Scalar f[] = {
            2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10
        };

        /*Numerical coefficients from MICHAELIDES for the enthalpy of brine*/
        static const Scalar a[4][3] = {
            { +9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }
        };

        const Scalar theta = T - 273.15;
        const Scalar salSat = f[0] + f[1]*theta + f[2]*theta*theta + f[3]*theta*theta*theta;

        /*Regularization*/
        using std::min;
        using std::max;
        const Scalar salinity = min(max(ThisType::salinity(),0.0), salSat);

        const Scalar hw = H2O::liquidEnthalpy(T, p)/1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/const Scalar h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
                        +((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        const Scalar m = (1E3/58.44)*(salinity/(1-salinity));

        using Dune::power;
        Scalar d_h = 0;
        for (int i = 0; i<=3; i++) {
            for (int j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * power(theta, i) * power(m, j);
            }
        }

        /* heat of dissolution for halite according to Michaelides 1971 */
        const Scalar delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without any dissolved gas */
        const Scalar h_ls1 =(1-salinity)*hw + salinity*h_NaCl + salinity*delta_h; /* kJ/kg */
        return h_ls1*1E3; /*J/kg*/
    }

    /*!
     * \brief Specific isobaric heat capacity of brine \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf  \cite IAPWS1997
     */
    static const Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        const Scalar eps = temperature*1e-8;
        return (liquidEnthalpy(temperature + eps, pressure)- liquidEnthalpy(temperature, pressure))/eps;
    }

    /*!
     * \brief Specific isobaric heat capacity of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf  \cite IAPWS1997
     */
    static const Scalar gasHeatCapacity(Scalar temperature,
                                        Scalar pressure)
    {
        return H2O::gasHeatCapacity(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return H2O::gasInternalEnergy(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of liquid brine \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    {
        return liquidEnthalpy(temperature, pressure) - pressure/liquidDensity(temperature, pressure);
    }

    /*!
     * \brief The density of steam at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { return H2O::gasDensity(temperature, pressure); }

    /*!
     *  \brief The molar density of steam in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *  We take the value of the H2O gas molar density here because salt is not in the gas phase.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return H2O::gasMolarDensity(temperature, pressure); }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return H2O::gasIsIdeal(); }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return H2O::gasIsCompressible(); }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return H2O::liquidIsCompressible(); }

    /*!
     * \brief The density of pure brine at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * \note The density is computed as a function of the salt mass fraction, pressure and temperature.
     * The used function is an empirical relationship fitted to experimental data.
     * It is presented by Batzle and Wang, 1992 (DOI: 10.1190/1.1443207) \cite batzle1992,
     * better description and comparison with other approaches in Adams and Bachu, 2002
     * (DOI: 10.1046/j.1468-8123.2002.00041.x) \cite adams2002.
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        using std::max;
        const Scalar TempC = temperature - 273.15;
        const Scalar pMPa = pressure/1.0E6;
        const Scalar salinity = max(0.0, ThisType::salinity());

        const Scalar rhow = H2O::liquidDensity(temperature, pressure);

        const Scalar  density = rhow +
                                1000*salinity*(
                                    0.668 +
                                    0.44*salinity +
                                    1.0E-6*(
                                        300*pMPa -
                                        2400*pMPa*salinity +
                                        TempC*(
                                            80.0 +
                                            3*TempC -
                                            3300*salinity -
                                            13*pMPa +
                                            47*pMPa*salinity)));
        assert(density > 0.0);
        return density;
    }

    /*!
     * \brief The molar density of brine in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
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
   { return H2O::gasPressure(temperature, density); }

   /*!
    * \brief The pressure of brine in \f$\mathrm{[Pa]}\f$ at a given density and
    *        temperature.
    *
    * \param temperature temperature of component in \f$\mathrm{[K]}\f$
    * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
    */
   static Scalar liquidPressure(Scalar temperature, Scalar density)
   {
       // We use the Newton method for this. For the initial value we
       // assume the pressure to be 10% higher than the vapor
       // pressure
       Scalar pressure = 1.1*vaporPressure(temperature);
       const Scalar eps = pressure*1e-7;

       Scalar deltaP = pressure*2;

       using std::abs;
       for (int i = 0; i < 5 && abs(pressure*1e-9) < abs(deltaP); ++i) {
           Scalar f = liquidDensity(temperature, pressure) - density;

           Scalar df_dp;
           df_dp = liquidDensity(temperature, pressure + eps);
           df_dp -= liquidDensity(temperature, pressure - eps);
           df_dp /= 2*eps;

           deltaP = - f/df_dp;

           pressure += deltaP;
       }
       assert(pressure > 0.0);
       return pressure;
   }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param temperature temperature of component
     * \param pressure pressure of component
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    { return H2O::gasViscosity(temperature, pressure); };

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure brine.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * \note The viscosity is computed as a function of the salt mass fraction and temperature.
     * The used function is an empirical relationship fitted to experimental data.
     * It is presented by Batzle and Wang, 1992 (DOI: 10.1190/1.1443207)  \cite batzle1992,
     * better description and comparison with other approaches in Adams and Bachu, 2002 (DOI: 10.1046/j.1468-8123.2002.00041.x) \cite adams2002.
     * However, the equation given in Adams and Bachu, 2002(DOI: 10.1046/j.1468-8123.2002.00041.x) \cite adams2002
     * is obviously wrong when compared to the original by Batzle and Wang, 1992 (DOI: 10.1190/1.1443207)  \cite batzle1992.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        // regularisation
        using std::max;
        temperature = max(temperature, 275.0);
        const Scalar salinity = max(0.0, ThisType::salinity());

        using std::pow;
        using Dune::power;
        using std::exp;
        const Scalar T_C = temperature - 273.15;
        const Scalar A = (0.42*power((pow(salinity, 0.8)-0.17), 2) + 0.045)*pow(T_C, 0.8);
        const Scalar mu_brine = 0.1 + 0.333*salinity + (1.65+91.9*salinity*salinity*salinity)*exp(-A); //[cP]
        assert(mu_brine > 0.0);
        return mu_brine/1000.0; //[Pa s]
    }

    /*!
     * \brief Thermal conductivity of a brine \f$\mathrm{[W/(m K)]}\f$.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * The thermal conductivity of brine is implemented based on the contribution of NaCl (\f$\lambda_{brine}\f$/\f$\lambda_{H_2O}\f$) of \cite Yusufova1975 https://link.springer.com/content/pdf/10.1007/BF00867119.pdf, also discussed in \cite Ozbek1980 https://docecity.com/thermal-conductivity-of-aqueous-sodium-chloride-acs-publicat-5f10766acba00.html
     */
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    {
            Scalar tempC = temperature-273.15;
            Scalar xNaCl = ThisType::salinity() * H2O::molarMass() / (ThisType::salinity() * H2O::molarMass() + (1-ThisType::salinity() )*Components::NaCl<Scalar>::molarMass()); // mole fraction of NaCl
            Scalar m = xNaCl/(H2O::molarMass()*(1- xNaCl)); // molality of NaCl
            Scalar S = 5844.3 * m / (1000 + 58.443 *m);
            Scalar contribNaClFactor = 1.0 - (2.3434e-3 - 7.924e-6*tempC + 3.924e-8*tempC*tempC)*S + (1.06e-5 - 2.0e-8*tempC + 1.2e-10*tempC*tempC)*S*S;
            return contribNaClFactor * H2O::liquidThermalConductivity(temperature, pressure);
    }
};

template <class Scalar, class H2O>
struct IsAqueous<Brine<Scalar, H2O>> : public std::true_type {};

} // end namespace Dumux::Components

#endif
