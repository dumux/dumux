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
 * \brief Properties of nDodecane.
 */
#ifndef DUMUX_nDodecane_HH
#define DUMUX_nDodecane_HH

#include <cmath>
#include <dumux/material/idealgas.hh>
#include <dumux/material/components/component.hh>
#include <dumux/material/constants.hh>


namespace Dumux
{
/*!
 * \ingroup Components
 * \brief nDodecane
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class nDodecane : public Component<Scalar, nDodecane<Scalar> >
{
    typedef Constants<Scalar> Consts;

public:
    /*!
     * \brief A human readable name for the nDodecane
     */
    static const char *name()
    { return "nDodecane"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of nDodecane
     */
    constexpr static Scalar molarMass()
    { return 0.17; }     //adil

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of nDodecane
     */
    constexpr static Scalar criticalTemperature()
    { return DUNE_THROW(Dune::NotImplemented, "criticalTemperature for nDodecane"); //adil }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of nDodecane
     */
    constexpr static Scalar criticalPressure()
    { return DUNE_THROW(Dune::NotImplemented, "criticalPressure for nDodecane"); //adil }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at nDodecane's boiling point (1 atm).
     */
    constexpr static Scalar boilingTemperature()
    { return 487; //adil http://en.wikipedia.org/wiki/Dodecane }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at nDodecane's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for nDodecane");
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at nDodecane's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for nDodecane");
    }

    /*!
     * \brief The saturation vapor pressure in \f$\mathrm{[Pa]}\f$ of pure nDodecane
     *        at a given temperature according to Antoine after Betz 1997 ->  Gmehling et al 1980
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */

    static Scalar vaporPressure(Scalar temperature) //adil
    {
        const Scalar A = 6.99765;
        const Scalar B = 1639.27;
        const Scalar C = 181.835;

        Scalar T = temperature - 273.15;
        Scalar psat = 1.334*std::pow(10.0, (A - (B/(T + C))));  // in [mbar]
        psat *= 100.0;  // in [Pa] (0.001*1.E5)

        return psat;
    }

    /*!
     * \brief Specific heat cap of liquid nDodecane \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Missenard group contrib. method (chap 5-7, Table 5-11, s. example 5-8)
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidHeatCapacity()
    {
        return 375.97/0.17;        //adil [J/kgK]
    }
    
        static Scalar gasHeatCapacity()
    {
        return 2732;        //adil [J/kgK]
    }


    /*!
     * \brief Specific enthalpy of liquid nDodecane \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidEnthalpy(const Scalar temperature,
                                 const Scalar pressure)
    {
      return DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy for nDodecane"); //adil http://en.wikipedia.org/wiki/Heptane;
    }

    /*!
     * \brief Latent heat of vaporization for nDodecane \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Chen method (chap. 7-11, Delta H_v = Delta H_v (T) according to chap. 7-12)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar heatVap()
    {
        return  323000;     //adil
    }
    
    /*Heat Conductivity [W/mK]*/
    
        static Scalar heatConductivityVap()
    {
        return 0.025;     //adil
    }


    /*Heat Conductivity [W/mK]*/
    
        static Scalar heatConductivityLiquid()
    {
        return 0.095;     //adil
    }

    /*!
     * \brief Specific enthalpy of nDodecane vapor \f$\mathrm{[J/kg]}\f$.
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
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of nDodecane gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        return 1/0.23;  //adil
    }

    /*!
     * \brief The density \f$\mathrm{[mol/m^3]}\f$ of nDodecane gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
    */
    static Scalar molarGasDensity(Scalar temperature, Scalar pressure)
    {
        return (gasDensity(temperature, pressure) / molarMass());
    }

    /*!
     * \brief The molar density of pure nDodecane at a given pressure and temperature
     * \f$\mathrm{[mol/m^3]}\f$.
     *
     * source : Reid et al. (fourth edition): Modified Racket technique (chap. 3-11, eq. 3-11.9)
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar molarLiquidDensity(Scalar temp, Scalar pressure)
    {  return DUNE_THROW(Dune::NotImplemented, "molarLiquidDensity for nDodecane"); //adil             // molar density [mol/m^3]
    }

    /*!
     * \brief The density of pure nDodecane at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity()
    {
        return 750; // [kg/m^3]     //adil
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
    { return false; }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of nDodecane vapor
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    static Scalar gasViscosity(Scalar temp, Scalar pressure, bool regularize=true)
    {
        return 7.4e-6 * (1/0.23); // kinem. viscosity [m^2/s]*gasDensity [kg/m^3]= Dynamic viscosity[Pa s]
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure nDodecane.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temp, Scalar pressure)
    {

        return 2.027e-6 * 750; // kinem. viscosity [m^2/s]*liquidDensity [kg/m^3]= Dynamic viscosity[Pa s]
    }
};

} // end namespace

#endif
