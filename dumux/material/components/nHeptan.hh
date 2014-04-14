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
 * \brief Properties of nHeptan.
 */
#ifndef DUMUX_nHeptan_HH
#define DUMUX_nHeptan_HH

#include <cmath>
#include <dumux/material/idealgas.hh>
#include <dumux/material/components/component.hh>
#include <dumux/material/constants.hh>


namespace Dumux
{
/*!
 * \ingroup Components
 * \brief nHeptan
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class nHeptan : public Component<Scalar, nHeptan<Scalar> >
{
    typedef Constants<Scalar> Consts;

public:
    /*!
     * \brief A human readable name for the nHeptan
     */
    static const char *name()
    { return "nHeptan"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of nHeptan
     */
    constexpr static Scalar molarMass()
    { return 0.1002; }     //adil

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of nHeptan
     */
    constexpr static Scalar criticalTemperature()
    { return DUNE_THROW(Dune::NotImplemented, "criticalTemperature for nHeptan"); //adil }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of nHeptan
     */
    constexpr static Scalar criticalPressure()
    { return DUNE_THROW(Dune::NotImplemented, "criticalPressure for nHeptan"); //adil }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at nHeptan's boiling point (1 atm).
     */
    constexpr static Scalar boilingTemperature()
    { return 371.2; //adil http://en.wikipedia.org/wiki/Heptane }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at nHeptan's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for nHeptan");
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at nHeptan's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for nHeptan");
    }

    /*!
     * \brief The saturation vapor pressure in \f$\mathrm{[Pa]}\f$ of pure nHeptan
     *        at a given temperature according to Antoine after Betz 1997 ->  Gmehling et al 1980
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */

    static Scalar vaporPressure(Scalar temperature) //adil
    {
        const Scalar A = 6.89385;
        const Scalar B = 1264.37;
        const Scalar C = 216.636;

        Scalar T = temperature - 273.15;
        Scalar psat = 1.334*std::pow(10.0, (A - (B/(T + C))));  // in [mbar]
        psat *= 100.0;  // in [Pa] (0.001*1.E5)

        return psat;
    }

    /*!
     * \brief Specific heat cap of liquid nHeptan \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Missenard group contrib. method (chap 5-7, Table 5-11, s. example 5-8)
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidHeatCapacity()
    {
        return 224.93/0.100;        //adil
    }
    
        static Scalar gasHeatCapacity()
    {
        return 2772;        //adil
    }


    /*!
     * \brief Specific enthalpy of liquid nHeptan \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidEnthalpy(const Scalar temperature,
                                 const Scalar pressure)
    {
      return DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy for nHeptan"); //adil http://en.wikipedia.org/wiki/Heptane;
    }

    /*!
     * \brief Latent heat of vaporization for nHeptan \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Chen method (chap. 7-11, Delta H_v = Delta H_v (T) according to chap. 7-12)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar heatVap()
    {
        return  316600;     //adil
    }
    
    /*Heat Conductivity [W/mK]*/
    
        static Scalar heatConductivityVap()
    {
        return 0.018;     //adil
    }


    /*Heat Conductivity [W/mK]*/
    
        static Scalar heatConductivityLiquid()
    {
        return 0.113;     //adil
    }

    /*!
     * \brief Specific enthalpy of nHeptan vapor \f$\mathrm{[J/kg]}\f$.
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
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of nHeptan gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        return 1/0.29;  //adil
    }

    /*!
     * \brief The density \f$\mathrm{[mol/m^3]}\f$ of nHeptan gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
    */
    static Scalar molarGasDensity(Scalar temperature, Scalar pressure)
    {
        return (gasDensity(temperature, pressure) / molarMass());
    }

    /*!
     * \brief The molar density of pure nHeptan at a given pressure and temperature
     * \f$\mathrm{[mol/m^3]}\f$.
     *
     * source : Reid et al. (fourth edition): Modified Racket technique (chap. 3-11, eq. 3-11.9)
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar molarLiquidDensity(Scalar temp, Scalar pressure)
    {  return DUNE_THROW(Dune::NotImplemented, "molarLiquidDensity for nHeptan"); //adil             // molar density [mol/m^3]
    }

    /*!
     * \brief The density of pure nHeptan at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity()
    {
        return 680; // [kg/m^3]     //adil
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
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of nHeptan vapor
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    static Scalar gasViscosity(Scalar temp, Scalar pressure, bool regularize=true)
    {
        return 7.584e-6 * 3.448; // kinem. viscosity [m^2/s]*gasDensity [kg/m^3]= Dynamic viscosity[Pa s]
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure nHeptan.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temp, Scalar pressure)
    {

        return 2.694e-7 * 680; // kinem. viscosity [m^2/s]*liquidDensity [kg/m^3]= Dynamic viscosity[Pa s]
    }
};

} // end namespace

#endif
