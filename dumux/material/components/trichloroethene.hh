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
 * \brief A simple implementation of Trichloroethene (TCE), a DNAPL.
 */
#ifndef DUMUX_TRICHLOROETHENE_HH
#define DUMUX_TRICHLOROETHENE_HH

#include <dumux/material/idealgas.hh>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A simple implementation of TCE as exemplary component for a dense NAPL.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Trichloroethene
: public Components::Base<Scalar, Trichloroethene<Scalar> >
, public Components::Liquid<Scalar, Trichloroethene<Scalar> >
, public Components::Gas<Scalar, Trichloroethene<Scalar> >
{
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for the dense NAPL TCE.
     */
    static std::string name()
    { return "Trichloroethene"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of TCE.
     */
    static constexpr Scalar molarMass()
    {
        return 131.39e-3; // [kg/mol]
    }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of TCE.
     */
    static Scalar criticalTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "criticalTemperature for TCE");
    }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of TCE.
     */
    static Scalar criticalPressure()
    {
        DUNE_THROW(Dune::NotImplemented, "criticalPressure for TCE");
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at TCE's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for TCE");
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at TCE's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for TCE");
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure TCE
     *        at a given temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    {
        return 3900; // [Pa] (at 20C)
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
     * \brief Returns true if the liquid phase viscostiy is constant
     */
    static constexpr bool liquidViscosityIsConstant()
    { return true; }

    /*!
     * \brief The density of steam at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
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
     * \brief The density of pure TCE at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1460.0; // [kg/m^3]
    }

    /*!
     * \brief The molar density of pure TCE in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure TCE.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 5.7e-4;// [Pa*s]
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
