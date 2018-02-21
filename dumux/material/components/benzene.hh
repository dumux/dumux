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
 * \ingroup Components
 * \brief A simple benzene component (LNAPL).
 */
#ifndef DUMUX_BENZENE_HH
#define DUMUX_BENZENE_HH

#include <dumux/material/idealgas.hh>
#include "component.hh"

#include <dune/common/deprecated.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A simple benzene component (LNAPL).
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Benzene : public Component<Scalar, Benzene<Scalar> >
{

public:
    /*!
     * \brief A human readable name for the benzene
     */
    static std::string name()
    { return "benzene"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of benzene
     */
    static Scalar molarMass()
    {
        DUNE_THROW(Dune::NotImplemented, "molar mass for benzene");
    }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of benzene.
     */
    static Scalar criticalTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "criticalTemperature for benzene");
    }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of benzene.
     */
    static Scalar criticalPressure()
    {
        DUNE_THROW(Dune::NotImplemented, "criticalPressure for benzene");
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at benzene's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for benzene");
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at benzene's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for benzene");
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure benzene
     *        at a given temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    {
        DUNE_THROW(Dune::NotImplemented, "vaporPressure for benzene");
    }

    /*!
     * \brief Specific enthalpy of benzene steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "gasEnthalpy for benzene");
    }

    /*!
     * \brief Specific enthalpy of liquid benzene \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy for benzene");
    }

    /*!
     * \brief The density of benzene steam at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
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
     * \brief The density of pure benzene at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 889.51; // [kg/m^3]
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of benzene steam.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure, bool regularize=true)
    {
        DUNE_THROW(Dune::NotImplemented, "gasViscosity for benzene");
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure benzene.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 1.12e-3;//[Pa s]
    }
};

} // end namespace Components

template<class Scalar>
using Benzene DUNE_DEPRECATED_MSG("Now in the namespace: Components") = Dumux::Components::Benzene<Scalar>;

} // end namespace Dumux

#endif
