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
 * \brief A ficitious component to be implemented in exercise 3.
 */
#ifndef DUMUX_MYCOMPRESSIBLECOMPONENT_HH
#define DUMUX_MYCOMPRESSIBLECOMPONENT_HH

#ifndef HEADERCHECK

#include <dumux/material/idealgas.hh>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>

#include "myincompressiblecomponent.hh"

namespace Dumux
{
/*!
 * \ingroup Components
 * \brief A ficitious component to be implemented in exercise 3.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class MyCompressibleComponent
: public Components::Base<Scalar, MyIncompressibleComponent<Scalar> >
, public Components::Liquid<Scalar, MyIncompressibleComponent<Scalar> >
{
public:
    /*!
     * \brief A human readable name for MyCompressibleComponent.
     */
    static std::string name()
    { return "MyCompressibleComponent"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of MyCompressibleComponent.
     */
    static Scalar molarMass()
    {
        return 131.39e-3; // [kg/mol]
    }

    /*!
     * \brief The density of MyCompressibleComponent at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        static const Scalar rho_min = 1440;
        static const Scalar rho_max = 1480;
        static const Scalar k = 5e-7;

        using std::exp;
        return rho_min + (rho_max - rho_min)/(1 + rho_min*exp(-1.0*k*(rho_max - rho_min)*pressure)); // [kg/m^3]
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of MyCompressibleComponent.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 5.7e-4;// [Pa*s]
    }

    /*****************************************************************
     * The function below is implemented in the scope of exercise 3b *
     *****************************************************************/

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of MyCompressibleComponent
     *        at a given temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    {
        return 3900; // [Pa] (at 20C)
    }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return true; }
};

} // end namespace

#endif
#endif
