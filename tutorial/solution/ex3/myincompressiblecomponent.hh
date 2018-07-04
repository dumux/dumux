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
 * \brief A fictitious component to be implemented in tutorial exercise 3.
 */
#ifndef DUMUX_MYINCOMPRESSIBLECOMPONENT_HH
#define DUMUX_MYINCOMPRESSIBLECOMPONENT_HH

#include <dumux/material/idealgas.hh>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>

namespace Dumux
{
/*!
 * \ingroup Components
 * \brief A fictitious component to be implemented in exercise 3.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class MyIncompressibleComponent
: public Components::Base<Scalar, MyIncompressibleComponent<Scalar> >
, public Components::Liquid<Scalar, MyIncompressibleComponent<Scalar> >
{
public:
    /*!
     * \brief A human readable name for MyIncompressibleComponent.
     */
    static std::string name()
    { return "MyIncompressibleComponent"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of MyIncompressibleComponent.
     */
    static Scalar molarMass()
    {
        return 131.39e-3; // [kg/mol]
    }

    /*!
     * \brief The density of MyIncompressibleComponent at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1460.0; // [kg/m^3]
    }

    /*!
     * \brief The molar density of MyIncompressibleComponent in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of MyIncompressibleComponent.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 5.7e-4;// [Pa*s]
    }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }
};

} // end namespace

#endif
