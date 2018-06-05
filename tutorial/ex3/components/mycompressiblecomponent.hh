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
 * \brief A fictitious component to be implemented in exercise 3.
 */
#ifndef DUMUX_MYCOMPRESSIBLECOMPONENT_HH
#define DUMUX_MYCOMPRESSIBLECOMPONENT_HH

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
class MyCompressibleComponent
: public Components::Base<Scalar, MyCompressibleComponent<Scalar> >
, public Components::Liquid<Scalar, MyCompressibleComponent<Scalar> >
{

public:
    /*!
     * \brief A human readable name for MyCompressibleComponent.
     */
    static std::string name()
    { return "MyCompressibleComponent"; }

    /*!
     * TODO: Copy the methods implemented in MyIncompressibleComponent and substitute
     *       the density calculation by the expression given in the exercise description.
     */

     /*!
      * \brief Returns true if the liquid phase is assumed to be compressible
      */
     static constexpr bool liquidIsCompressible()
     { return true; }

     /*!
      * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
      */
     static Scalar molarMass()
     {
        // TODO: replace the line below by a meaningful return statement
        DUNE_THROW(Dune::NotImplemented, "Todo: implement molar mass");
     }

     /*!
      * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the liquid component at a given pressure in
      *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
      *
      * \param temperature temperature of component in \f$\mathrm{[K]}\f$
      * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
      */
     static Scalar liquidDensity(Scalar temperature, Scalar pressure)
     {
         // TODO: replace the line below by a meaningful return statement
         DUNE_THROW(Dune::NotImplemented, "Todo: implement liquid density");
     }

     /*!
      * \brief The dynamic liquid viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component.
      *
      * \param temperature temperature of component in \f$\mathrm{[K]}\f$
      * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
      */
     static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
     {
         // TODO: replace the line below by a meaningful return statement
         DUNE_THROW(Dune::NotImplemented, "Todo: implement liquid viscosity");
     }

     /*!
      * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
      *        temperature in \f$\mathrm{[K]}\f$.
      *
      * \param T temperature of the component in \f$\mathrm{[K]}\f$
      */
     static Scalar vaporPressure(Scalar t)
     {
         // TODO: replace the line below by a meaningful return statement
         DUNE_THROW(Dune::NotImplemented, "Todo: implement vapour pressure");
     }
};

} // end namespace

#endif
