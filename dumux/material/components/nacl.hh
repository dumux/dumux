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
 * \brief Material properties of pure salt \f$NaCl\f$.
 */
#ifndef DUMUX_NACL_HH
#define DUMUX_NACL_HH

#include <dumux/common/exceptions.hh>

#include <cmath>
#include <iostream>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the NaCl properties
 */
template <class Scalar>
class NaCl
: public Components::Base<Scalar, NaCl<Scalar> >
, public Components::Solid<Scalar, NaCl<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the NaCl.
     */
    static std::string name()
    {
        return "NaCl";
    }

    /*!
     * \brief The molar mass of NaCl in \f$\mathrm{[kg/mol]}\f$.
     */
    static constexpr Scalar molarMass()
    {
        return 58.4428e-3 ;
    }

    /*!
     * \brief The mass density \f$\mathrm{[kg/m^3]}\f$ of NaCl.
     */
    static Scalar solidDensity(Scalar temperature)
    {
        return 2165.0;
    }

    /*!
     * \brief The mass density \f$\mathrm{[kg/m^3]}\f$ of NaCl.
     */
    static Scalar solidMolarDensity(Scalar temperature)
    {
        return solidDensity(temperature)/molarMass();
    }

    /*!
     * \brief The specific heat capacity \f$\mathrm{[J/molK]}\f$ of NaCl.
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    {
        return 50.50;
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    {
        return 6.49;
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
