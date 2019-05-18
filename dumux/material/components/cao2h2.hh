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
 * \brief Material properties of pure Calciumhydroxide \f$CaO2H2\f$.
 */
#ifndef DUMUX_CAO2H2_HH
#define DUMUX_CAO2H2_HH

#include <dumux/common/exceptions.hh>

#include <cmath>
#include <iostream>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the CaO2H2 properties
 */
template <class Scalar>
class CaO2H2
: public Components::Base<Scalar, CaO2H2<Scalar> >
, public Components::Solid<Scalar, CaO2H2<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the CaO2H2.
     */
    static const char *name()
    {
        return "CaO2H2";
    }

    /*!
     * \brief The molar mass of CaOH2 in \f$\mathrm{[kg/mol]}\f$.
     */
    static constexpr Scalar molarMass()
    {
        return 74.093e-3 ;
    }

    /*!
     * \brief The mass density \f$\mathrm{[kg/m^3]}\f$ of CaO2H2.
     */
    static Scalar solidDensity(Scalar temperature)
    {
        return 2200.0; //at 293 K ; Shao et al. (2013)
    }

    /*!
     * \brief The molar density \f$\mathrm{[mol/m^3]}\f$ of CaO2H2.
     * Molar density at 293 K. Literature value from Shao et al. (2013).
     */
    static Scalar solidMolarDensity(Scalar temperature)
    {
        return solidDensity(temperature)/molarMass();
    }

    /*!
     * \brief The specific heat capacity \f$\mathrm{[J/kgK]}\f$ of CaO2H2.
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    {
        return 1530;  //Nagel et al. (2014) : 1530 J/kgK
    }

    /*!
     * \brief The thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    {
        return  0.4;  //Nagel et al. (2014)
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
