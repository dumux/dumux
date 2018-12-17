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
 * \brief A class for the Ca2+ (Calcium ion) component properties
 */
#ifndef DUMUX_CA_ION_HH
#define DUMUX_CA_ION_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Ca2+ (Calcium ion) component properties.
 */
template <class Scalar>
class CalciumIon
: public Components::Base<Scalar, CalciumIon<Scalar> >
, public Components::Ion<Scalar, CalciumIon<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the Ca ion.
     */
    static std::string name()
    { return "Ca2+"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the Ca ion.
     */
    static constexpr Scalar molarMass()
    { return 40.078e-3; } // kg/mol

    /*!
     * \brief The charge balance of the Ca ion.
     */
    static constexpr int charge()
    { return +2; }

};

} // end namespace Components
} // end namespace Dumux

#endif
