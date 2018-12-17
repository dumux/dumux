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
 * \brief A class for the CO3 ion properties.
 */
#ifndef DUMUX_CO3_ION_HH
#define DUMUX_CO3_ION_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the CO3 fluid properties.
 */
template <class Scalar>
class CarbonateIon
: public Components::Base<Scalar, CarbonateIon<Scalar> >
, public Components::Ion<Scalar, CarbonateIon<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the CO3 ion.
     */
    static std::string name()
    { return "CO3-"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the CO3 ion.
     */
    static constexpr Scalar molarMass()
    { return 60.0092e-3; } // kg/mol

    /*!
     * \brief The charge balance of the CO3 ion.
     */
    static constexpr int charge()
    { return -2; }

};

} // end namespace Components
} // end namespace Dumux

#endif
