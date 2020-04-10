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
 * \brief A class for the Ammonia (NH3) component properties
 */
#ifndef DUMUX_MATERIAL_COMPONENTS_NH3_HH
#define DUMUX_MATERIAL_COMPONENTS_NH3_HH

#include <dumux/material/components/base.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Ammonia (NH3) component properties
 */
template <class Scalar>
class Ammonia
: public Components::Base<Scalar, Ammonia<Scalar> >
{
public:

    /*!
     * \brief A human readable name for NH3.
     */
    static std::string name()
    { return "NH3"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of NH3.
     */
    static Scalar molarMass()
    { return 0.017031; } // kg/mol

};

} // end namespace Components
} // end namespace Dumux

#endif
