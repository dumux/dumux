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
 * \ingroup Discretization
 * \brief The available discretization methods in Dumux
 */
#ifndef DUMUX_DISCRETIZARION_METHOD_HH
#define DUMUX_DISCRETIZARION_METHOD_HH

#include <ostream>
#include <string>

namespace Dumux {

/*!
 * \brief The available discretization methods in Dumux
 * \ingroup Discretization
 * \note Use none if specifying a discretization method is required but
 *       the class in question is not specific to a a discretization method
 *       or the classification is non-applicable
 */
enum class DiscretizationMethod
{
    none, box, cctpfa, ccmpfa, staggered, fem
};

/*!
 * \brief Convert discretization method to string
 * \ingroup Discretization
 */
inline std::string toString(DiscretizationMethod m)
{
    switch (m)
    {
        case DiscretizationMethod::box: return "box";
        case DiscretizationMethod::cctpfa: return "cctpfa";
        case DiscretizationMethod::ccmpfa: return "ccmpfa";
        case DiscretizationMethod::staggered: return "staggered";
        case DiscretizationMethod::fem: return "fem";
        default: return "none";
    }
}

/*!
 * \brief Write discretization method to stream
 * \ingroup Discretization
 */
inline std::ostream& operator<<(std::ostream& stream, DiscretizationMethod m)
{ stream << toString(m); return stream; }

} // end namespace Dumux

#endif
