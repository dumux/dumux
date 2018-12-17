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
 * \brief Interface for components that are ions.
 */
#ifndef DUMUX_COMPONENT_ION_HH
#define DUMUX_COMPONENT_ION_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Interface for components that are ions.
 */
template<class Scalar, class Component>
class Ion
{
public:
    /*!
     * \brief Returns the charge of the ion.
     */
    template<class C = Component>
    static constexpr int charge()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: charge()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }
};

} // end namespace Components
} // end namespace Dumux

#endif
