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
 * \brief Component traits, i.e. information extracted from components
 */
#ifndef DUMUX_COMPONENT_TRAITS_HH
#define DUMUX_COMPONENT_TRAITS_HH

#include <type_traits>

#include <dumux/material/components/solid.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {

/*!
 * \ingroup Components
 * \brief Component traits, i.e. information extracted from components
 */
template<class Component>
struct ComponentTraits
{
    using Scalar = typename Component::Scalar;

    //! if the component implements a solid state
    static constexpr bool hasSolidState = std::is_base_of<Components::Solid<Scalar, Component>, Component>::value;

    //! if the component implements a liquid state
    static constexpr bool hasLiquidState = std::is_base_of<Components::Liquid<Scalar, Component>, Component>::value;

    //! if the component implements a gaseous state
    static constexpr bool hasGasState = std::is_base_of<Components::Gas<Scalar, Component>, Component>::value;

    //! if the component implements an ion
    static constexpr bool isIon = std::is_base_of<Components::Ion<Scalar, Component>, Component>::value;
};

} // end namespace Dumux

#endif
