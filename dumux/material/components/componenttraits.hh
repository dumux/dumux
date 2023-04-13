// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

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
