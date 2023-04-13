// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Defines the flux traits.
 */
#ifndef DUMUX_FLUX_TRAITS_HH
#define DUMUX_FLUX_TRAITS_HH

#include <type_traits>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Trait of an advection type stating whether it implements a stationary velocity field
 */
template<class AdvectionType>
struct HasStationaryVelocityField : public std::false_type {};

/*!
 * \ingroup Flux
 * \brief Traits of a flux variables type
 */
template<class FluxVariables>
struct FluxTraits
{
    static constexpr bool hasStationaryVelocityField()
    {
        return HasStationaryVelocityField<typename FluxVariables::AdvectionType>::value;
    }
};

} // namespace Dumux

#endif
