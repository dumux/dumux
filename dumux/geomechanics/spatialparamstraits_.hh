// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsModels
 * \brief helper struct detecting if the user-defined spatial params class has a lameParamsAtPos function
 */
#ifndef DUMUX_GEOMECHANICS_SPATIAL_PARAMS_TRAITS__HH
#define DUMUX_GEOMECHANICS_SPATIAL_PARAMS_TRAITS__HH

#include <utility>

#ifndef DOXYGEN
namespace Dumux::Detail {

// helper struct detecting if the user-defined spatial params class has a lameParamsAtPos function
template<class GlobalPosition>
struct hasLameParamsAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.lameParamsAtPos(std::declval<GlobalPosition>()))
    {}
};

} // end namespace Dumux::Detail
#endif  // DOXYGEN
#endif
