// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsModels
 * \brief helper struct detecting if the user-defined spatial params class has a lameParamsAtPos function
 */
#ifndef DUMUX_POROMECHANICS_SPATIAL_PARAMS_TRAITS__HH
#define DUMUX_POROMECHANICS_SPATIAL_PARAMS_TRAITS__HH

#include <utility>
#include <dumux/solidmechanics/elastic/spatialparamstraits_.hh>

#ifndef DOXYGEN
namespace Dumux::Detail {

// helper struct detecting if the user-defined problem class has an effectiveFluidDensityAtPos function
template<class GlobalPosition>
struct hasEffFluidDensityAtPos
{
    template<class Problem>
    auto operator()(const Problem& a)
    -> decltype(a.effectiveFluidDensityAtPos(std::declval<GlobalPosition>()))
    {}
};

// helper struct detecting if the user-defined problem class has an effectivePorePressureAtPos function
template<class GlobalPosition>
struct hasEffPorePressureAtPos
{
    template<class Problem>
    auto operator()(const Problem& a)
    -> decltype(a.effectivePorePressureAtPos(std::declval<GlobalPosition>()))
    {}
};

// helper struct detecting if the user-defined spatial params class has a reactiveVolumeFractionAtPos function
template<class GlobalPosition, class SolidSystem>
struct hasReactiveVolumeFractionAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.template reactiveVolumeFractionAtPos<SolidSystem>(std::declval<GlobalPosition>(), 0))
    {}
};

// helper struct detecting if the user-defined spatial params class has a biotCoefficientAtPos function
template<class GlobalPosition>
struct hasBiotCoeffAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.biotCoefficientAtPos(std::declval<GlobalPosition>()))
    {}
};

} // end namespace Detail
#endif  // DOXYGEN
#endif
