// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Python wrapper for volume variables (finite volume schemes)
 */

#ifndef DUMUX_PYTHON_COMMON_VOLUMEVARIABLES_HH
#define DUMUX_PYTHON_COMMON_VOLUMEVARIABLES_HH

#include <dune/common/std/type_traits.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dune/python/common/typeregistry.hh>
#include <dune/common/classname.hh>

namespace Dumux::Python::Impl {

// helper structs and functions for detecting member functions in volVars
template <class VolumeVariables>
using PhaseTemperatureDetector = decltype(std::declval<VolumeVariables>().temperature(0));

template<class VolumeVariables>
static constexpr bool hasPhaseTemperature()
{ return Dune::Std::is_detected<PhaseTemperatureDetector, VolumeVariables>::value; }

template <class VolumeVariables>
using MoleFractionDetector = decltype(std::declval<VolumeVariables>().moleFraction(0, 0));

template<class VolumeVariables>
static constexpr bool hasMoleFraction()
{ return Dune::Std::is_detected<MoleFractionDetector, VolumeVariables>::value; }

template <class VolumeVariables>
using MassFractionDetector = decltype(std::declval<VolumeVariables>().massFraction(0, 0));

template<class VolumeVariables>
static constexpr bool hasMassFraction()
{ return Dune::Std::is_detected<MassFractionDetector, VolumeVariables>::value; }

template <class VolumeVariables>
using SaturationDetector = decltype(std::declval<VolumeVariables>().saturation(0));

template<class VolumeVariables>
static constexpr bool hasSaturation()
{ return Dune::Std::is_detected<SaturationDetector, VolumeVariables>::value; }

template <class VolumeVariables>
using PermeabilityDetector = decltype(std::declval<VolumeVariables>().permeability());

template<class VolumeVariables>
static constexpr bool hasPermeability()
{ return Dune::Std::is_detected<PermeabilityDetector, VolumeVariables>::value; }

template<class VolumeVariables>
void registerVolumeVariables(pybind11::handle scope)
{
    using namespace Dune::Python;

    auto [cls, addedToRegistry] = insertClass<VolumeVariables>(
        scope, "VolumeVariables",
        GenerateTypeName(Dune::className<VolumeVariables>()),
        IncludeFiles{}
    );

    if (addedToRegistry)
    {
        using pybind11::operator""_a;

        cls.def("pressure", &VolumeVariables::pressure, "phaseIdx"_a=0);
        cls.def("density", &VolumeVariables::density, "phaseIdx"_a=0);
        cls.def("temperature", &VolumeVariables::temperature);

        if constexpr(hasSaturation<VolumeVariables>())
            cls.def("saturation", &VolumeVariables::saturation, "saturation"_a=0);
        if constexpr(hasPhaseTemperature<VolumeVariables>())
            cls.def("temperature", &VolumeVariables::temperature, "phaseIdx"_a=0);
        if constexpr(hasMoleFraction<VolumeVariables>())
            cls.def("moleFraction", &VolumeVariables::moleFraction, "phaseIdx"_a, "compIdx"_a);
        if constexpr(hasMassFraction<VolumeVariables>())
            cls.def("massFraction", &VolumeVariables::massFraction, "phaseIdx"_a, "compIdx"_a);
        if constexpr(hasPermeability<VolumeVariables>())
            cls.def("permeability", &VolumeVariables::permeability);
    }
}

} // namespace Dumux::Python

#endif
