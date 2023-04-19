// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Python wrapper for fluid/solid system components
 */

#ifndef DUMUX_PYTHON_MATERIAL_COMPONENT_HH
#define DUMUX_PYTHON_MATERIAL_COMPONENT_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/material/components/componenttraits.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/typetraits/typetraits.hh>


namespace Dumux::Python::Detail {
//! Helper struct to deactivate static assertions in component's base classes.
struct DisableStaticAssert {};

} // end namespace Dumux::Python::Detail

namespace Dumux {
/*!
 * \brief Specialization of Dumux::AlwaysFalse for the struct defined
 *        above. This is done in order to deactivate the static_assert in
 *        the base classes of components. If the base class function is compiled
 *        we do not call it (see below).
 */
template<>
struct AlwaysFalse<Dumux::Python::Detail::DisableStaticAssert> : public std::true_type {};

} // end namespace Dumux

namespace Dumux::Python::Detail {

struct Name { template<class C> auto operator()(C&& c) -> decltype(C::template name<DisableStaticAssert>()) {} };
struct MolarMass { template<class C> auto operator()(C&& c) -> decltype(C::template molarMass<DisableStaticAssert>()) {} };
struct VaporPressure { template<class C> auto operator()(C&& c) -> decltype(C::template vaporPressure<DisableStaticAssert>(0.0)) {} };

struct LiquidIsCompressible { template<class C> auto operator()(C&& c) -> decltype(C::template liquidIsCompressible<DisableStaticAssert>()) {} };
struct LiquidDensity { template<class C> auto operator()(C&& c) -> decltype(C::template liquidDensity<DisableStaticAssert>(0.0, 0.0)) {} };
struct LiquidMolarDensity { template<class C> auto operator()(C&& c) -> decltype(C::template liquidMolarDensity<DisableStaticAssert>(0.0, 0.0)) {} };
struct LiquidViscosity { template<class C> auto operator()(C&& c) -> decltype(C::template liquidViscosity<DisableStaticAssert>(0.0, 0.0)) {} };
struct LiquidEnthalpy { template<class C> auto operator()(C&& c) -> decltype(C::template liquidEnthalpy<DisableStaticAssert>(0.0, 0.0)) {} };
struct LiquidInternalEnergy { template<class C> auto operator()(C&& c) -> decltype(C::template liquidInternalEnergy<DisableStaticAssert>(0.0, 0.0)) {} };
struct LiquidHeatCapacity { template<class C> auto operator()(C&& c) -> decltype(C::template liquidHeatCapacity<DisableStaticAssert>(0.0, 0.0)) {} };
struct LiquidThermalConductivity { template<class C> auto operator()(C&& c) -> decltype(C::template liquidThermalConductivity<DisableStaticAssert>(0.0, 0.0)) {} };

struct GasIsIdeal { template<class C> auto operator()(C&& c) -> decltype(C::template gasIsIdeal<DisableStaticAssert>()) {} };
struct GasIsCompressible { template<class C> auto operator()(C&& c) -> decltype(C::template gasIsCompressible<DisableStaticAssert>()) {} };
struct GasDensity { template<class C> auto operator()(C&& c) -> decltype(C::template gasDensity<DisableStaticAssert>(0.0, 0.0)) {} };
struct GasMolarDensity { template<class C> auto operator()(C&& c) -> decltype(C::template gasMolarDensity<DisableStaticAssert>(0.0, 0.0)) {} };
struct GasViscosity { template<class C> auto operator()(C&& c) -> decltype(C::template gasViscosity<DisableStaticAssert>(0.0, 0.0)) {} };
struct GasEnthalpy { template<class C> auto operator()(C&& c) -> decltype(C::template gasEnthalpy<DisableStaticAssert>(0.0, 0.0)) {} };
struct GasInternalEnergy { template<class C> auto operator()(C&& c) -> decltype(C::template gasInternalEnergy<DisableStaticAssert>(0.0, 0.0)) {} };
struct GasHeatCapacity { template<class C> auto operator()(C&& c) -> decltype(C::template gasHeatCapacity<DisableStaticAssert>(0.0, 0.0)) {} };
struct GasThermalConductivity { template<class C> auto operator()(C&& c) -> decltype(C::template gasThermalConductivity<DisableStaticAssert>(0.0, 0.0)) {} };

struct SolidIsCompressible { template<class C> auto operator()(C&& c) -> decltype(C::template solidIsCompressible<DisableStaticAssert>()) {} };
struct SolidDensity { template<class C> auto operator()(C&& c) -> decltype(C::template solidDensity<DisableStaticAssert>(0.0)) {} };
struct SolidThermalConductivity { template<class C> auto operator()(C&& c) -> decltype(C::template solidThermalConductivity<DisableStaticAssert>(0.0)) {} };
struct SolidHeatCapacity { template<class C> auto operator()(C&& c) -> decltype(C::template solidHeatCapacity<DisableStaticAssert>(0.0)) {} };

struct Charge { template<class C> auto operator()(C&& c) -> decltype(C::template charge<DisableStaticAssert>()) {} };
struct CriticalTemperature { template<class C> auto operator()(C&& c) -> decltype(C::template criticalTemperature<DisableStaticAssert>()) {} };
struct CriticalPressure { template<class C> auto operator()(C&& c) -> decltype(C::template criticalPressure<DisableStaticAssert>()) {} };

} // end namespace Dumux::Python::Detail


namespace Dumux::Python {

template <class Comp, class... options>
void registerComponent(pybind11::handle scope,
                       pybind11::class_<Comp, options...> cls)
{
    using pybind11::operator""_a;

    cls.def(pybind11::init());

    // helper lambda to check if the component implements a certain function
    constexpr auto implements = [](const auto& c)
    {
        return !decltype(isValid(c)(std::declval<Comp>()))::value;
    };

    if constexpr (implements(Detail::Name{}))
        cls.def_property_readonly_static("name", &Comp::name);
    if constexpr (implements(Detail::MolarMass{}))
        cls.def_property_readonly_static("molarMass", &Comp::molarMass);

    if constexpr (ComponentTraits<Comp>::hasLiquidState)
    {
        if constexpr (implements(Detail::LiquidIsCompressible{}))
            cls.def_property_readonly_static("liquidIsCompressible", &Comp::liquidIsCompressible);
        if constexpr (implements(Detail::LiquidDensity{}))
            cls.def_static("liquidDensity", &Comp::liquidDensity, "temperature"_a, "pressure"_a);
        if constexpr (implements(Detail::LiquidMolarDensity{}))
            cls.def_static("liquidMolarDensity", &Comp::liquidMolarDensity, "temperature"_a, "pressure"_a);
        if constexpr (implements(Detail::LiquidViscosity{}))
            cls.def_static("liquidViscosity", &Comp::liquidViscosity, "temperature"_a, "pressure"_a);
        if constexpr (implements(Detail::LiquidEnthalpy{}))
            cls.def_static("liquidEnthalpy", &Comp::liquidEnthalpy, "temperature"_a, "pressure"_a);
        if constexpr (implements(Detail::LiquidInternalEnergy{}))
            cls.def_static("liquidInternalEnergy", &Comp::liquidInternalEnergy, "temperature"_a, "pressure"_a);
        if constexpr (implements(Detail::LiquidHeatCapacity{}))
            cls.def_static("liquidHeatCapacity", &Comp::liquidHeatCapacity, "temperature"_a, "pressure"_a);
        if constexpr (implements(Detail::LiquidThermalConductivity{}))
            cls.def_static("liquidThermalConductivity", &Comp::liquidThermalConductivity, "temperature"_a, "pressure"_a);
        if constexpr(implements(Detail::VaporPressure{}))
            cls.def_static("vaporPressure", &Comp::vaporPressure, "temperature"_a);
    }

    if constexpr (ComponentTraits<Comp>::hasGasState)
    {
        if constexpr(implements(Detail::GasDensity{}))
            cls.def_static("gasDensity", &Comp::gasDensity, "temperature"_a, "pressure"_a);
        if constexpr(implements(Detail::GasMolarDensity{}))
            cls.def_static("gasMolarDensity", &Comp::gasMolarDensity, "temperature"_a, "pressure"_a);
        if constexpr(implements(Detail::GasIsCompressible{}))
            cls.def_property_readonly_static("gasIsCompressible", &Comp::gasIsCompressible);
        if constexpr(implements(Detail::GasIsIdeal{}))
            cls.def_property_readonly_static("gasIsIdeal", &Comp::gasIsIdeal);
        if constexpr(implements(Detail::GasViscosity{}))
            cls.def_static("gasViscosity", &Comp::gasViscosity, "temperature"_a, "pressure"_a);
        if constexpr(implements(Detail::GasEnthalpy{}))
            cls.def_static("gasEnthalpy", &Comp::gasEnthalpy, "temperature"_a, "pressure"_a);
        if constexpr(implements(Detail::GasInternalEnergy{}))
            cls.def_static("gasInternalEnergy", &Comp::gasInternalEnergy, "temperature"_a, "pressure"_a);
        if constexpr(implements(Detail::GasHeatCapacity{}))
            cls.def_static("gasHeatCapacity", &Comp::gasHeatCapacity, "temperature"_a, "pressure"_a);
        if constexpr(implements(Detail::GasThermalConductivity{}))
            cls.def_static("gasThermalConductivity", &Comp::gasThermalConductivity, "temperature"_a, "pressure"_a);
    }

    if constexpr (ComponentTraits<Comp>::hasSolidState)
    {
        if constexpr(implements(Detail::SolidIsCompressible{}))
            cls.def_property_readonly_static("solidIsCompressible", &Comp::solidIsCompressible);
        if constexpr(implements(Detail::SolidDensity{}))
            cls.def_property_readonly_static("solidDensity", &Comp::solidDensity, "temperature_a");
        if constexpr(implements(Detail::SolidThermalConductivity{}))
            cls.def_static("solidThermalConductivity", &Comp::solidThermalConductivity, "temperature_a");
        if constexpr(implements(Detail::SolidHeatCapacity{}))
            cls.def_static("solidHeatCapacity", &Comp::solidHeatCapacity, "temperature_a");
    }

    if constexpr (ComponentTraits<Comp>::isIon)
    {
        if constexpr(implements(Detail::Charge{}))
            cls.def_property_readonly_static("charge", &Comp::charge);
    }

    if constexpr (ComponentTraits<Comp>::hasLiquidState || ComponentTraits<Comp>::hasGasState)
    {
        if constexpr(implements(Detail::CriticalTemperature{}))
            cls.def_property_readonly_static("criticalTemperature", &Comp::criticalTemperature);
        if constexpr(implements(Detail::CriticalPressure{}))
            cls.def_property_readonly_static("criticalPressure", &Comp::criticalPressure);
    }
}


} // namespace Dumux::Python

#endif
