// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Common
 * \ingroup Tests
 * \brief A property system example using property type aliases
 */
#include <iostream>
#include <type_traits>
#include <dumux/common/properties/propertysystem.hh>

namespace Dumux::Properties {
DUMUX_DEFINE_PROPERTY(GasUsage)
DUMUX_DEFINE_PROPERTY(TopSpeed)
DUMUX_DEFINE_PROPERTY(NumSeats)
DUMUX_DEFINE_PROPERTY(AutomaticTransmission)
DUMUX_DEFINE_PROPERTY(CannonCaliber)
DUMUX_DEFINE_PROPERTY(Payload)
} // end namespace Dumux::Properties

namespace Dumux::Properties::TTag {

struct CompactCar {
    template<class TypeTag>
    using TopSpeed = std::integral_constant<int, 30*getPropValue<TypeTag, Properties::GasUsage>()>;
    using NumSeats = std::integral_constant<int, 5>;
    using GasUsage = std::integral_constant<int, 4>;
};

struct Truck {
    using TopSpeed = std::integral_constant<int, 100>;
    using NumSeats = std::integral_constant<int, 2>;
    using GasUsage = std::integral_constant<int, 18>;
    using Payload = std::integral_constant<int, 35>;
};

struct Tank {
    using TopSpeed = std::integral_constant<int, 60>;
    using GasUsage = std::integral_constant<int, 18>;
    using CannonCaliber = std::integral_constant<int, 120>;
};

struct Sedan {
    using InheritsFrom = std::tuple<CompactCar>;
    using GasUsage = std::integral_constant<int, 7>;
    using AutomaticTransmission = std::true_type;
};

struct Pickup {
    using InheritsFrom = std::tuple<Truck, Sedan>;
    using TopSpeed = std::integral_constant<int, 120>;
    using Payload = std::integral_constant<int, 5>;
};

struct HummerH1 {
    using InheritsFrom = std::tuple<Tank, Pickup>;
    using TopSpeed = std::integral_constant<int, getPropValue<Pickup, Properties::TopSpeed>()>;
};

} // end namespace Dumux::Properties::TTag

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace Properties;

    // test assumptions on the compact car
    static_assert(getPropValue<TTag::CompactCar, Properties::TopSpeed>() == 30*getPropValue<TTag::CompactCar, Properties::GasUsage>());
    static_assert(getPropValue<TTag::CompactCar, Properties::NumSeats>() == 5);
    static_assert(getPropValue<TTag::CompactCar, Properties::GasUsage>() == 4);

    // test assumptions on the truck
    static_assert(getPropValue<TTag::Truck, Properties::TopSpeed>() == 100);
    static_assert(getPropValue<TTag::Truck, Properties::NumSeats>() == 2);
    static_assert(getPropValue<TTag::Truck, Properties::Payload>() == 35);

    // test assumptions on tank
    static_assert(getPropValue<TTag::Tank, Properties::TopSpeed>() == 60);
    static_assert(getPropValue<TTag::Tank, Properties::GasUsage>() == 18);
    static_assert(getPropValue<TTag::Tank, Properties::CannonCaliber>() == 120);

    // test assumptions on sedan
    static_assert(getPropValue<TTag::Sedan, Properties::AutomaticTransmission>());
    static_assert(getPropValue<TTag::Sedan, Properties::GasUsage>() == 7);
    static_assert(getPropValue<TTag::Sedan, Properties::TopSpeed>() == 30*getPropValue<TTag::Sedan, Properties::GasUsage>());
    static_assert(getPropValue<TTag::Sedan, Properties::NumSeats>() == getPropValue<TTag::CompactCar, Properties::NumSeats>());

    // test assumption on pickup
    static_assert(getPropValue<TTag::Pickup, Properties::TopSpeed>() == 120);
    static_assert(getPropValue<TTag::Pickup, Properties::Payload>() == 5);
    static_assert(getPropValue<TTag::Pickup, Properties::AutomaticTransmission>() == getPropValue<TTag::Sedan, Properties::AutomaticTransmission>());
    static_assert(getPropValue<TTag::Pickup, Properties::NumSeats>() == getPropValue<TTag::Truck, Properties::NumSeats>());

    // test assumption on hummer
    static_assert(getPropValue<TTag::HummerH1, Properties::TopSpeed>() == getPropValue<TTag::Pickup, Properties::TopSpeed>());
    static_assert(getPropValue<TTag::HummerH1, Properties::Payload>() == getPropValue<TTag::Pickup, Properties::Payload>());
    static_assert(getPropValue<TTag::HummerH1, Properties::AutomaticTransmission>() == getPropValue<TTag::Pickup, Properties::AutomaticTransmission>());
    static_assert(getPropValue<TTag::HummerH1, Properties::NumSeats>() == getPropValue<TTag::Pickup, Properties::NumSeats>());
    static_assert(getPropValue<TTag::HummerH1, Properties::GasUsage>() == getPropValue<TTag::Tank, Properties::GasUsage>());
    static_assert(getPropValue<TTag::HummerH1, Properties::CannonCaliber>() == getPropValue<TTag::Tank, Properties::CannonCaliber>());

    std::cout << "-- Top speed of sedan: " << getPropValue<Properties::TTag::Sedan, Properties::TopSpeed>() << "\n";
    std::cout << "-- Top speed of truck: " << getPropValue<Properties::TTag::Truck, Properties::TopSpeed>() << "\n";

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
