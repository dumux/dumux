@addtogroup Properties

The DuMu<sup>x</sup> property system is based on the concept of type traits
with added inheritance. It is implemented using template metaprogramming.

In the context of the DuMu<sup>x</sup> property system, a property is an arbitrary
class which may contain type definitions, values and methods.
Just like normal classes, properties can be arranged in hierarchies. In
the context of the DuMu<sup>x</sup> property system, nodes of the inheritance
hierarchy are called **type tags**.

It also supports property **nesting**. Property nesting means that the definition of
a property can depend on the value of other properties which may be
defined for arbitrary levels of the inheritance hierarchy.

This section gives a high level overview over the property system's design and principle ideas
illustrated by self-contained examples.

# How to use the property system

All source files which use the property system should include
the header file [dumux/common/properties.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/common/properties.hh).
Declaration of type tags and property tags as well as defining properties must be done inside the
namespace Dumux::Properties.

## Defining type tags

New nodes in the type tag hierarchy can be defined in the Dumux::Properties::TTag namespace using

```cpp
// Create new type tags
namespace TTag {
struct NewTypeTagName { using InheritsFrom = std::tuple<BaseTagName1, BaseTagName2, ...>; };
} // end namespace TTag
```

where the `InheritsFrom` alias is optional.
To avoid inconsistencies in the hierarchy, each type tag may be defined only
once for a program. If you call Dumux::GetProp the property system will first look for the properties defined in
`BaseTagName1` in the `InheritsFrom` list.
If a defined property is found this property is returned.
If no defined property is found the search will continue in the ancestors of `BaseTagName1`.
If again no defined property is found the search will continue in the second `BaseTagName2` in the list, and so on.
If no defined property is found at all, a compiler error is triggered.

Example:
```cpp
namespace Dumux::Properties::TTag {
struct MyBaseTypeTag1 {};
struct MyBaseTypeTag2 {};

struct MyDerivedTypeTag
{
    using InheritsFrom = std::tuple<
        MyBaseTypeTag1, MyBaseTypeTag2
    >;
};
} // end namespace Dumux::Properties::TTag
```

## Defining new properties
New property tags are defined using the macro DUMUX_DEFINE_PROPERTY, e.g.

```cpp
namespace Dumux::Properties {
DUMUX_DEFINE_PROPERTY(MyPropertyTag)
} // end namespace Dumux::Properties
```

Essentially this corresponds to the following code

```cpp
namespace Dumux::Properties {
template<class TypeTag, class MyTypeTag>
struct MyPropertyTag { using type = UndefinedProperty; };
} // end namespace Dumux::Properties
```

If you need to forward declare a property you can use

```cpp
// forward declaration
namespace Dumux::Properties {
template<class TypeTag, class MyTypeTag>
struct NewPropTagName;
} // end namespace Dumux::Properties
```

### Specializing properties

The value of a property on a given node of the type tag hierarchy is
defined by means of [partial template specialization](https://en.cppreference.com/w/cpp/language/partial_specialization)

```cpp
template<class TypeTag>
struct PropertyTagName<TypeTag, TTag::TypeTagName>
{
  // arbitrary body of a struct
};
```

where here the property `PropertyTagName` is specialized
for the tag `TTag::TypeTagName`.
The body typically contains either the [type alias](https://en.cppreference.com/w/cpp/language/type_alias)
`type`, or a `static constexpr` data member `value`.
However, you can of course write in the body whatever you like.

```cpp
template<class TypeTag>
struct PropertyTagName<TypeTag, TTag::TypeTagName> { using type = <type>; };

template<class TypeTag>
struct PropertyTagName<TypeTag, TTag::TypeTagName> { static constexpr bool value = <booleanValue>; };

template<class TypeTag>
struct PropertyTagName<TypeTag, TTag::TypeTagName> { static constexpr int value = <integerValue>; };
```

Here is an example including a type tag, property definitions and specializations:

```cpp
namespace Dumux::Properties {

// Create new type tag
namespace TTag {
struct MyTypeTag {};
} // end namespace TTag

// Define some properties
DUMUX_DEFINE_PROPERTY(MyCustomProperty)
DUMUX_DEFINE_PROPERTY(MyType)
DUMUX_DEFINE_PROPERTY(MyBoolValue)
DUMUX_DEFINE_PROPERTY(MyIntValue)
DUMUX_DEFINE_PROPERTY(MyScalarValue)

// Set the properties for the new type tag
template<class TypeTag>
struct MyCustomProperty<TypeTag, TTag::MyTypeTag>
{
    static void print()
    { std::cout << "Hello, World!\n"; }
};

template<class TypeTag>
struct MyType<TypeTag, TTag::MyTypeTag> { using type = unsigned int; };

template<class TypeTag>
struct MyBoolValue<TypeTag, TTag::MyTypeTag> { static constexpr bool value = true; };

template<class TypeTag>
struct MyIntValue<TypeTag, TTag::MyTypeTag> { static constexpr int value = 12345; };

template<class TypeTag>
struct MyScalarValue<TypeTag, TTag::MyTypeTag> { static constexpr double value = 12345.67890; };

} // end namespace Dumux::Properties
```

## Retrieving property values

The type of a property can be retrieved using

```cpp
using Prop = GetProp<TypeTag, Properties::PropertyTag>;
```

There is a helper struct and a helper function to retrieve the `type` and `value` members of a property

```cpp
using PropType = GetPropType<TypeTag, Properties::PropertyTag>;
constexpr auto propValue = getPropValue<TypeTag, Properties::PropertyTag>();
```

Example:

```cpp
template <TypeTag>
class MyClass {
    // retrieve the ::value attribute of the 'UseMoles' property
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr bool useMoles2 = GetProp<TypeTag, Properties::UseMoles>::value; // equivalent

    // retrieve the ::type attribute of the 'Scalar' property
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Scalar2 = GetProp<TypeTag, Properties::Scalar>::type; // equivalent
};
```

## Nesting property definitions

Inside property definitions there is access to all other properties
which are defined somewhere on the type tag hierarchy. The node for
which the current property is requested is available via the template argument
`TypeTag`. Inside property class bodies `GetPropType` can be used to
retrieve other properties and create aliases.

Example:

```cpp
template<class TypeTag>
struct Vector<TypeTag, TTag::MyModelTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = std::vector<Scalar>;
};
```

# A self-contained property example

As a concrete example, let us consider some kinds of cars: Compact
cars, sedans, trucks, pickups, military tanks and the Hummer-H1 sports
utility vehicle. Since all these cars share some characteristics, it
makes sense to inherit those from the closest matching car type and
only specify the properties which are different. Thus, an inheritance
diagram for the car types above might look like outlined in the figure.

@mermaid{prop_example}

## Defining type tags, properties, and specializations

Using the DuMu<sup>x</sup> property system,
this type hierarchy with inheritance is
defined by:

```cpp
#include <dumux/common/propertysystem.hh>
#include <iostream>

namespace Dumux::Properties::TTag {
struct CompactCar {};
struct Truck {};
struct Tank {};

struct Sedan { using InheritsFrom = std::tuple<CompactCar>; };
struct Pickup { using InheritsFrom = std::tuple<Truck, Sedan>; };
struct HummerH1 { using InheritsFrom = std::tuple<Tank, Pickup>; };
} // end namespace Dumux::Properties::TTag
```

The Figure lists a few property names which
make sense for at least one of the nodes.
These property names can be defined as
follows:

```cpp
template<class TypeTag, class MyTypeTag>
struct GasUsage { using type = UndefinedProperty; }; // [l/100km]
template<class TypeTag, class MyTypeTag>
struct TopSpeed { using type = UndefinedProperty; }; // [km/h]
template<class TypeTag, class MyTypeTag>
struct NumSeats { using type = UndefinedProperty; }; // []
template<class TypeTag, class MyTypeTag>
struct AutomaticTransmission { using type = UndefinedProperty; }; // true/false
template<class TypeTag, class MyTypeTag>
struct CannonCaliber { using type = UndefinedProperty; }; // [mm]
template<class TypeTag, class MyTypeTag>
struct Payload { using type = UndefinedProperty; }; // [t]
```

So far, the inheritance hierarchy and the property names are completely
separate. What is missing is setting some values for the property
names on specific nodes of the inheritance hierarchy. Let us assume
the following:

* For a compact car, the top speed is the gas usage per 100 km
  times 30, the number of seats is 5 and the gas usage is 4 l/100km.
* A truck is by law limited to 100 km/h top speed, the number
  of seats is 2, it uses 18 l/100km and has a cargo payload of 35 tons.
* A tank exhibits a top speed of 60 km/h, uses 65 l/100km
  and features a 120 mm diameter cannon
* A sedan has a gas usage of 7 l/100km, as well as an automatic
  transmission. In every other aspect it is like a compact car.
* A pick-up truck has a top speed of 120 km/h and a payload of
  5 tons. In every other aspect it is like a sedan or a truck but if in
  doubt, it is more like a truck.
* The Hummer-H1 SUV exhibits the same top speed as a pick-up
  truck. In all other aspects it is similar to a pickup and a tank,
  but, if in doubt, more like a tank.

Using the DuMu<sup>x</sup> property system, these assumptions are formulated
using

```cpp
template<class TypeTag>
struct TopSpeed<TypeTag, TTag::CompactCar>
{
    static constexpr int value
        = getPropValue<TypeTag, Properties::GasUsage>() * 30;
};

template<class TypeTag>
struct NumSeats<TypeTag, TTag::CompactCar>
{ static constexpr int value = 5; };

template<class TypeTag>
struct GasUsage<TypeTag, TTag::CompactCar>
{ static constexpr int value = 4; };

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::Truck>
{ static constexpr int value = 100; };

template<class TypeTag>
struct NumSeats<TypeTag, TTag::Truck>
{ static constexpr int value = 2; };

template<class TypeTag>
struct GasUsage<TypeTag, TTag::Truck>
{ static constexpr int value = 18; };

template<class TypeTag>
struct Payload<TypeTag, TTag::Truck>
{ static constexpr int value = 35; };

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::Tank>
{ static constexpr int value = 60; };

template<class TypeTag>
struct GasUsage<TypeTag, TTag::Tank>
{ static constexpr int value = 65; };

template<class TypeTag>
struct CannonCaliber<TypeTag, TTag::Tank>
{ static constexpr int value = 120; };

template<class TypeTag>
struct GasUsage<TypeTag, TTag::Sedan>
{ static constexpr int value = 7; };

template<class TypeTag>
struct AutomaticTransmission<TypeTag, TTag::Sedan>
{ static constexpr bool value = true; };

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::Pickup>
{ static constexpr int value = 120; };

template<class TypeTag>
struct Payload<TypeTag, TTag::Pickup>
{ static constexpr int value = 5; };

template<class TypeTag>
struct TopSpeed<TypeTag, TTag::HummerH1>
{
    static constexpr int value
        = getPropValue<TypeTag, TTag::Pickup::TopSpeed<TypeTag>>();
};
```

## Property type alias version

The above hierarchy can also be written in a more terse notation using property type aliases.
For this to work, the properties _must_ be defined with the DUMUX_DEFINE_PROPERTY macro.

```cpp
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
```

## Retrieving properties

The property values can be retrieved with Dumux::getPropValue and some diagnostic messages can
be generated. For example

```cpp
int main()
{
    std::cout << "-- Top speed of sedan: " << getPropValue<Properties::TTag::Sedan, Properties::TopSpeed>() << "\n";
    std::cout << "-- Top speed of truck: " << getPropValue<Properties::TTag::Truck, Properties::TopSpeed>() << "\n";
}
```

will yield the following output:

```sh
-- Top speed of sedan: 210
-- Top speed of truck: 100
```
