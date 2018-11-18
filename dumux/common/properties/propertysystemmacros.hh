// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Provides the magic behind the DuMuX property system.
 * \ingroup Properties
 *
 * Properties allow to associate arbitrary data types to
 * identifiers. A property is always defined on a pair (TypeTag,
 * PropertyTag) where TypeTag is the identifier for the object the
 * property is defined for and PropertyTag is an unique identifier of
 * the property.
 *
 * Type tags are hierarchic and inherit properties defined on their
 * ancesters. At each level, properties defined on lower levels can be
 * overwritten or even made undefined. It is also possible to define
 * defaults for properties if it makes sense.
 *
 * Properties may make use other properties for the respective type
 * tag and these properties can also be defined on an arbitrary level
 * of the hierarchy.
 */
#ifndef DUMUX_PROPERTY_SYSTEM_MACROS_HH
#define DUMUX_PROPERTY_SYSTEM_MACROS_HH

#include <dumux/common/properties/propertysystem.hh>

namespace Dumux {
namespace Properties {

namespace TTag {}

/*!
 * \brief Makes a type out of a type tag name
 */
#define TTAG(TypeTagName) ::Dumux::Properties::TTag::TypeTagName

/*!
 * \brief Makes a type out of a property tag name
 */
//#define PTAG(PropTagName) PropTagName

/*!
 * \brief Makes a type out of a property tag name
 */
#define PTAG_(PropTagName) ::Dumux::Properties::PropTagName


// in the old property system the order in inherit_from was the other way around
// this flips the order of a tuple to restore old behaviour when using the macro.
// when you are using non-macro version make sure to flip the order.
template<class Tuple, class IndexSequence>
struct ReverseTupleImpl;

template<class Tuple, size_t... I>
struct ReverseTupleImpl<Tuple, std::index_sequence<I...>>
{
  using type = std::tuple<std::tuple_element_t<sizeof...(I) - 1 - I, Tuple>...>;
};

// revert tuple argument order
template<class Tuple>
using ReverseTuple = typename ReverseTupleImpl<Tuple, std::make_index_sequence<std::tuple_size<Tuple>::value>>::type;

// a temporary hack to make the macro still work, we set using InheritsFrom = void,
// which gets picked up by the new property as non inheritance, this can be removed
// once all macros are gone
namespace Detail {
template<class TypeTagTuple>
struct GetTypeTagInheritance;

template<class OneTypeTag>
struct GetTypeTagInheritance<std::tuple<OneTypeTag>>
{
    using type = void;
};

template<class FirstTypeTag, class ...OtherTypeTags>
struct GetTypeTagInheritance<std::tuple<FirstTypeTag, OtherTypeTags...>>
{
    // reverse order to restore old behaviour
    using type = ReverseTuple<std::tuple<OtherTypeTags...>>;
};
} // end namespace Detail

/*!
 * \brief Define a new type tag.
 *
 * A type tag can inherit the properties defined on up to five parent
 * type tags. Examples:
 *
 * \code
 * // The type tag doesn't inherit any properties from other type tags
 * NEW_TYPE_TAG(FooTypeTag);
 *
 * // BarTypeTag inherits all properties from FooTypeTag
 * NEW_TYPE_TAG(BarTypeTag, INHERITS_FROM(FooTypeTag));
 *
 * // FooBarTypeTag inherits the properties of FooTypeTag as well as
 * // those of BarTypeTag. Properties defined on BarTypeTag have
 * // preceedence over those defined for FooTypeTag:
 * NEW_TYPE_TAG(FooBarTypeTag, INHERITS_FROM(FooTypeTag, BarTypeTag));
 * \endcode
 */
#define DUMUX_GET_HEAD_(Arg1, ...) Arg1

#define NEW_TYPE_TAG(...)                                                          \
    namespace TTag {                                                               \
    struct DUMUX_GET_HEAD_(__VA_ARGS__)                                            \
    { using InheritsFrom = Detail::GetTypeTagInheritance<std::tuple<__VA_ARGS__>>::type; };   \
    } extern int semicolonHack_

/*!
 * \brief Syntactic sugar for NEW_TYPE_TAG.
 *
 * See the documentation for NEW_TYPE_TAG.
 */
#define INHERITS_FROM(...) __VA_ARGS__

/*!
 * \brief Define a property tag.
 *
 * A property tag is the unique identifier for a property. It may only
 * be declared once in your program. There is also no hierarchy of
 * property tags as for type tags.
 *
 * Examples:
 *
 * \code
 * NEW_PROP_TAG(blubbPropTag);
 * NEW_PROP_TAG(blabbPropTag);
 * \endcode
 */
#define NEW_PROP_TAG(PTagName)                             \
    template<class TypeTag, class MyTypeTag>               \
    struct PTagName { using type = UndefinedProperty; };   \
    extern int semicolonHack_

/*!
 * \brief Set a property for a specific type tag.
 *
 * After this macro, you must to specify a complete body of a class
 * template, including the trailing semicolon. If you need to retrieve
 * another property within the class body, you can use TypeTag as the
 * argument for the type tag for the GET_PROP macro.
 *
 * Example:
 *
 * \code
 * SET_PROP(FooTypeTag, blubbPropTag)
 * {
 *    static int value = 10;
 *    static int calculate(int arg)
 *    { calculateInternal_(arg); }
 *
 * private:
 *    // retrieve the blabbProp property for the real TypeTag the
 *    // property is defined on. Note that blabbProb does not need to
 *    // be defined on FooTypeTag, but can also be defined for some
 *    // derived type tag.
 *    using blabb = typename GET_PROP(TypeTag, blabbProp);
 *
 *    static int calculateInternal_(int arg)
 *    { return arg * blabb::value; };
 * \endcode
 * };
 */
#define SET_PROP(EffTypeTagName, PropTagName)                   \
    template <class TypeTag>                                    \
    struct PropTagName<TypeTag, TTAG(EffTypeTagName)>

/*!
 * \brief Set a property to a simple constant integer value.
 *
 * The constant can be accessed by the 'value' attribute.
 */
#define SET_INT_PROP(EffTypeTagName, PropTagName, /*Value*/...)    \
    template <class TypeTag>                                    \
    struct PropTagName<TypeTag, TTAG(EffTypeTagName)>     \
    {                                                           \
        using type = int;                                       \
        static constexpr int value = __VA_ARGS__;               \
    }

/*!
 * \brief Set a property to a simple constant boolean value.
 *
 * The constant can be accessed by the 'value' attribute.
 */
#define SET_BOOL_PROP(EffTypeTagName, PropTagName, /*Value*/...)    \
    template <class TypeTag>                                    \
    struct PropTagName<TypeTag, TTAG(EffTypeTagName)>     \
    {                                                           \
        using type = bool;                                       \
        static constexpr bool value = __VA_ARGS__;               \
    }

/*!
 * \brief Set a property which defines a type.
 *
 * The type can be accessed by the 'type' attribute.
 */
#define SET_TYPE_PROP(EffTypeTagName, PropTagName, /*Value*/...)  \
    template <class TypeTag>                                      \
    struct PropTagName<TypeTag, TTAG(EffTypeTagName)>             \
    {                                                             \
        using type = __VA_ARGS__;                                 \
    }

/*!
 * \brief Set a property to a simple constant scalar value.
 *
 * The constant can be accessed by the 'value' attribute. In order to
 * use this macro, the property tag "Scalar" needs to be defined for
 * the real type tag.
 */
#define SET_SCALAR_PROP(EffTypeTagName, PropTagName, ...)               \
    template <class TypeTag>                                    \
    struct PropTagName<TypeTag, TTAG(EffTypeTagName)>     \
    {                                                           \
        using Scalar = Dumux::GetPropType<TypeTag, Scalar>;            \
    public:                                                             \
        using type = Scalar;                                            \
        static const Scalar value;                                      \
    };                                                                  \
    template <class TypeTag>                                            \
    const typename PropTagName<TypeTag, TTAG(EffTypeTagName)>::type   \
    PropTagName<TypeTag, TTAG(EffTypeTagName)>::value(__VA_ARGS__)

/*!
 * \brief Set a property to a simple constant string value.
 *
 * The constant can be accessed by the 'value' attribute and is of
 * type std::string.
 */
#define SET_STRING_PROP(EffTypeTagName, PropTagName, ...)               \
    template <class TypeTag>                                    \
    struct PropTagName<TypeTag, TTAG(EffTypeTagName)>     \
    {                                                           \
    public:                                                             \
        using type = std::string;                                            \
        static const std::string value;                                      \
    };                                                                  \
    template <class TypeTag>                                            \
    const typename PropTagName<TypeTag, TTAG(EffTypeTagName)>::type   \
    PropTagName<TypeTag, TTAG(EffTypeTagName)>::value(__VA_ARGS__)


// getters
#define GET_PROP(TypeTag, PropTagName) ::Dumux::Properties::Detail::GetPropImpl<TypeTag, PTAG_(PropTagName)>::type
#define GET_PROP_VALUE(TypeTag, PropTagName) ::Dumux::Properties::Detail::GetPropImpl<TypeTag, PTAG_(PropTagName)>::type::value
#define GET_PROP_TYPE(TypeTag, PropTagName) ::Dumux::Properties::Detail::GetPropImpl<TypeTag, PTAG_(PropTagName)>::type::type

} // namespace Properties
} // namespace Dumux

#endif // DUMUX_PROPERTY_SYSTEM_HH
