// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Properties
 * \ingroup Typetraits
 * \author Timo Koch
 * \brief The Dumux property system, traits with inheritance
 */
#ifndef DUMUX_PROPERTY_SYSTEM_HH
#define DUMUX_PROPERTY_SYSTEM_HH

#include <tuple>
#include <type_traits>
#include <dune/common/std/type_traits.hh>

namespace Dumux::Properties {

/*!
 * \ingroup Properties
 * \brief a tag to mark properties as undefined
 */
struct UndefinedProperty {};

/*!
 * \ingroup Properties
 * \brief a tag to specify a direct alias for property extraction
 */
template<class P> struct PropertyAlias;

} // end namespace Dumux::Properties


// hide from doxygen
#ifndef DOXYGEN

//! implementation details for template meta programming
namespace Dumux::Properties::Detail {

//! check if a property P is defined
template<class P>
constexpr auto isDefinedProperty(int)
-> decltype(std::integral_constant<bool, !std::is_same_v<typename P::type, UndefinedProperty>>{})
{ return {}; }

//! fall back if a Property is defined
template<class P>
constexpr std::true_type isDefinedProperty(...) { return {}; }

//! check if a TypeTag inherits from other TypeTags
template<class T>
constexpr auto hasParentTypeTag(int)
-> decltype(std::declval<typename T::InheritsFrom>(), std::true_type{})
{ return {}; }

//! fall back if a TypeTag doesn't inherit
template<class T>
constexpr std::false_type hasParentTypeTag(...)
{ return {}; }


//! detect if the tag T has an alias with the same name as that of the property P
template<class P, class T>
using TypeAliasPropertyDetector = typename PropertyAlias<P>::template Alias<T>;

//! a property imitating the actual property P that extracts type and value
//! from alias members of tag T instead of from a property specialization for tag T
template<class P, class T>
struct TypeAliasProperty
{ using type = Dune::Std::detected_or_t<UndefinedProperty, TypeAliasPropertyDetector, P, T>; };

//! detect if the tag T has a template alias with the same name as that of the property P
template<class P, class T, class TypeTag>
using TemplateAliasPropertyDetector = typename PropertyAlias<P>::template TemplateAlias<T, TypeTag>;

//! a property imitating the actual property P that extracts type
//! from template alias members of tag T instead of from a property specialization for tag T
template<class P, class T, class TypeTag>
struct TemplateAliasProperty
{ using type = Dune::Std::detected_or_t<UndefinedProperty, TemplateAliasPropertyDetector, P, T, TypeTag>; };


//! detector for value members
template<class T>
using ValueMemberDetector = decltype(T::value);

//! specialization if there is no value member
template<class T, bool hasValueMember = Dune::Std::is_detected_v<ValueMemberDetector, T>>
struct ValueMember { static constexpr bool value = false; };

//! specialization if there is a value member
template<class T>
struct ValueMember<T, true> { static constexpr auto value = T::value; };

//! extract values from properties specializations
template<class PropertySpecialization>
struct GetPropValue { static constexpr auto value = ValueMember<PropertySpecialization>::value; };

//! extract values from type alias properties
template<class P, class T>
struct GetPropValue<TypeAliasProperty<P, T>>
{ static constexpr auto value = ValueMember<typename TypeAliasProperty<P, T>::type>::value; };

//! extract values from template alias properties
template<class P, class T, class TypeTag>
struct GetPropValue<TemplateAliasProperty<P, T, TypeTag>>
{ static constexpr auto value = ValueMember<typename TemplateAliasProperty<P, T, TypeTag>::type>::value; };


//! helper alias to concatenate multiple tuples
template<class ...Tuples>
using ConCatTuples = decltype(std::tuple_cat(std::declval<Tuples>()...));

//! helper struct to get the first property that is defined in the TypeTag hierarchy
template<class TypeTag, template<class,class> class Property, class TTagList>
struct GetDefined;

//! helper struct to iteratre over the TypeTag hierarchy
template<class TypeTag, template<class,class> class Property, class TTagList, class Enable>
struct GetNextTypeTag;

template<class TypeTag, template<class,class> class Property, class LastTypeTag>
struct GetNextTypeTag<TypeTag, Property, std::tuple<LastTypeTag>, std::enable_if_t<hasParentTypeTag<LastTypeTag>(int{}), void>>
{ using type = typename GetDefined<TypeTag, Property, typename LastTypeTag::InheritsFrom>::type; };

template<class TypeTag, template<class,class> class Property, class LastTypeTag>
struct GetNextTypeTag<TypeTag, Property, std::tuple<LastTypeTag>, std::enable_if_t<!hasParentTypeTag<LastTypeTag>(int{}), void>>
{ using type = UndefinedProperty; };

template<class TypeTag, template<class,class> class Property, class FirstTypeTag, class ...Args>
struct GetNextTypeTag<TypeTag, Property, std::tuple<FirstTypeTag, Args...>, std::enable_if_t<hasParentTypeTag<FirstTypeTag>(int{}), void>>
{ using type = typename GetDefined<TypeTag, Property, ConCatTuples<typename FirstTypeTag::InheritsFrom, std::tuple<Args...>>>::type; };

template<class TypeTag, template<class,class> class Property, class FirstTypeTag, class ...Args>
struct GetNextTypeTag<TypeTag, Property, std::tuple<FirstTypeTag, Args...>, std::enable_if_t<!hasParentTypeTag<FirstTypeTag>(int{}), void>>
{ using type = typename GetDefined<TypeTag, Property, std::tuple<Args...>>::type; };

template<class TypeTag, template<class,class> class Property, class LastTypeTag>
struct GetDefined<TypeTag, Property, std::tuple<LastTypeTag>>
{
// For clang, the following alias triggers compiler warnings if instantiated
// from something like `GetPropType<..., DeprecatedProperty>`, even if that is
// contained in a diagnostic pragma construct that should prevent these warnings.
// As a workaround, also add the pragmas around this line.
// See the discussion in MR 1647 for more details.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
    using LastType = Property<TypeTag, LastTypeTag>;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
    using DirectType = TypeAliasProperty<LastType, LastTypeTag>;
    using DirectTemplateType = TemplateAliasProperty<LastType, LastTypeTag, TypeTag>;
    // See below for an explanation of this
    using type = std::conditional_t<
        isDefinedProperty<LastType>(int{}), LastType,
        std::conditional_t<
            isDefinedProperty<DirectTemplateType>(int{}), DirectTemplateType,
            std::conditional_t<
                isDefinedProperty<DirectType>(int{}), DirectType,
                typename GetNextTypeTag<TypeTag, Property, std::tuple<LastTypeTag>, void>::type
            >
        >
    >;
};

template<class TypeTag, template<class,class> class Property, class FirstTypeTag, class ...Args>
struct GetDefined<TypeTag, Property, std::tuple<FirstTypeTag, Args...>>
{
// See the comment above.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
    using FirstType = Property<TypeTag, FirstTypeTag>;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
    using DirectType = TypeAliasProperty<FirstType, FirstTypeTag>;
    using DirectTemplateType = TemplateAliasProperty<FirstType, FirstTypeTag, TypeTag>;
    // First we check if the property is specialized for the current type tag
    // If yes, we found the correct specialization, if no we keep searching.
    // Second, we check if the type tag contains an alias with the property name.
    // If yes, we found the correct specialization, if no we keep searching.
    // Third, we check if the type tag contains a template alias with the property name (template argument is TypeTag)
    // If yes, we found the correct specialization, if no we the property is undefined (default definition)
    using type = std::conditional_t<
        isDefinedProperty<FirstType>(int{}), FirstType,
        std::conditional_t<
            isDefinedProperty<DirectTemplateType>(int{}), DirectTemplateType,
            std::conditional_t<
                isDefinedProperty<DirectType>(int{}), DirectType,
                typename GetNextTypeTag<TypeTag, Property, std::tuple<FirstTypeTag, Args...>, void>::type
            >
        >
    >;
};

//! helper struct to extract get the Property specialization given a TypeTag, asserts that the property is defined
template<class TypeTag, template<class,class> class Property>
struct GetPropImpl
{
    using type = typename Detail::GetDefined<TypeTag, Property, std::tuple<TypeTag>>::type;
    static_assert(!std::is_same_v<type, UndefinedProperty>, "Property is undefined!");
};

template<class TypeTag, template<class,class> class Property, class T>
struct GetPropOrImpl
{
    using PT = typename Detail::GetDefined<TypeTag, Property, std::tuple<TypeTag>>::type;
    struct WrapperT { using type = T; }; // fake property wrapper
    using type = std::conditional_t<std::is_same_v<PT, UndefinedProperty>, WrapperT, PT>;
};

template<class ParentTag, class TypeTag, bool hasParents = hasParentTypeTag<TypeTag>(int{})>
struct InheritsFrom;

template<class ParentTag, class TypeTag>
struct InheritsFrom<ParentTag, TypeTag, false> {
    static constexpr bool value = std::is_same_v<ParentTag, TypeTag>;
};

template<class ParentTag, class TypeTag>
struct InheritsFrom<ParentTag, TypeTag, true> {
    static constexpr bool value = std::is_same_v<ParentTag, TypeTag>
        || InheritsFrom<ParentTag, typename TypeTag::InheritsFrom, false>::value;
};

template<class ParentTag, class... TypeTags>
struct InheritsFrom<ParentTag, std::tuple<TypeTags...>, false> {
    static constexpr bool value = (InheritsFrom<ParentTag, TypeTags>::value || ...);
};

} // end namespace Dumux::Properties::Detail

#endif // DOXYGEN

namespace Dumux::Properties {

/*!
 * \ingroup Properties
 * \brief whether the property is defined/specialized for TypeTag
 */
template<class TypeTag, template<class,class> class Property>
inline constexpr bool hasDefinedType()
{
    using type = typename Detail::GetDefined<TypeTag, Property, std::tuple<TypeTag>>::type;
    return !std::is_same_v<type, UndefinedProperty>;
}

/*!
 * \ingroup Properties
 * \brief Return true if the given type tag inherits from the given parent type tag
 */
template<class ParentTypeTag, class TypeTag>
inline constexpr bool inheritsFrom()
{
    return Detail::InheritsFrom<ParentTypeTag, TypeTag>::value;
}

} // end namespace Dumux::Properties

namespace Dumux {

/*!
 * \ingroup Properties
 * \brief get the type of a property
 */
template<class TypeTag, template<class,class> class Property>
using GetProp = typename Properties::Detail::GetPropImpl<TypeTag, Property>::type;

/*!
 * \ingroup Properties
 * \brief get the type of a property or the type T if the property is undefined
 */
template<class TypeTag, template<class,class> class Property, class T>
using GetPropOr = typename Properties::Detail::GetPropOrImpl<TypeTag, Property, T>::type;

// See the comment above.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

/*!
 * \ingroup Properties
 * \brief get the type alias defined in the property
 */
template<class TypeTag, template<class,class> class Property>
using GetPropType = typename GetProp<TypeTag, Property>::type;

/*!
 * \ingroup Properties
 * \brief get the type alias defined in the property or the type T if the property is undefined
 */
template<class TypeTag, template<class,class> class Property, class T>
using GetPropTypeOr = typename GetPropOr<TypeTag, Property, T>::type;

/*!
 * \ingroup Properties
 * \brief get the value data member of a property
 */
template<class TypeTag, template<class,class> class Property>
inline constexpr auto getPropValue() { return Properties::Detail::GetPropValue<GetProp<TypeTag, Property>>::value; }
#ifdef __clang__
#pragma clang diagnostic pop
#endif

} // end namespace Dumux

/*!
 * \ingroup Properties
 * \brief A preprocessor macro to define properties
 * \note Every property can only be defined once (names have to be unique to the program)
 * \note Properties should be defined in the namespace Dumux::Properties
 * \details The macro defines two components for each property.
 * The first is the definition of the property. For example for
 * a property Scalar we get
 *
 * \code{.cpp}
   template<class TypeTag, class MyTypeTag> \
   struct Scalar { using type = UndefinedProperty; };
   \endcode
 *
 * The second is the specialization of the PropertyAlias template
 * for the newly defined property. For a property Scalar, we get
 *
 * \code{.cpp}
   template<class ...Args> // specialization for property Scalar
   struct PropertyAlias<Scalar<Args...>> {
       template <class MyTypeTag> using Alias = typename MyTypeTag::Scalar;
       template <class MyTypeTag, class TypeTag> using TemplateAlias = typename MyTypeTag::template Scalar<TypeTag>;
   };
   \endcode
 *
 * The specialization contains the template alias "Alias" that
 * can be used to check if a given type tag "MyTypeTag" has an alias member
 * "Scalar", and a template alias "TemplateAlias" that can be used to check
 * if a given type tag "MyTypeTag" has a template alias Scalar and can be
 * instantiated with a given type tag "TypeTag" (this will be the user-
 * end type tag).
 */
#define DUMUX_DEFINE_PROPERTY(Prop) \
    template<class TypeTag, class MyTypeTag> \
    struct Prop { \
        using type = UndefinedProperty; \
    }; \
    template<class ...Args> \
    struct PropertyAlias<Prop<Args...>> { \
        template <class MyTypeTag> using Alias = typename MyTypeTag::Prop; \
        template <class MyTypeTag, class TypeTag> using TemplateAlias = typename MyTypeTag::template Prop<TypeTag>; \
    };

#endif
