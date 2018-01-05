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
#ifndef DUMUX_PROPERTY_SYSTEM_HH
#define DUMUX_PROPERTY_SYSTEM_HH

#include <cstring>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
// For is_base_of
#include <type_traits>

// Integral Constant Expressions
#include <dune/common/classname.hh>

namespace Dumux
{
namespace Properties
{
#define DUMUX_GET_HEAD_(Arg1, ...) Arg1
#define DUMUX_GET_TAIL_(Arg1, ...) Blubber // __VA_ARGS__

#define DUMUX_STRINGIGY_HEAD_(Arg1, ...) #Arg1
#define DUMUX_STRINGIGY_TAIL_(Arg1, ...) #__VA_ARGS__

#if !defined NO_PROPERTY_INTROSPECTION

//! Internal macro which is only required if the property introspection is enabled
#define PROP_INFO_(EffTypeTagName, PropKind, PropTagName, ...)          \
    template <>                                                         \
    struct PropertyInfo<TTAG(EffTypeTagName), PTAG_(PropTagName)>       \
    {                                                                   \
    static int init() {                                                 \
        PropertyRegistryKey key(                                        \
            /*effTypeTagName=*/ Dune::className<TTAG(EffTypeTagName)>(), \
            /*kind=*/PropKind,                                          \
            /*name=*/#PropTagName,                                      \
            /*value=*/#__VA_ARGS__,                                     \
            /*file=*/__FILE__,                                          \
            /*line=*/__LINE__);                                         \
        PropertyRegistry::addKey(key);                                  \
        return 0;                                                       \
    };                                                                  \
    static int foo;                                                     \
    };                                                                  \
    int PropertyInfo<TTAG(EffTypeTagName), PTAG_(PropTagName)>::foo =   \
    PropertyInfo<TTAG(EffTypeTagName), PTAG_(PropTagName)>::init();

//! Internal macro which is only required if the property introspection is enabled
#define TTAG_INFO_(...)                                                 \
    template <>                                                         \
    struct TypeTagInfo<DUMUX_GET_HEAD_(__VA_ARGS__)>               \
    {                                                                   \
        static int init() {                                             \
            TypeTagRegistry::addChildren<__VA_ARGS__>();                \
            return 0;                                                   \
        }                                                              \
        static int foo;                                                 \
    };                                                                  \
    int TypeTagInfo<DUMUX_GET_HEAD_(__VA_ARGS__)>::foo =           \
        TypeTagInfo<DUMUX_GET_HEAD_(__VA_ARGS__)>::init();

#else
//! Don't do anything if introspection is disabled
#define PROP_INFO_(EffTypeTagName, PropKind, PropTagName, ...)
#define TTAG_INFO_(EffTypeTagName, ...)
#endif

// some macros for simplification

/*!
 * \brief Makes a type out of a type tag name
 */
#define TTAG(TypeTagName) ::Dumux::Properties::TTag::TypeTagName

/*!
 * \brief Makes a type out of a property tag name
 */
#define PTAG(PropTagName) PropTagName

/*!
 * \brief Makes a type out of a property tag name
 */
#define PTAG_(PropTagName) ::Dumux::Properties::PTag::PropTagName

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
#define NEW_TYPE_TAG(...)                                               \
    namespace TTag {                                                    \
    struct DUMUX_GET_HEAD_(__VA_ARGS__, blubb)                          \
        : public TypeTag<__VA_ARGS__>                                   \
    { };                                                                \
    TTAG_INFO_(__VA_ARGS__, void)                                       \
    }                                                                   \
    extern int semicolonHack_

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
    namespace PTag {                                       \
    struct PTagName; } extern int semicolonHack_

//! Internal macro
#define SET_PROP_(EffTypeTagName, PropKind, PropTagName, ...)       \
    template <class TypeTag>                                        \
    struct Property<TypeTag,                                        \
                    TTAG(EffTypeTagName),                           \
                    PTAG_(PropTagName)>;                            \
    PROP_INFO_(EffTypeTagName,                                      \
               /*kind=*/PropKind,                                   \
               PropTagName,                                         \
               /*value=*/__VA_ARGS__)                               \
    template <class TypeTag>                                        \
    struct Property<TypeTag,                                        \
                    TTAG(EffTypeTagName),                           \
                    PTAG_(PropTagName) >

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
    struct Property<TypeTag,                                    \
                    TTAG(EffTypeTagName),                       \
                    PTAG_(PropTagName)>;                        \
    PROP_INFO_(EffTypeTagName,                                  \
               /*kind=*/"opaque",                               \
               PropTagName,                                     \
               /*value=*/"<opaque>")                            \
    template <class TypeTag>                                    \
    struct Property<TypeTag,                                    \
                    TTAG(EffTypeTagName),                       \
                    PTAG_(PropTagName) >

/*!
 * \brief Explicitly unset a property for a type tag.
 *
 * This means that the property will not be inherited from the type
 * tag's parents and that no default will be used.
 *
 * Example:
 *
 * \code
 * // make the blabbPropTag property undefined for the BarTypeTag.
 * UNSET_PROP(BarTypeTag, blabbPropTag);
 * \endcode
 */
#define UNSET_PROP(EffTypeTagName, PropTagName)                 \
    template <>                                                 \
    struct PropertyUnset<TTAG(EffTypeTagName),                  \
                         PTAG_(PropTagName) >;                  \
    PROP_INFO_(EffTypeTagName,                                  \
               /*kind=*/"withdraw",                             \
               PropTagName,                                     \
               /*value=*/<none>)                                \
    template <>                                                 \
    struct PropertyUnset<TTAG(EffTypeTagName),                  \
                         PTAG_(PropTagName) >                   \
        : public PropertyExplicitlyUnset                        \
        {}

/*!
 * \brief Set a property to a simple constant integer value.
 *
 * The constant can be accessed by the 'value' attribute.
 */
#define SET_INT_PROP(EffTypeTagName, PropTagName, /*Value*/...)    \
    SET_PROP_(EffTypeTagName,                                   \
              /*kind=*/"int   ",                                \
              PropTagName,                                      \
              /*value=*/__VA_ARGS__)                            \
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
    SET_PROP_(EffTypeTagName,                                       \
              /*kind=*/"bool  ",                                    \
              PropTagName,                                          \
              /*value=*/__VA_ARGS__)                                \
    {                                                               \
        using type = bool;                                          \
        static constexpr bool value = __VA_ARGS__;                  \
    }

/*!
 * \brief Set a property which defines a type.
 *
 * The type can be accessed by the 'type' attribute.
 */
#define SET_TYPE_PROP(EffTypeTagName, PropTagName, /*Value*/...)  \
    SET_PROP_(EffTypeTagName,                                     \
              /*kind=*/"type  ",                                  \
              PropTagName,                                        \
              /*value=*/__VA_ARGS__)                              \
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
    SET_PROP_(EffTypeTagName,                                           \
              /*kind=*/"scalar",                                        \
              PropTagName,                                              \
              /*value=*/__VA_ARGS__)                                    \
    {                                                                   \
        using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);         \
    public:                                                             \
        using type = Scalar;                                            \
        static const Scalar value;                                      \
    };                                                                  \
    template <class TypeTag>                                            \
    const typename Property<TypeTag, TTAG(EffTypeTagName), PTAG_(PropTagName)>::type   \
    Property<TypeTag, TTAG(EffTypeTagName), PTAG_(PropTagName)>::value(__VA_ARGS__)

/*!
 * \brief Set a property to a simple constant string value.
 *
 * The constant can be accessed by the 'value' attribute and is of
 * type std::string.
 */
#define SET_STRING_PROP(EffTypeTagName, PropTagName, ...)               \
    SET_PROP_(EffTypeTagName,                                           \
              /*kind=*/"string",                                        \
              PropTagName,                                              \
              /*value=*/__VA_ARGS__)                                    \
    {                                                                   \
    public:                                                             \
        using type = std::string;                                       \
        static const std::string value;                                 \
    };                                                                  \
    template <class TypeTag>                                            \
    const typename Property<TypeTag, TTAG(EffTypeTagName), PTAG_(PropTagName)>::type \
    Property<TypeTag, TTAG(EffTypeTagName), PTAG_(PropTagName)>::value(__VA_ARGS__)

/*!
 * \brief Get the property for a type tag.
 *
 * If you use GET_PROP within a template and want to refer to some
 * type (including the property itself), GET_PROP must be preceeded by
 * the 'typename' keyword.
 */
#define GET_PROP(TypeTag, PropTagName) \
    ::Dumux::Properties::GetProperty<TypeTag, PTAG_(PropTagName)>::p
#define GET_PROP_(TypeTag, PropTag) \
    ::Dumux::Properties::GetProperty<TypeTag, PropTag>::p

/*!
 * \brief Access the 'value' attribute of a property for a type tag.
 *
 * This is just for convenience and equivalent to GET_PROP(TypeTag,
 * PropTag) :: value.  If the property doesn't have an attribute named
 * 'value', this yields a compiler error.
 */
#define GET_PROP_VALUE(TypeTag, PropTagName)                            \
    ::Dumux::Properties::GetProperty<TypeTag, PTAG_(PropTagName)>::p::value
#define GET_PROP_VALUE_(TypeTag, PropTag)                               \
    ::Dumux::Properties::GetProperty<TypeTag, PropTag>::p::value

/*!
 * \brief Access the 'type' attribute of a property for a type tag.
 *
 * This is just for convenience and equivalent to GET_PROP(TypeTag,
 * PropTag) :: type.  If the property doesn't have an attribute named
 * 'type', this yields a compiler error. Also, if you use this macro
 * within a template, it must be preceeded by the 'typename' keyword.
 */
#define GET_PROP_TYPE(TypeTag, PropTagName) \
    ::Dumux::Properties::GetProperty<TypeTag, PTAG_(PropTagName)>::p::type
#define GET_PROP_TYPE_(TypeTag, PropTag) \
    ::Dumux::Properties::GetProperty<TypeTag, PropTag>::p::type

#if !defined NO_PROPERTY_INTROSPECTION
/*!
 * \brief Return a human readable diagnostic message how exactly a
 *        property was defined.
 *
 * This is only enabled if the NO_PROPERTY_INTROSPECTION macro is not
 * defined.
 *
 * Example:
 *
 * \code
 * int main()
 * {
 *    std::cout << PROP_DIAGNOSTIC(FooBarTypeTag, blabbPropTag) << "\n";
 * };
 * \endcode
 */
#define PROP_DIAGNOSTIC(TypeTag, PropTagName) \
    ::Dumux::Properties::getDiagnostic<TypeTag>(#PropTagName)

#else
/*!
 * \brief Return a human readable diagnostic message how exactly a
 *        property was defined.
 *
 * This is only enabled if the NO_PROPERTY_INTROSPECTION macro is not
 * defined.
 *
 * Example:
 *
 * \code
 * int main()
 * {
 *    std::cout << PROP_DIAGNOSTIC(FooBarTypeTag, blabbPropTag) << "\n";
 * };
 * \endcode
 */
#define PROP_DIAGNOSTIC(TypeTag, PropTagName) "Property introspection disabled by NO_PROPERTY_INTROSPECTION"
#endif


//////////////////////////////////////////////
// some serious template kung fu. Don't look at it too closely, it
// might damage your brain!
//////////////////////////////////////////////

//! \cond false

namespace PTag {}
namespace TTag {}

#if !defined NO_PROPERTY_INTROSPECTION

namespace TTag
{
template <class EffTypeTag>
struct TypeTagInfo
{};
}

template <class EffTypeTagName, class PropTagName>
struct PropertyInfo
{};
class PropertyRegistryKey
{
public:
    PropertyRegistryKey()
    {};

    PropertyRegistryKey(const std::string &effTypeTagName,
                        const std::string &propertyKind,
                        const std::string &propertyName,
                        const std::string &propertyValue,
                        const std::string &fileDefined,
                        int lineDefined)
        : effTypeTagName_(effTypeTagName)
        , propertyKind_(propertyKind)
        , propertyName_(propertyName)
        , propertyValue_(propertyValue)
        , fileDefined_(fileDefined)
        , lineDefined_(lineDefined)
    {
    };

    // copy constructor
    PropertyRegistryKey(const PropertyRegistryKey &v)
        : effTypeTagName_(v.effTypeTagName_)
        , propertyKind_(v.propertyKind_)
        , propertyName_(v.propertyName_)
        , propertyValue_(v.propertyValue_)
        , fileDefined_(v.fileDefined_)
        , lineDefined_(v.lineDefined_)
    {};

    const std::string &effTypeTagName() const
    { return effTypeTagName_; }
    const std::string &propertyKind() const
    { return propertyKind_; }
    const std::string &propertyName() const
    { return propertyName_; }
    const std::string &propertyValue() const
    { return propertyValue_; }
    const std::string &fileDefined() const
    { return fileDefined_; }
    int lineDefined() const
    { return lineDefined_; }

private:
    std::string effTypeTagName_;
    std::string propertyKind_;
    std::string propertyName_;
    std::string propertyValue_;
    std::string fileDefined_;
    int lineDefined_;
};

class PropertyRegistry
{
public:
    using KeyList = std::map<std::string, PropertyRegistryKey>;
    using KeyListMap = std::map<std::string, KeyList>;

    static void addKey(const PropertyRegistryKey &key)
    {
        keys_[key.effTypeTagName()][key.propertyName()] = key;
    }

    static const PropertyRegistryKey &getKey(const std::string &effTypeTagName,
                                             const std::string &propertyName)
    {
        return keys_[effTypeTagName][propertyName];
    };

    static const KeyList &getKeys(const std::string &effTypeTagName)
    {
        return keys_[effTypeTagName];
    };

private:
    static KeyListMap keys_;
};
PropertyRegistry::KeyListMap PropertyRegistry::keys_;

class TypeTagRegistry
{
public:
    using ChildrenList = std::list<std::string>;
    using ChildrenListMap = std::map<std::string, ChildrenList>;

    template <class TypeTag, class Child1, class Child2, class Child3, class Child4, class Child5, class Dummy>
    static void addChildren()
    {
        std::string typeTagName = Dune::className<TypeTag>();
        if (typeid(Child1) != typeid(void))
            keys_[typeTagName].push_front(Dune::className<Child1>());
        if (typeid(Child2) != typeid(void))
            keys_[typeTagName].push_front(Dune::className<Child2>());
        if (typeid(Child3) != typeid(void))
            keys_[typeTagName].push_front(Dune::className<Child3>());
        if (typeid(Child4) != typeid(void))
            keys_[typeTagName].push_front(Dune::className<Child4>());
        if (typeid(Child5) != typeid(void))
            keys_[typeTagName].push_front(Dune::className<Child5>());
    }

    template <class TypeTag, class Child1, class Child2, class Child3, class Child4, class Dummy>
    static void addChildren()
    { addChildren<TypeTag, Child1, Child2, Child3, Child4, void, Dummy>(); }

    template <class TypeTag, class Child1, class Child2, class Child3, class Dummy>
    static void addChildren()
    { addChildren<TypeTag, Child1, Child2, Child3, void, void, Dummy>(); }

    template <class TypeTag, class Child1, class Child2, class Dummy>
    static void addChildren()
    { addChildren<TypeTag, Child1, Child2, void, void,  void, Dummy>(); }

    template <class TypeTag, class Child1, class Dummy>
    static void addChildren()
    { addChildren<TypeTag, Child1, void, void, void, void, Dummy>(); }

    template <class TypeTag, class Dummy>
    static void addChildren()
    { addChildren<TypeTag, void, void, void, void, void, Dummy>(); }

    static const ChildrenList &children(const std::string &typeTagName)
    {
        return keys_[typeTagName];
    };

private:
    static ChildrenListMap keys_;
};

TypeTagRegistry::ChildrenListMap TypeTagRegistry::keys_;

#endif // !defined NO_PROPERTY_INTROSPECTION

using std::is_void;
using std::is_base_of;

// logical AND, OR and NOT operations to be used for template meta programming
template <bool b1, bool b2, bool b3 = true, bool b4 = true, bool b5 = true, bool b6 = true, bool b7 = true>
struct ice_and
{
    static const bool value = false;
};

template <>
struct ice_and<true, true, true, true, true, true, true>
{
    static const bool value = true;
};

template <bool b1, bool b2, bool b3 = false, bool b4 = false, bool b5 = false, bool b6 = false, bool b7 = false>
struct ice_or
{
    static const bool value = true;
};

template <>
struct ice_or<false, false, false, false, false, false, false>
{
    static const bool value = false;
};

template <bool b>
struct ice_not
{
    static const bool value = false;
};

template <>
struct ice_not<false>
{
    static const bool value = true;
};

//! \internal
class PropertyUndefined {};
//! \internal
class PropertyExplicitlyUnset {};

//! \internal
template <class RealTypeTag,
          class EffectiveTypeTag,
          class PropertyTag>
struct Property : public PropertyUndefined
{
};

//! \internal
template <class EffectiveTypeTag,
          class PropertyTag>
struct PropertyUnset : public PropertyUndefined
{
};

//! \internal
template <class RealTypeTag,
          class PropertyTag>
struct DefaultProperty : public PropertyUndefined
{
};

//! \internal
template <class Tree, class PropertyTag>
struct propertyExplicitlyUnset
{
    const static bool value =
        is_base_of<PropertyExplicitlyUnset,
                   PropertyUnset<typename Tree::SelfType,
                                 PropertyTag>
                   >::value;
};

//! \internal
template <class Tree, class PropertyTag>
class propertyExplicitlyUnsetOnTree
{
    static const bool explicitlyUnset = propertyExplicitlyUnset<Tree, PropertyTag>::value;

    static const bool isLeaf = ice_and<is_void<typename Tree::Child1>::value,
                                       is_void<typename Tree::Child2>::value,
                                       is_void<typename Tree::Child3>::value,
                                       is_void<typename Tree::Child4>::value,
                                       is_void<typename Tree::Child5>::value >::value;

public:
    static const bool value =
        ice_or<explicitlyUnset,
               ice_and<ice_not<isLeaf>::value,
                       propertyExplicitlyUnsetOnTree<typename Tree::Child1, PropertyTag>::value,
                       propertyExplicitlyUnsetOnTree<typename Tree::Child2, PropertyTag>::value,
                       propertyExplicitlyUnsetOnTree<typename Tree::Child3, PropertyTag>::value,
                       propertyExplicitlyUnsetOnTree<typename Tree::Child4, PropertyTag>::value,
                       propertyExplicitlyUnsetOnTree<typename Tree::Child5, PropertyTag>::value
                       >::value
               >::value;
};

//! \internal
template <class PropertyTag>
struct propertyExplicitlyUnsetOnTree<void, PropertyTag>
{
    const static bool value = std::true_type::value;
};

//! \internal
template <class RealTypeTag, class Tree, class PropertyTag>
struct propertyDefinedOnSelf
{
    const static bool value =
        ice_not<is_base_of<PropertyUndefined,
                           Property<RealTypeTag,
                                    typename Tree::SelfType,
                                    PropertyTag>
                           >::value
                           >::value;
};

//! \internal
template <class RealTypeTag, class Tree, class PropertyTag>
class propertyDefinedOnTree
{
    static const bool notExplicitlyUnset =
        ice_not<propertyExplicitlyUnsetOnTree<Tree,
                                              PropertyTag>::value >::value;

public:
    static const bool value =
        ice_and<notExplicitlyUnset,
                ice_or<propertyDefinedOnSelf<RealTypeTag, Tree, PropertyTag>::value,
                       propertyDefinedOnTree<RealTypeTag, typename Tree::Child1, PropertyTag>::value,
                       propertyDefinedOnTree<RealTypeTag, typename Tree::Child2, PropertyTag>::value,
                       propertyDefinedOnTree<RealTypeTag, typename Tree::Child3, PropertyTag>::value,
                       propertyDefinedOnTree<RealTypeTag, typename Tree::Child4, PropertyTag>::value,
                       propertyDefinedOnTree<RealTypeTag, typename Tree::Child5, PropertyTag>::value
                       >::value >::value;
};

//! \internal
template <class RealTypeTag, class PropertyTag>
class propertyDefinedOnTree<RealTypeTag, void, PropertyTag>
{
public:
    static const bool value = std::false_type::value;
};

//! \internal
template <class RealTypeTag, class PropertyTag>
struct defaultPropertyDefined
{
    const static bool value =
        ice_not<is_base_of<PropertyUndefined,
                           DefaultProperty<RealTypeTag,
                                           PropertyTag>
                           >::value
                           >::value;
};

//! \internal
template <class RealTypeTag, class Tree, class PropertyTag>
class defaultPropertyDefinedOnTree
{
    static const bool isLeaf = ice_and<is_void<typename Tree::Child1>::value,
                                       is_void<typename Tree::Child2>::value,
                                       is_void<typename Tree::Child3>::value,
                                       is_void<typename Tree::Child4>::value,
                                       is_void<typename Tree::Child5>::value >::value;

    static const bool explicitlyUnset =
        propertyExplicitlyUnsetOnTree<Tree, PropertyTag>::value;

public:
    static const bool value =
        ice_and<ice_not<explicitlyUnset>::value,
                ice_or<ice_and<isLeaf,defaultPropertyDefined<RealTypeTag, PropertyTag>::value >::value,
                       defaultPropertyDefinedOnTree<RealTypeTag,typename Tree::Child1, PropertyTag>::value,
                       defaultPropertyDefinedOnTree<RealTypeTag,typename Tree::Child2, PropertyTag>::value,
                       defaultPropertyDefinedOnTree<RealTypeTag,typename Tree::Child3, PropertyTag>::value,
                       defaultPropertyDefinedOnTree<RealTypeTag,typename Tree::Child4, PropertyTag>::value,
                       defaultPropertyDefinedOnTree<RealTypeTag,typename Tree::Child5, PropertyTag>::value
                       >::value >::value;
};

//! \internal
template <class RealTypeTag, class PropertyTag>
struct defaultPropertyDefinedOnTree<RealTypeTag,void, PropertyTag>
{
    static const bool value = std::false_type::value;
};

//! \internal
template <class RealTypeTag, class Tree, class PropertyTag>
class propertyDefined
{
public:
    static const bool onSelf = propertyDefinedOnSelf<RealTypeTag,Tree,PropertyTag>::value;

    static const bool onChild1 = propertyDefinedOnTree<RealTypeTag,typename Tree::Child1,PropertyTag>::value;
    static const bool onChild2 = propertyDefinedOnTree<RealTypeTag,typename Tree::Child2,PropertyTag>::value;
    static const bool onChild3 = propertyDefinedOnTree<RealTypeTag,typename Tree::Child3,PropertyTag>::value;
    static const bool onChild4 = propertyDefinedOnTree<RealTypeTag,typename Tree::Child4,PropertyTag>::value;
    static const bool onChild5 = propertyDefinedOnTree<RealTypeTag,typename Tree::Child5,PropertyTag>::value;

    static const bool asDefault =
        defaultPropertyDefinedOnTree<RealTypeTag, Tree,PropertyTag>::value;

    static const bool onChildren =
        ice_or<onChild1,
               onChild2,
               onChild3,
               onChild4,
               onChild5
               >::value;

    static const bool value =
        ice_or<onSelf ,
               onChildren>::value;


};

//! \internal
template <class RealTypeTag, class Tree, class PropertyTag>
class propertyTagIndex
{
    using definedWhere = propertyDefined<RealTypeTag, Tree, PropertyTag>;

public:
    static const int value =
        definedWhere::onSelf ? 0 :
        ( definedWhere::onChild5 ? 5 :
          ( definedWhere::onChild4 ? 4 :
            ( definedWhere::onChild3 ? 3 :
              ( definedWhere::onChild2 ? 2 :
                ( definedWhere::onChild1 ? 1 :
                  ( definedWhere::asDefault ? -1 :
                    -1000))))));
};


//! \internal
template <class SelfT,
          class Child1T = void,
          class Child2T = void,
          class Child3T = void,
          class Child4T = void,
          class Child5T = void>
class TypeTag
{
public:
    using SelfType = SelfT;

    using Child1 = Child1T;
    using Child2 = Child2T;
    using Child3 = Child3T;
    using Child4 = Child4T;
    using Child5 = Child5T;
};

NEW_TYPE_TAG(__Default);

//! \internal
template <class EffectiveTypeTag,
          class PropertyTag,
          class RealTypeTag=EffectiveTypeTag,
          int tagIndex = propertyTagIndex<RealTypeTag, EffectiveTypeTag, PropertyTag>::value >
struct GetProperty
{
};

// property not defined, but a default property is available
//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, -1>
{
    using p = DefaultProperty<RealTypeTag, PropertyTag>;
};

// property defined on self
//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 0>
{
    using p = Property<RealTypeTag, TypeTag, PropertyTag>;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 1>
{
    using p = typename GetProperty<typename TypeTag::Child1, PropertyTag, RealTypeTag>::p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 2>
{
    using p = typename GetProperty<typename TypeTag::Child2, PropertyTag, RealTypeTag>::p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 3>
{
    using p = typename GetProperty<typename TypeTag::Child3, PropertyTag, RealTypeTag>::p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 4>
{
    using p = typename GetProperty<typename TypeTag::Child4, PropertyTag, RealTypeTag>::p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 5>
{
    using p = typename GetProperty<typename TypeTag::Child5, PropertyTag, RealTypeTag>::p;
};

#if !defined NO_PROPERTY_INTROSPECTION
std::string canonicalTypeTagNameToName_(const std::string &canonicalName)
{
    std::string result(canonicalName);
    result.replace(0, strlen("Dumux::Properties::TTag::"), "");
    return result;
}

inline bool getDiagnostic_(const std::string &typeTagName,
                           const std::string &propTagName,
                           std::string &result,
                           const std::string &indent)
{
    const PropertyRegistryKey *key = 0;

    const PropertyRegistry::KeyList &keys =
        PropertyRegistry::getKeys(typeTagName);
    PropertyRegistry::KeyList::const_iterator it = keys.begin();
    for (; it != keys.end(); ++it) {
        if (it->second.propertyName() == propTagName) {
            key = &it->second;
            break;
        }
    }

    if (key) {
        std::ostringstream oss;
        oss << indent
            << key->propertyKind() << " "
            << key->propertyName() << " defined on '"
            << canonicalTypeTagNameToName_(key->effTypeTagName()) << "' at "
            << key->fileDefined() << ":" << key->lineDefined() << "\n";
        result = oss.str();
        return true;
    }

    // print properties defined on children
    using ChildrenList = TypeTagRegistry::ChildrenList;
    const ChildrenList &children = TypeTagRegistry::children(typeTagName);
    ChildrenList::const_iterator ttagIt = children.begin();
    std::string newIndent = indent + "  ";
    for (; ttagIt != children.end(); ++ttagIt) {
        if (getDiagnostic_(*ttagIt, propTagName, result, newIndent)) {
            result.insert(0, indent + "Inherited from " + canonicalTypeTagNameToName_(typeTagName) + "\n");
            return true;
        }
    }

    return false;
}

template <class TypeTag>
const std::string getDiagnostic(std::string propTagName)
{
    std::string result;

    std::string TypeTagName(Dune::className<TypeTag>());

    propTagName.replace(0, strlen("PTag("), "");
    int n = propTagName.length();
    propTagName.replace(n - 1, 1, "");
    //TypeTagName.replace(0, strlen("Dumux::Properties::TTag::"), "");

    if (!getDiagnostic_(TypeTagName, propTagName, result, "")) {
        // check whether the property is a default property
        const PropertyRegistry::KeyList &keys =
            PropertyRegistry::getKeys(Dune::className<TTAG(__Default)>());
        PropertyRegistry::KeyList::const_iterator it = keys.begin();
        for (; it != keys.end(); ++it) {
            const PropertyRegistryKey &key = it->second;
            if (key.propertyName() != propTagName)
                continue; // property already printed

            std::ostringstream oss;
            oss << "fallback " << key.propertyName()
                << " defined at " << key.fileDefined()
                << ":" << key.lineDefined()
                <<"\n";
            result = oss.str();
        }
    }


    return result;
}

inline void print_(const std::string &typeTagName,
                   std::ostream &os,
                   const std::string &indent,
                   std::set<std::string> &printedProperties)
{
    if (indent == "")
        os << indent << "Properties for " << canonicalTypeTagNameToName_(typeTagName) << ":";
    else
        os << indent << "Inherited from " << canonicalTypeTagNameToName_(typeTagName) << ":";
    const PropertyRegistry::KeyList &keys =
        PropertyRegistry::getKeys(typeTagName);
    PropertyRegistry::KeyList::const_iterator it = keys.begin();
    bool somethingPrinted = false;
    for (; it != keys.end(); ++it) {
        const PropertyRegistryKey &key = it->second;
        if (printedProperties.count(key.propertyName()) > 0)
            continue; // property already printed
        if (!somethingPrinted) {
            os << "\n";
            somethingPrinted = true;
        }
        os << indent << "  "
           << key.propertyKind() << " " << key.propertyName();
        if (key.propertyKind() != "opaque")
            os << " = '" << key.propertyValue() << "'";
        os << " defined at " << key.fileDefined()
           << ":" << key.lineDefined()
           << "\n";
        printedProperties.insert(key.propertyName());
    };
    if (!somethingPrinted)
        os << " (none)\n";
    // print properties defined on children
    using ChildrenList = TypeTagRegistry::ChildrenList;
    const ChildrenList &children = TypeTagRegistry::children(typeTagName);
    ChildrenList::const_iterator ttagIt = children.begin();
    std::string newIndent = indent + "  ";
    for (; ttagIt != children.end(); ++ttagIt) {
        print_(*ttagIt, os, newIndent, printedProperties);
    }
}

//! \internal
template <class TypeTag>
void print(std::ostream &os = std::cout)
{
    std::set<std::string> printedProps;
    print_(Dune::className<TypeTag>(), os, "", printedProps);

    // print the default properties
    const PropertyRegistry::KeyList &keys =
        PropertyRegistry::getKeys(Dune::className<TTAG(__Default)>());
    PropertyRegistry::KeyList::const_iterator it = keys.begin();
    for (; it != keys.end(); ++it) {
        const PropertyRegistryKey &key = it->second;
        if (printedProps.count(key.propertyName()) > 0)
            continue; // property already printed
        os << "  default " << key.propertyName()
           << " (" << key.fileDefined()
           << ":" << key.lineDefined()
           << ")\n";
        printedProps.insert(key.propertyName());
    };
}
#else // !defined NO_PROPERTY_INTROSPECTION
template <class TypeTag>
void print(std::ostream &os = std::cout)
{
    std::cout <<
        "The Dumux property system was compiled with the macro\n"
        "NO_PROPERTY_INTROSPECTION defined.\n"
        "No diagnostic messages this time, sorry.\n";
}

template <class TypeTag>
const std::string getDiagnostic(std::string propTagName)
{
    std::string result;
    result =
        "The Dumux property system was compiled with the macro\n"
        "NO_PROPERTY_INTROSPECTION defined.\n"
        "No diagnostic messages this time, sorry.\n";
    return result;
};

#endif // !defined NO_PROPERTY_INTROSPECTION


template <class TypeTag, class Ancestor>
class AncestorAdder;

template <class TypeTag>
class TypeTagAncestors
{
public:
    using AncestorMatrix = std::vector<std::vector<std::string> >;

    static void addAncestors()
    {
        AncestorAdder<TypeTag, TypeTag>::add();
    }

    static AncestorMatrix& ancestors()
    {
        static AncestorMatrix ancestors_;
        return ancestors_;
    };

    static int& row()
    {
        static int row_ = 0;
        return row_;
    }

    static void print(std::ostream& os)
    {
        const auto& a = ancestors();

        std::vector<size_t> colSizes;
        for (size_t i = 0; i < a.size(); ++i)
        {
            auto numCols = a[i].size();
            colSizes.resize(numCols);

            std::vector<bool> bifurcation(numCols, false);
            for (size_t j = 0; j < numCols; ++j)
            {
                for (int below = i+1; below < a.size(); ++below)
                {
                    if (j > 0 && a[below].size() > j && !a[below][j].empty())
                    {
                        if (a[below][j-1].empty())
                        {
                            bifurcation[j] = true;
                            break;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                if (a[i][j].empty())
                {
                    if (bifurcation[j])
                        os << " | ";
                    else if (j > 0)
                        os << "   ";
                    os << std::string(colSizes[j], ' ');
                }
                else
                {
                    if (j > 0)
                        os << " - ";
                    os << a[i][j];
                    colSizes[j] = a[i][j].size();
                }
            }
            os << std::endl;

            for (size_t j = 0; j < numCols; ++j)
            {
                if (bifurcation[j])
                    os << " | ";
                else if (j > 0)
                    os << "   ";
                os << std::string(colSizes[j], ' ');
            }
            os << std::endl;
        }
    }

    static bool contains(const std::string& str)
    {
        const auto& a = ancestors();

        for (size_t i = 0; i < a.size(); ++i)
        {
            for (size_t j = 0; j < a[i].size(); ++j)
            {
                auto found = a[i][j].find(str);
                if (found != std::string::npos)
                    return true;
            }
        }

        return false;
    }
};

template <class TypeTag, class Ancestor>
class AncestorAdder
{
    using TTAncestors = TypeTagAncestors<TypeTag>;

public:
    static void add(unsigned col = 0)
    {
        auto& ancestors = TTAncestors::ancestors();

        auto ancestorFull = Dune::className<Ancestor>();
        // 25 is the size of "Dumux::Properties::TTag"
        auto ancestorName = ancestorFull.substr(25, ancestorFull.size() - 25);

        if (ancestors.size() < TTAncestors::row()+1)
        {
            ancestors.resize(TTAncestors::row()+1);
        }

        if (ancestors[TTAncestors::row()].size() < col+1)
        {
            ancestors[TTAncestors::row()].resize(col+1);
        }
        ancestors[TTAncestors::row()][col] = ancestorName;

        AncestorAdder<TypeTag, typename Ancestor::Child1>::add(++col);

        if (Dune::className<typename Ancestor::Child2>() != "void")
        {
            TTAncestors::row()++;
        }
        AncestorAdder<TypeTag, typename Ancestor::Child2>::add(col);

        if (Dune::className<typename Ancestor::Child3>() != "void")
        {
            TTAncestors::row()++;
        }
        AncestorAdder<TypeTag, typename Ancestor::Child3>::add(col);

        if (Dune::className<typename Ancestor::Child4>() != "void")
        {
            TTAncestors::row()++;
        }
        AncestorAdder<TypeTag, typename Ancestor::Child4>::add(col);

        if (Dune::className<typename Ancestor::Child5>() != "void")
        {
            TTAncestors::row()++;
        }
        AncestorAdder<TypeTag, typename Ancestor::Child5>::add(col);
    }
};

template <class TypeTag>
class AncestorAdder<TypeTag, void>
{
public:
    static void add(unsigned col = 0){}
};

//! \endcond

} // namespace Properties
} // namespace Dumux

#endif // DUMUX_PROPERTY_SYSTEM_HH
