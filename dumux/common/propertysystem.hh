// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

/*!
 * \file
 * \brief Provides the magic behind the DuMuX property system.
 *
 * \defgroup Properties Property System
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
#ifndef DUMUX_PROPERTIES_HH
#define DUMUX_PROPERTIES_HH

// For is_base_of
#include <boost/type_traits.hpp>

// Integral Constant Expressions

// string formating
#include <boost/format.hpp>


namespace Dumux
{
namespace Properties
{
#if !defined NO_PROPERTY_INTROSPECTION
//! Internal macro which is only required if the property introspection is enabled
#define PROP_INFO_(EffTypeTagName, PropTagName) \
    template <>                                                         \
    struct PropertyInfo< TTAG(EffTypeTagName), PTAG(PropTagName) >      \
    {                                                                   \
        static std::string propertyName() { return #PropTagName ; }     \
        static std::string fileDefined() { return __FILE__; }           \
        static int lineDefined() { return __LINE__; }                   \
    };
#else
//! Internal macro which is only required if the property introspection is enabled
#define PROP_INFO_(EffTypeTagName, PropTagName)
#endif

// some macros for simplification

/*!
 * \brief Makes a type out of a type tag name
 */
#define TTAG(TypeTagName) ::Dumux::Properties::TTag::TypeTagName

/*!
 * \brief Makes a type out of a property tag name
 */
#define PTAG(PropTagName) ::Dumux::Properties::PTag::PropTagName

/*!
 * \brief Define a new type tag.
 *
 * A type tag can inherit the properties defined on up to five parent
 * type tags. Examples:
 *
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
 */
#define NEW_TYPE_TAG(TypeTagName, ...) \
    namespace TTag {                                                \
    struct TypeTagName : public TypeTag<TypeTagName, ##__VA_ARGS__> \
    {                                                               \
        static const std::string name() { return #TypeTagName ; }   \
    };                                                              \
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
 * NEW_PROP_TAG(blubbPropTag);
 * NEW_PROP_TAG(blabbPropTag);
 */
#define NEW_PROP_TAG(PTagName) \
    namespace PTag {                                       \
    struct PTagName; } extern int semicolonHack_
/*                                                         \
    {                                                      \
        static const char *name() { return #PTagName ;};   \
    };                                                     \
    } extern int semicolonHack_
*/


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
 *    typedef typename GET_PROP(TypeTag, blabbProp) blabb;
 *
 *    static int calculateInternal_(int arg)
 *    { return arg * blabb::value; };
 * };
 */
#define SET_PROP(EffTypeTagName, PropTagName) \
    template <class TypeTag>                                    \
    struct Property<TypeTag, \
                    TTAG(EffTypeTagName), \
                    PTAG(PropTagName)>;                         \
    PROP_INFO_(EffTypeTagName, PropTagName) \
    template <class TypeTag>                                    \
    struct Property<TypeTag, \
                    TTAG(EffTypeTagName), \
                    PTAG(PropTagName) >

/*!
 * \brief Set the default for a property.
 *
 * SET_PROP_DEFAULT works exactly like SET_PROP, except that it does
 * not require an effective type tag. Defaults are used whenever a
 * property was not explicitly set or explicitly unset for a type tag.
 *
 * Example:
 *
 * // set a default for the blabbPropTag property tag
 * SET_PROP_DEFAULT(blabbPropTag)
 * {
 *    static const int value = 3;
 * };
 */
#define SET_PROP_DEFAULT(PropTagName) \
    template <class TypeTag>                                            \
    struct DefaultProperty<TypeTag, PTAG(PropTagName)>;                 \
    template <>                                                         \
    struct PropertyInfo< void, PTAG(PropTagName) >                      \
    {                                                                   \
        static std::string propertyName() { return #PropTagName ; }     \
        static std::string typeTagName() { return "*" ; }               \
        static std::string fileDefined() { return __FILE__; }           \
        static int lineDefined() { return __LINE__; }                   \
    };                                                                  \
    template <class TypeTag>                                            \
    struct DefaultProperty<TypeTag, PTAG(PropTagName) >

/*!
 * \brief Explicitly unset a property for a type tag.
 *
 * This means that the property will not be inherited from the type
 * tag's parents and that no default will be used.
 *
 * Example:
 *
 * // make the blabbPropTag property undefined for the BarTypeTag.
 * UNSET_PROP(BarTypeTag, blabbPropTag);
 */
#define UNSET_PROP(EffTypeTagName, PropTagName) \
    template <>                                                 \
    struct PropertyUnset<TTAG(EffTypeTagName), \
                         PTAG(PropTagName) >;                   \
    PROP_INFO_(EffTypeTagName, PropTagName) \
    template <>                                                 \
    struct PropertyUnset<TTAG(EffTypeTagName), \
                         PTAG(PropTagName) >                    \
        : public PropertyExplicitlyUnset \
        {}

/*!
 * \brief Set a property to a simple constant integer value.
 *
 * The constant can be accessed by the 'value' attribute.
 */
#define SET_INT_PROP(EffTypeTagName, PropTagName, Value) \
    SET_PROP(EffTypeTagName, PropTagName) \
    {                                                           \
        typedef int type;                                       \
        static const int value = Value;                         \
    }

/*!
 * \brief Set a property to a simple constant boolean value.
 *
 * The constant can be accessed by the 'value' attribute.
 */
#define SET_BOOL_PROP(EffTypeTagName, PropTagName, Value) \
    SET_PROP(EffTypeTagName, PropTagName) \
    {                                                      \
        typedef bool type;                                 \
        static const bool value = Value;                   \
    }

/*!
 * \brief Set a property which defines a type.
 *
 * The type can be accessed by the 'type' attribute.
 */
#define SET_TYPE_PROP(EffTypeTagName, PropTagName, Type) \
    SET_PROP(EffTypeTagName, PropTagName) \
    {                                                      \
        typedef Type type;                                 \
    }

/*!
 * \brief Set a property to a simple constant scalar value.
 *
 * The constant can be accessed by the 'value' attribute. In order to
 * use this macro, the property tag "Scalar" needs to be defined for
 * the real type tag.
 */
#define SET_SCALAR_PROP(EffTypeTagName, PropTagName, Value) \
    SET_PROP(EffTypeTagName, PropTagName) \
    {                                                                   \
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;   \
    public:                                                             \
        typedef Scalar type;                                            \
        static const Scalar value = Value;                              \
    }

/*!
 * \brief Get the property for a type tag.
 *
 * If you use GET_PROP within a template and want to refer to some
 * type (including the property itself), GET_PROP must be preceeded by
 * the 'typename' keyword.
 */
#define GET_PROP(TypeTag, PropTag) \
    ::Dumux::Properties::GetProperty< TypeTag, PropTag>::p

/*!
 * \brief Access the 'value' attribute of a property for a type tag.
 *
 * This is just for convenience and equivalent to GET_PROP(TypeTag,
 * PropTag) :: value.  If the property doesn't have an attribute named
 * 'value', this yields a compiler error.
 */
#define GET_PROP_VALUE(TypeTag, PropTag) \
    ::Dumux::Properties::GetProperty< TypeTag, PropTag>::p::value

/*!
 * \brief Access the 'type' attribute of a property for a type tag.
 *
 * This is just for convenience and equivalent to GET_PROP(TypeTag,
 * PropTag) :: type.  If the property doesn't have an attribute named
 * 'type', this yields a compiler error. Also, if you use this macro
 * within a template, it must be preceeded by the 'typename' keyword.
 */
#define GET_PROP_TYPE(TypeTag, PropTag) \
    ::Dumux::Properties::GetProperty< TypeTag, \
                                     PropTag>::p::type

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
 * int main()
 * {
 *    std::cout << PROP_DIAGNOSTIC(FooBarTypeTag, blabbPropTag) << "\n";
 * };
 */
#define PROP_DIAGNOSTIC(TypeTag, PropTag) \
    ::Dumux::Properties::propertyDiagnostic< TypeTag, \
                                            TypeTag, \
                                            PropTag>::message()
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
 * int main()
 * {
 *    std::cout << PROP_DIAGNOSTIC(FooBarTypeTag, blabbPropTag) << "\n";
 * };
 */
#define PROP_DIAGNOSTIC(TypeTag, PropTag) "Property introspection disabled by NO_PROPERTY_INTROSPECTION"
#endif


//////////////////////////////////////////////
// some serious template kung fu. Don't look at it too closely, it
// might damage your brain!
//////////////////////////////////////////////

//! \cond 0

namespace PTag {}
namespace TTag {}

using namespace boost::type_traits;
using boost::is_void;
using boost::is_base_of;

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
template <class EffectiveTypeTag, class PropertyTag>
struct PropertyInfo
{
    static std::string typeTagName()
    { return "Unknown"; }

    static std::string propertyName()
    { return "Unknown (PropertyInfo is not specialized for this property)"; }

    static std::string fileDefined()
    { return __FILE__; }

    static int lineDefined()
    { return __LINE__; }
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
    const static bool value = boost::true_type::value;
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
    static const bool value = boost::false_type::value;
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
    static const bool value = boost::false_type::value;
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
    typedef propertyDefined<RealTypeTag, Tree, PropertyTag> definedWhere;

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
    typedef SelfT SelfType;

    typedef Child1T Child1;
    typedef Child2T Child2;
    typedef Child3T Child3;
    typedef Child4T Child4;
    typedef Child5T Child5;
};

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
    typedef DefaultProperty<RealTypeTag, PropertyTag>  p;
};

// property defined on self
//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 0>
{
    typedef Property<RealTypeTag, TypeTag, PropertyTag>   p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 1>
{
    typedef typename GetProperty<typename TypeTag::Child1, PropertyTag, RealTypeTag>::p p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 2>
{
    typedef typename GetProperty<typename TypeTag::Child2, PropertyTag, RealTypeTag>::p p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 3>
{
    typedef typename GetProperty<typename TypeTag::Child3, PropertyTag, RealTypeTag>::p p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 4>
{
    typedef typename GetProperty<typename TypeTag::Child4, PropertyTag, RealTypeTag>::p p;
};

//! \internal
template <class TypeTag, class PropertyTag, class RealTypeTag>
struct GetProperty<TypeTag, PropertyTag, RealTypeTag, 5>
{
    typedef typename GetProperty<typename TypeTag::Child5, PropertyTag, RealTypeTag>::p p;
};

#if !defined NO_PROPERTY_INTROSPECTION
//! \internal
template <class RealTypeTag, class Tree, class PropertyTag>
struct propertyDiagnostic
{
    typedef typename Tree::SelfType EffTypeTag;

    static const int inheritedFrom = propertyTagIndex<RealTypeTag, EffTypeTag, PropertyTag>::value;

    static const bool explicitlyUnset = propertyExplicitlyUnsetOnTree<Tree, PropertyTag>::value;
    static const bool explicitlyUnsetOnSelf = propertyExplicitlyUnset<Tree, PropertyTag>::value;

    typedef propertyDefined<RealTypeTag, Tree, PropertyTag> definedWhere;

    typedef PropertyInfo<EffTypeTag, PropertyTag>           propInfo;
    typedef PropertyInfo<void, PropertyTag>                 propInfoDefault;

    static const std::string typeTagName()
    {
        return Tree::SelfType::name();
    };

    static const std::string propertyName()
    {
        switch (inheritedFrom) {
        case -1: return propInfoDefault::propertyName();
        case -1000:
        case 0: return propInfo::propertyName();
        case 1: return Child1Diagnostic::propertyName();
        case 2: return Child2Diagnostic::propertyName();
        case 3: return Child3Diagnostic::propertyName();
        case 4: return Child4Diagnostic::propertyName();
        case 5: return Child5Diagnostic::propertyName();
        }
    };

    typedef propertyDiagnostic<RealTypeTag, typename Tree::Child1, PropertyTag> Child1Diagnostic;
    typedef propertyDiagnostic<RealTypeTag, typename Tree::Child2, PropertyTag> Child2Diagnostic;
    typedef propertyDiagnostic<RealTypeTag, typename Tree::Child3, PropertyTag> Child3Diagnostic;
    typedef propertyDiagnostic<RealTypeTag, typename Tree::Child4, PropertyTag> Child4Diagnostic;
    typedef propertyDiagnostic<RealTypeTag, typename Tree::Child5, PropertyTag> Child5Diagnostic;


    static std::string message(const std::string &indent="", bool topLevel=true)
    {
        std::string result;
        if (topLevel) {
            result =
                (boost::format("%sProperty '%s' for type tag '%s'\n")
                 %indent
                 %propertyName()
                 %typeTagName()
                    ).str();
        }
        std::string newIndent = indent + "  ";

        if (explicitlyUnsetOnSelf) {
            result += (boost::format("%sexplicitly unset at %s:%d\n")
                       %newIndent
                       %propInfo::fileDefined()
                       %propInfo::lineDefined()).str();
            return result;
        }
        else if (explicitlyUnset) {
            result += explicitlyUnsetMsg(indent);
            return result;
        };

        switch (inheritedFrom) {
        case -1:
            result += newIndent;
            result += (boost::format("default from %s:%d\n")
                       %propInfoDefault::fileDefined()
                       %propInfoDefault::lineDefined()).str();
            break;

        case 0:
            result += newIndent;
            result += (boost::format("defined at %s:%d\n")
                       %propInfo::fileDefined()
                       %propInfo::lineDefined()).str();
            break;

        case 1:
            result += explicitlyUnsetOnChildMsg(indent, 5);
            result += explicitlyUnsetOnChildMsg(indent, 4);
            result += explicitlyUnsetOnChildMsg(indent, 3);
            result += explicitlyUnsetOnChildMsg(indent, 2);
            result += newIndent;
            result += (boost::format("inherited from '%s'\n"
                                     "%s") %
                       Child1Diagnostic::typeTagName() %
                       Child1Diagnostic::message(newIndent, false)).str();
            break;

        case 2:
            result += explicitlyUnsetOnChildMsg(indent, 5);
            result += explicitlyUnsetOnChildMsg(indent, 4);
            result += explicitlyUnsetOnChildMsg(indent, 3);
            result += newIndent;
            result += (boost::format("inherited from '%s'\n"
                                     "%s") %
                       Child2Diagnostic::typeTagName() %
                       Child2Diagnostic::message(newIndent, false)).str();
            break;

        case 3:
            result += explicitlyUnsetOnChildMsg(indent, 5);
            result += explicitlyUnsetOnChildMsg(indent, 4);
            result += newIndent;
            result += (boost::format("inherited from '%s'\n"
                                     "%s") %
                       Child3Diagnostic::typeTagName() %
                       Child3Diagnostic::message(newIndent, false)).str();
            break;

        case 4:
            result += explicitlyUnsetOnChildMsg(indent, 5);
            result += newIndent;
            result += (boost::format("inherited from '%s'\n"
                                     "%s") %
                       Child4Diagnostic::typeTagName() %
                       Child4Diagnostic::message(newIndent, false)).str();
            break;

        case 5:
            result += newIndent;
            result += (boost::format("inherited from '%s'\n"
                                     "%s") %
                       Child5Diagnostic::typeTagName() %
                       Child5Diagnostic::message(newIndent, false)).str();
            break;

        case -1000:
            result += boost::format("was not set anywhere\n").str();
        };

        return result;
    };

    static std::string explicitlyUnsetMsg(const std::string &indent)
    {
        if (explicitlyUnsetOnSelf) {
            return (boost::format("%sexplicitly unset for '%s' at %s:%d\n")
                    %indent
                    %typeTagName()
                    %propInfo::fileDefined()
                    %propInfo::lineDefined()).str();
        }
        else if (explicitlyUnset) {
            std::string result = (boost::format("%sexplicitly unset for all parents of '%s':\n")
                                  %indent
                                  %typeTagName()).str();
            result += explicitlyUnsetOnChildMsg(indent, 5);
            result += explicitlyUnsetOnChildMsg(indent, 4);
            result += explicitlyUnsetOnChildMsg(indent, 3);
            result += explicitlyUnsetOnChildMsg(indent, 2);
            result += explicitlyUnsetOnChildMsg(indent, 1);
            return result;
        };
        return "";
    }

    static std::string explicitlyUnsetOnChildMsg(const std::string &indent, int unsetChild)
    {
        switch (unsetChild) {
        case 5:
            return Child5Diagnostic::explicitlyUnsetMsg(indent + "  ");
        case 4:
            return Child4Diagnostic::explicitlyUnsetMsg(indent + "  ");
        case 3:
            return Child3Diagnostic::explicitlyUnsetMsg(indent + "  ");
        case 2:
            return Child2Diagnostic::explicitlyUnsetMsg(indent + "  ");
        case 1:
            return Child1Diagnostic::explicitlyUnsetMsg(indent + "  ");
        }
        return "";
    }
};

//! \internal
template <class RealTypeTag, class PropertyTag>
struct propertyDiagnostic<RealTypeTag, void, PropertyTag>
{
    typedef PropertyInfo<void, PropertyTag>           propInfo;

    static const std::string propertyName()
    { return propInfo::propertyName(); }

    static const std::string typeTagName()
    { return "* (Default Value)"; }

    static std::string message(const std::string &ident="", bool topLevel=true)
    {
        return "";
    };

    static std::string explicitlyUnsetMsg(const std::string &indent, int numUnsetChilds=5)
    {
        return "";
    };
};
#endif // !defined NO_PROPERTY_INTROSPECTION

//! \endcond

} // namespace Properties
} // namespace Dumux

#endif
