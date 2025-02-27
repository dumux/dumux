# Changing property name

Suppose that we want to change the name of the property `OldProperty` to `NewProperty`. In order to stay backwards compatible, a user should be able to use both the new name and the old name. If he uses the old name, he should get a deprecation warning, but his code should still do the same as before.

The only working possibility is to have two properties `OldProperty` and `NewProperty` and to _set by default the new one to the old one_. The other way around, a user could only invoke `OldProperty` in his using declarations if he doesn't specialize `NewProperty` in his own code.

## Actions to be done in the module `dumux` before release

#### Deprecate the old property

```c++
template<class TypeTag, class MyTypeTag>
struct [[deprecated("Use OldProperty instead.")]] OldProperty { using type = UndefinedProperty; };
```
Properties are usually defined in `dumux/common/properties.hh`.

#### Define the new property

```c++
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

template<class TypeTag, class T>
struct NewPropertyHelper
{ using type = GetPropType<TypeTag, Properties::OldProperty>; };

template<class TypeTag>
struct NewPropertyHelper<TypeTag, UndefinedProperty>
{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct NewProperty
{
    using type = typename NewPropertyHelper<TypeTag, typename OldProperty<TypeTag, MyTypeTag>::type>::type;
};

#pragma GCC diagnostic pop
```
The ignore pragmas avoid the emission of deprecation warnings for using the old property inside.
The indirection via the helper class is required to allow specializing the new property with an undefined/unspecialized old one.

#### Treat specializations

```c++
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

struct OldProperty<TypeTag, TTag::SpecialModel> {
    using type = SpecialType<...>;
};

#pragma GCC diagnostic pop
```
Specializations still have to be performed in terms of the old property, otherwise
user code can't employ the old name. Only the ignore pragmas have to be added.

#### Use the new name

Change all
```c++
using MyOldFavoriteName = GetPropType<TypeTag, Properties::OldProperty>;
```
to
```c++
using MyNewFavoriteName = GetPropType<TypeTag, Properties::NewProperty>;
```

#### Export the new name

Very often, `MyOldFavoriteName` in the snippet above coincides with `OldProperty`. In this case, it is advisable
to set `MyNewFavoriteName` to `NewProperty`. While this is unproblematic if the alias is local,
take special care for exported names:
```c++
public:
    using OldProperty = GetPropType<TypeTag, Properties::OldProperty>;
```
should become
```c++
public:
    using NewProperty = GetPropType<TypeTag, Properties::NewProperty>;
    using OldProperty [[deprecated("Use NewProperty instead.")]] = NewProperty;
```

## Actions to be done in derived modules and user code

#### Change using declarations

Change all
```c++
using MyOldFavoriteName = GetPropType<TypeTag, Properties::OldProperty>;
```
to
```c++
using MyNewFavoriteName = GetPropType<TypeTag, Properties::NewProperty>;
```
If `MyOldFavoriteName` is exported and you care about backwards compatibility, consider
deprecating it as outlined above.

#### Treat specializations

Change all
```c++
struct OldProperty<TypeTag, TTag::SpecialModel> {
    using type = SpecialType<...>;
};
```
to
```c++
struct NewProperty<TypeTag, TTag::SpecialModel> {
    using type = SpecialType<...>;
};
```
This should work because everywhere we now get the new property, so it is possible to overwrite it without breaking the code.

## Actions to be done in the module `dumux` after release

#### Only define the new property

Replace the two snippets "Deprecate the old property" and "Define the new property" by
```c++
template<class TypeTag, class MyTypeTag>
struct NewProperty { using type = UndefinedProperty; };
```

#### Only specialize the new property

Replace the dumux code base snippet "Treat specializations" by
```c++
struct NewProperty<TypeTag, TTag::SpecialModel> {
    using type = SpecialType<...>;
};
```

#### Remove exports of the old name

Remove the line containing the deprecation from the snippet "Export the new name".

### Known limitation

Unfortunately, clang emits false positives when invoking the outlined strategy.
A workaround is to prevent clang from emitting deprecation warnings triggered
by using `GetPropType` and `getPropValue()`. To this end, the line
```c++
     using LastType = Property<TypeTag, LastTypeTag>;
```
in `GetDefined` in `dumux/common/properties/propertysystem.hh` has to be augmented to
```c++
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
     using LastType = Property<TypeTag, LastTypeTag>;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
```
The lines
```c++
     using FirstType = Property<TypeTag, FirstTypeTag>;
```
and
```c++
template<class TypeTag, template<class,class> class Property>
using GetPropType = typename Properties::Detail::GetPropImpl<TypeTag, Property>::type::type;

template<class TypeTag, template<class,class> class Property>
constexpr auto getPropValue() { return Properties::Detail::GetPropImpl<TypeTag, Property>::type::value; }
```
have to be treated in the same way.

This should be augmented by a general warning message that is put best before the deprecation of the old property in `dumux/common/properties.hh`:
```c++
#if defined(__clang__) && !defined(DONT_EMIT_CLANG_NEWPROPERTY_WARNING)
#warning "The property `OldProperty` is deprecated in favor of `NewProperty` \
and will be removed after release X.Y. \
If clang is used, no deprecation warnings are emitted. \
We recommend to use gcc for getting rid of the warnings. \
You can suppress this message by defining the preprocessor variable \
DONT_EMIT_CLANG_NEWPROPERTY_WARNING."
#endif
```
