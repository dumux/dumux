#!/usr/bin/env python3

"""
Properties
"""

class Property:
    def __init__(self, **kwargs):
        if 'object' in kwargs:
            object = kwargs.get('object')
            assert(hasattr(object, '_typeName'))
            assert('type' not in kwargs and 'value' not in kwargs)
            includes = object._includes if hasattr(object, '_includes') else []
            requiredPropertyTypes = object._requiredPropertyTypes if hasattr(object, '_requiredPropertyTypes') else []
            self._typeName = object._typeName
            self._includes = includes
            self._requiredPropertyTypes = requiredPropertyTypes
        elif 'type' in kwargs:
            assert('object' not in kwargs and 'value' not in kwargs)
            self._typeName = kwargs.get('type')
            self._includes = kwargs.get('includes', [])
            self._requiredPropertyTypes = kwargs.get('requiredProperties', [])
        else:
            assert('value' in kwargs)
            assert('object' not in kwargs and 'type' not in kwargs)
            self._value =  kwargs.get('value')


# Converts a Type Property to a string
def typePropertyToString(propertyName, typeTagName, typeArg):
    propertyString = 'template<class TypeTag>\n'
    propertyString += 'struct {}<TypeTag, TTag::{}>\n{{\n'.format(propertyName, typeTagName)

    if isinstance(typeArg, (float)):
        propertyString += '    using type = {};\n'.format('double')
    elif isinstance(typeArg, (int)):
        propertyString += '    using type = {};\n'.format('int')
    elif isinstance(typeArg, Property) or not isinstance(typeArg, (str)):
        if hasattr(typeArg, '_requiredPropertyTypes') or hasattr(typeArg, '_requiredPropertyValues'):
            propertyString += 'private:\n'
        if hasattr(typeArg, '_requiredPropertyTypes'):
            for reqProp in typeArg._requiredPropertyTypes:
                propertyString += '    using {} = {};\n'.format(reqProp, 'GetPropType<TypeTag, Properties::{}>'.format(reqProp))

        if hasattr(typeArg, '_requiredPropertyValues'):
            for reqProp in typeArg._requiredPropertyValues:
                reqPropLowerCase = reqProp[0].lower() + reqProp[1:]
                propertyString += '    static constexpr auto {} = {};\n'.format(reqPropLowerCase, 'getPropValue<TypeTag, Properties::{}>()'.format(reqProp))

        propertyString += 'public:\n'
        propertyString += '    using type = {};\n'.format(typeArg._typeName)

    propertyString += '};'

    return propertyString

# Converts a Value Property to a string
def valuePropertyToString(propertyName, typeTagName, value):
    propertyString = 'template<class TypeTag>\n'
    propertyString += 'struct {}<TypeTag, TTag::{}>\n{{'.format(propertyName, typeTagName)

    # make sure to get the correct C++ types and values
    if isinstance(value, bool):
        value = str(value).lower()
        type = 'bool'
    elif isinstance(value, int):
        type = 'int'
    else:
        type = 'Scalar'
        propertyString += '\nprivate:\n'
        propertyString += '    using Scalar = GetPropType<TypeTag, Properties::Scalar>;\n'
        propertyString += 'public:'

    propertyString += '\n    static constexpr {} value = {};\n'.format(type, value)
    propertyString += '};'

    return propertyString

existingTypeTags = {'CCTpfaModel':'dumux/discretization/cctpfa.hh',
                    'BoxModel':'dumux/discretization/box.hh',
                    'OneP':'dumux/porousmediumflow/1p/model.hh'}

def listTypeTags():
    print("The following TypeTags are availabe:")
    print(existingTypeTags.keys())

def getKnownProperties():
    with open('../../dumux/common/properties.hh') as f:
        result = []
        for line in f:
            if line.startswith('struct'):
                result.append(line.split(' ')[1])
        return result

class TypeTag:
    knownProperties = getKnownProperties()

    def __init__(self, name, *, inheritsFrom=None):
        self.name = name
        self.inheritsFrom = inheritsFrom
        self.includes = []
        self.properties = {}
        self.newPropertyDefinitions = []

        if name in existingTypeTags.keys():
            if inheritsFrom is not None:
                raise ValueError("Existing TypeTag {} cannot inherit from other TypeTags. Use TypeTag({}) only.".format(name, name))
            self.isExistingTypeTag = True
            self.includes = [existingTypeTags[name]]
        else:
            self.isExistingTypeTag = False

        if self.inheritsFrom is not None:
            # treat existing TypeTags by converting the given string to a real TypeTag object
            for idx, parentTypeTag in enumerate(self.inheritsFrom):
                if not isinstance(parentTypeTag, TypeTag):
                    if not isinstance(parentTypeTag, str):
                        raise ValueError("Unknown parent TypeTag {}. Use either argument of type TypeTag "
                                         "or a string for an existing TypeTag. List of existing TypeTags: {}".format(parentTypeTag, existingTypeTags.keys()))
                    if parentTypeTag not in existingTypeTags.keys():
                        raise ValueError("Unknown TypeTag {}. List of existing TypeTags: {}".format(parentTypeTag, existingTypeTags.keys()))
                    self.inheritsFrom[idx] = TypeTag(parentTypeTag)

            # pick up the properties and includes of the parent TypeTag
            for parentTypeTag in reversed(self.inheritsFrom):
                self.newPropertyDefinitions += parentTypeTag.newPropertyDefinitions
                for key in parentTypeTag.properties:
                    self.properties[key] = parentTypeTag.properties[key]
                if parentTypeTag.includes is not None:
                    for include in parentTypeTag.includes:
                        self.includes.append(include)

        self._typeName = 'Dumux::Properties::TTag::' + name

    # the [] operator for setting values
    def __setitem__(self, key, value):
        assert(isinstance(value, Property))
        if key not in self.knownProperties:
            print('Adding', key, 'as new property')
            self.newPropertyDefinitions += [key]

        self.properties[key] = value
        if hasattr(value, '_includes'):
            for include in value._includes:
                self.includes.append(include)

    # the [] operator for getting values
    def __getitem__(self, key):
        return self.properties[key]

    # returns the TypeTag as a string
    def getTypeTag(self):
        return 'Dumux::Properties::TTag::' + self.name

    # creates a string resembling a properties.hh file
    def getProperties(self):
        file = '#ifndef DUMUX_{}_PROPERTIES_HH\n'.format(self.name.upper())
        file += '#define DUMUX_{}_PROPERTIES_HH\n\n'.format(self.name.upper())

        file += '#include <dumux/common/properties.hh>\n'

        for include in self.includes:
            assert('<' not in include)
            file += '#include <{}>\n'.format(include)

        file += '\n'

        file += 'namespace Dumux::Properties {\n\n'
        file += 'namespace TTag {\n'

        if self.inheritsFrom is not None:
            for otherTypeTag in self.inheritsFrom:
                if not otherTypeTag.isExistingTypeTag:
                    file += 'struct {}{{}}; \n'.format(otherTypeTag.name)

            args = ",".join(x.name for x in self.inheritsFrom)
        else:
            args = ""

        file += 'struct {} {{ using InheritsFrom = std::tuple<{}>; }}; \n'.format(self.name, args)
        file += '} // end namespace TTag\n\n'

        for newDef in self.newPropertyDefinitions:
            file += 'template<class TypeTag, class MyTypeTag>\n'
            file += 'struct ' + newDef + ' { using type =  UndefinedProperty; };\n\n'

        for prop in self.properties:
            if hasattr(self[prop], '_value'):
                file += valuePropertyToString(prop, self.name, self[prop]._value) +'\n\n'
            else:
                file += typePropertyToString(prop, self.name, self[prop]) +'\n\n'

        file += '} // end namespace Dumux::Properties \n\n'
        file += '#endif'

        return file
