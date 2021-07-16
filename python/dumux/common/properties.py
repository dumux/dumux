import os
import string
import random
from dune.common.hashit import hashIt

"""
Properties
"""


class Property:
    """"Properties are used to construct a model"""

    def __init__(self, object=None, value=None, type=None, includes=[], requiredProperties=[]):
        if object is not None:
            assert(hasattr(object, '_typeName'))
            if type is not None or value is not None:
                raise ValueError("The Property constructor expects exactly one of the following arguments: object, type, or value.")
            if includes or requiredProperties:
                raise ValueError("The arguments includes and requiredProperties are ignored if the object argument is specified.")
            self._typeName = object._typeName
            self._includes = object._includes if hasattr(object, '_includes') else []
            self._requiredPropertyTypes = object._requiredPropertyTypes if hasattr(object, '_requiredPropertyTypes') else []
        elif value is not None:
            if object is not None or type is not None:
                raise ValueError("The Property constructor expects exactly one of the following arguments: object, type, or value.")
            if includes or requiredProperties:
                raise ValueError("The arguments includes and requiredProperties are ignored if the value argument is specified.")
            self._value = value
        elif type is not None:
            if object is not None or value is not None:
                raise ValueError("The Property constructor expects exactly one of the following arguments: object, type, or value.")
            self._typeName = type
            self._includes = includes
            self._requiredPropertyTypes = requiredProperties
        else:
            raise ValueError("The Property constructor expects exactly one of the following arguments: object, type, or value.")


def typePropertyToString(propertyName, typeTagName, typeArg):
    """Converts a Type Property to a string"""

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


def valuePropertyToString(propertyName, typeTagName, value):
    """Converts a Value Property to a string"""

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


TYPETAGS = {
    'CCTpfaModel': {
        'include': 'dumux/discretization/cctpfa.hh',
        'description': 'A cell-centered two-point flux finite volume discretization scheme.'
    },
    'BoxModel': {
        'include': 'dumux/discretization/box.hh',
        'description': 'A node-centered two-point flux finite volume discretization scheme'
    },
    'OneP': {
        'include': 'dumux/porousmediumflow/1p/model.hh',
        'description': 'A model for single-phase flow in porous media.'
    },
}


def listTypeTags():
    """List all available TypeTags/Models that can be inherited from"""

    print("\n**********************************\n")
    print("The following TypeTags are availabe:")
    for key in TYPETAGS.keys():
        print(key, ":", TYPETAGS[key]['description'])
    print("\n**********************************")


def getKnownProperties():
    filepath = os.path.abspath(os.path.dirname(__file__) + '/../../../../dumux/common/properties.hh')
    with open(filepath) as f:
        result = []
        for line in f:
            if line.startswith('struct'):
                result.append(line.split(' ')[1])
        return result


class TypeTag:
    knownProperties = getKnownProperties()

    def __init__(self, name=None, *, inheritsFrom=None, gridGeometry=None, scalar='double'):
        self.inheritsFrom = inheritsFrom
        self.includes = []
        self.properties = {}
        self.newPropertyDefinitions = []
        self.gridGeometry = gridGeometry

        if name is not None:
            self.name = name
        else:
            if gridGeometry is None and inheritsFrom is None:
                self.name = "typetag_" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))
            else:
                self.name = "typetag_" + hashIt("".join(inheritsFrom) + gridGeometry._typeName + scalar)

        if self.name in TYPETAGS.keys():
            if inheritsFrom is not None:
                raise ValueError(f"Existing TypeTag {name} cannot inherit from other TypeTags. Use TypeTag({name}) only.")
            self.isExistingTypeTag = True
            self.includes = [TYPETAGS[self.name]['include']]
        else:
            self.isExistingTypeTag = False

        if self.gridGeometry is not None:
            discretizationMethod = self.gridGeometry.discMethod
            map = {
                "box": "BoxModel",
                "cctpfa": "CCTpfaModel",
            }
            if discretizationMethod in map:
                self.inheritsFrom += [map[discretizationMethod]]

        if self.inheritsFrom is not None:
            # treat existing TypeTags by converting the given string to a real TypeTag object
            for idx, parentTypeTag in enumerate(self.inheritsFrom):
                if not isinstance(parentTypeTag, TypeTag):
                    if not isinstance(parentTypeTag, str):
                        raise ValueError("Unknown parent TypeTag {}. Use either argument of type TypeTag "
                                         "or a string for an existing TypeTag. List of existing TypeTags: {}".format(parentTypeTag, TYPETAGS.keys()))
                    if parentTypeTag not in TYPETAGS.keys():
                        raise ValueError("Unknown TypeTag {}. List of existing TypeTags: {}".format(parentTypeTag, TYPETAGS.keys()))
                    self.inheritsFrom[idx] = TypeTag(parentTypeTag)

            # pick up the properties and includes of the parent TypeTag
            for parentTypeTag in reversed(self.inheritsFrom):
                self.newPropertyDefinitions += parentTypeTag.newPropertyDefinitions
                for key in parentTypeTag.properties:
                    self.properties[key] = parentTypeTag.properties[key]
                if parentTypeTag.includes is not None:
                    for include in parentTypeTag.includes:
                        self.includes.append(include)

        self._typeName = 'Dumux::Properties::TTag::' + self.name

        # set the scalar type
        self.__setitem__('Scalar', Property(type=scalar))

    # the [] operator for setting values
    def __setitem__(self, key, value):
        if not isinstance(value, Property):
            raise ValueError('Only values of type Property can be assigned to a model')

        if key not in self.knownProperties:
            print(f'Adding {key} as new property')
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

        file += f'struct {self.name}\n{{\n'
        file += f'    using InheritsFrom = std::tuple<{args}>;\n'
        if self.gridGeometry is not None:
           file += f'    using GridGeometry = {self.gridGeometry._typeName};\n'
        file += '};\n} // end namespace TTag\n\n'

        for newDef in self.newPropertyDefinitions:
            file += 'template<class TypeTag, class MyTypeTag>\n'
            file += 'struct ' + newDef + ' { using type =  UndefinedProperty; };\n\n'

        for prop in self.properties:
            if hasattr(self[prop], '_value'):
                file += valuePropertyToString(prop, self.name, self[prop]._value) +'\n\n'
            else:
                file += typePropertyToString(prop, self.name, self[prop]) +'\n\n'

        if self.gridGeometry is not None:
            file += typePropertyToString('Grid', self.name, Property(type='typename TypeTag::GridGeometry::Grid')) +'\n\n'

        file += '} // end namespace Dumux::Properties \n\n'
        file += '#endif'

        return file


class Model(TypeTag):
    # TODO maybe rename TypeTag to Model and remove this class here
    pass
