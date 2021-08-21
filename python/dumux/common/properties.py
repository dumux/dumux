"""
The DuMux property system in Python consisting of
Property, TypeTag and Model
"""

import os
from dataclasses import dataclass
from typing import List, Union
from dune.common.hashit import hashIt


@dataclass
class Property:
    """
    Properties are used to construct a model
    Instances of Property are created with
    Property.fromInstance(...), Property.fromCppType(...),
    or  Property.fromValue(...).
    """

    cppType: str = None
    cppIncludes: List[str] = None
    requiredPropertyTypes: List[str] = None
    requiredPropertyValues: List[str] = None
    value: Union[bool, int, float] = None

    @classmethod
    def fromInstance(cls, inst):
        """Create a Property from an instance of a wrapper"""
        if not hasattr(inst, "_typeName"):
            raise TypeError(
                "The given instance {inst} does not have an attribute _typeName. "
                "Only generated C++ objects (with Python bindings) are accepted."
            )

        return cls(cppType=inst._typeName, cppIncludes=inst._includes)

    @classmethod
    def fromCppType(
        cls,
        cppType: str,
        *,
        cppIncludes: List[str] = None,
        requiredPropertyTypes: List[str] = None,
        requiredPropertyValues: List[str] = None,
    ):
        """Create a Property from a given C++ type and includes"""
        return cls(
            cppType=cppType,
            cppIncludes=cppIncludes,
            requiredPropertyTypes=requiredPropertyTypes,
            requiredPropertyValues=requiredPropertyValues,
        )

    @classmethod
    def fromValue(cls, value):
        """Create a Property from a given value"""
        return cls(value=value)


def typePropertyToString(propertyName, typeTagName, typeArg):
    """Converts a Type Property to a string"""

    propertyString = "template<class TypeTag>\n"
    propertyString += "struct {}<TypeTag, TTag::{}>\n{{\n".format(propertyName, typeTagName)

    if isinstance(typeArg, (float)):
        propertyString += "    using type = {};\n".format("double")
    elif isinstance(typeArg, (int)):
        propertyString += "    using type = {};\n".format("int")
    elif isinstance(typeArg, Property):
        if typeArg.requiredPropertyTypes or typeArg.requiredPropertyValues:
            propertyString += "private:\n"
        if typeArg.requiredPropertyTypes is not None:
            for reqProp in typeArg.requiredPropertyTypes:
                propertyString += "    using {} = {};\n".format(
                    reqProp, "GetPropType<TypeTag, Properties::{}>".format(reqProp)
                )
        if typeArg.requiredPropertyValues is not None:
            for reqProp in typeArg.requiredPropertyValues:
                reqPropLowerCase = reqProp[0].lower() + reqProp[1:]
                propertyString += "    static constexpr auto {} = {};\n".format(
                    reqPropLowerCase, "getPropValue<TypeTag, Properties::{}>()".format(reqProp)
                )

        propertyString += "public:\n"
        propertyString += "    using type = {};\n".format(typeArg.cppType)

    propertyString += "};"

    return propertyString


def valuePropertyToString(propertyName, typeTagName, value):
    """Converts a Value Property to a string"""

    propertyString = "template<class TypeTag>\n"
    propertyString += "struct {}<TypeTag, TTag::{}>\n{{".format(propertyName, typeTagName)

    # make sure to get the correct C++ types and values
    if isinstance(value, bool):
        value = str(value).lower()
        cppType = "bool"
    elif isinstance(value, int):
        cppType = "int"
    elif isinstance(value, float):
        cppType = "Scalar"
        propertyString += "\nprivate:\n"
        propertyString += "    using Scalar = GetPropType<TypeTag, Properties::Scalar>;\n"
        propertyString += "public:"
    else:
        raise ValueError(f"Invalid argument {value}. Expects bool, int or float.")

    propertyString += "\n"
    propertyString += f"    static constexpr {cppType} value = {value};"
    propertyString += "\n"
    propertyString += "};"

    return propertyString


_typeTags = {
    "CCTpfaModel": {
        "include": "dumux/discretization/cctpfa.hh",
        "description": "A cell-centered two-point flux finite volume discretization scheme.",
    },
    "BoxModel": {
        "include": "dumux/discretization/box.hh",
        "description": "A node-centered two-point flux finite volume discretization scheme",
    },
    "OneP": {
        "include": "dumux/porousmediumflow/1p/model.hh",
        "description": "A model for single-phase flow in porous media.",
    },
}


def listTypeTags():
    """List all available TypeTags/Models that can be inherited from"""

    print("\n**********************************\n")
    print("The following TypeTags are availabe:")
    for key, value in _typeTags.items():
        print(key, ":", value["description"])
    print("\n**********************************")


def predefinedProperties():
    """Create a list of properties defined in properties.hh"""

    propertiesHeader = os.path.abspath(
        os.path.dirname(__file__) + "/../../../../dumux/common/properties.hh"
    )
    with open(propertiesHeader, encoding="utf-8") as header:
        properties = []
        for line in header:
            if line.startswith("struct"):
                properties.append(line.split(" ")[1])
        return properties


class TypeTag:
    """TypeTags are inheritable collections of properties"""

    knownProperties = predefinedProperties()

    def __init__(self, name, *, inheritsFrom=None, gridGeometry=None):
        self.inheritsFrom = inheritsFrom
        self.includes = []
        self.properties = {}
        self.newPropertyDefinitions = []
        self.gridGeometry = gridGeometry
        self.name = name

        if self.name in _typeTags.keys():
            if inheritsFrom is not None:
                raise ValueError(
                    f"Existing TypeTag {name} cannot inherit from other TypeTags."
                    f" Use TypeTag({name}) only."
                )
            self.isExistingTypeTag = True
            self.includes = [_typeTags[self.name]["include"]]
        else:
            self.isExistingTypeTag = False

        if self.inheritsFrom is not None:
            # treat existing TypeTags by converting the given string to a real TypeTag instance
            for idx, parentTypeTag in enumerate(self.inheritsFrom):
                if not isinstance(parentTypeTag, TypeTag):
                    if not isinstance(parentTypeTag, str):
                        raise ValueError(
                            "Unknown parent TypeTag {}. Use either argument of type TypeTag "
                            "or a string for an existing TypeTag. "
                            "List of existing TypeTags: {}".format(parentTypeTag, _typeTags.keys())
                        )
                    if parentTypeTag not in _typeTags.keys():
                        raise ValueError(
                            "Unknown TypeTag {}. List of existing TypeTags: {}".format(
                                parentTypeTag, _typeTags.keys()
                            )
                        )
                    self.inheritsFrom[idx] = TypeTag(parentTypeTag)

            # pick up the properties and includes of the parent TypeTag
            for parentTypeTag in reversed(self.inheritsFrom):
                self.newPropertyDefinitions += parentTypeTag.newPropertyDefinitions
                for key in parentTypeTag.properties:
                    self.properties[key] = parentTypeTag.properties[key]
                if parentTypeTag.includes is not None:
                    for include in parentTypeTag.includes:
                        self.includes.append(include)

    def __setitem__(self, key, value: Property):
        """the [] operator for setting values"""
        if not isinstance(value, Property):
            raise ValueError("Only values of type Property can be assigned to a model")

        if key not in self.knownProperties:
            print(f"Adding {key} as new property")
            self.newPropertyDefinitions += [key]

        self.properties[key] = value
        if value.cppIncludes is not None:
            for include in value.cppIncludes:
                self.includes.append(include)

    def __getitem__(self, key):
        """the [] operator for getting values"""
        return self.properties[key]

    @property
    def cppType(self):
        """Returns the TypeTag as a string"""
        return "Dumux::Properties::TTag::" + self.name

    @property
    def cppHeader(self):
        """creates a string resembling a properties.hh file"""

        file = "#ifndef DUMUX_{}_PROPERTIES_HH\n".format(self.name.upper())
        file += "#define DUMUX_{}_PROPERTIES_HH\n\n".format(self.name.upper())

        file += "#include <dumux/common/properties.hh>\n"

        for include in self.includes:
            assert "<" not in include
            file += "#include <{}>\n".format(include)

        file += "\n"

        file += "namespace Dumux::Properties {\n\n"
        file += "namespace TTag {\n"

        if self.inheritsFrom is not None:
            for otherTypeTag in self.inheritsFrom:
                if not otherTypeTag.isExistingTypeTag:
                    file += "struct {}{{}}; \n".format(otherTypeTag.name)

            args = ",".join(x.name for x in self.inheritsFrom)
        else:
            args = ""

        file += f"struct {self.name}\n{{\n"
        file += f"    using InheritsFrom = std::tuple<{args}>;\n"
        if self.gridGeometry is not None:
            file += f"    using GridGeometry = {self.gridGeometry._typeName};\n"
        file += "};\n} // end namespace TTag\n\n"

        for newDef in self.newPropertyDefinitions:
            file += "template<class TypeTag, class MyTypeTag>\n"
            file += "struct " + newDef + " { using type =  UndefinedProperty; };\n\n"

        for prop in self.properties:
            if self[prop].value is not None:
                file += valuePropertyToString(prop, self.name, self[prop].value) + "\n\n"
            else:
                file += typePropertyToString(prop, self.name, self[prop]) + "\n\n"

        if self.gridGeometry is not None:
            file += (
                typePropertyToString(
                    "Grid", self.name, Property.fromCppType("typename TypeTag::GridGeometry::Grid")
                )
                + "\n\n"
            )

        file += "} // end namespace Dumux::Properties \n\n"
        file += "#endif"

        return file


class Model(TypeTag):
    """A DuMux model specifies all properties necessary to build a DuMux simulator"""

    def __init__(self, *, inheritsFrom: List[str], gridGeometry, scalar: str = "double"):
        # generate a generic name
        genericName = "TypeTag" + hashIt("".join(inheritsFrom) + gridGeometry._typeName + scalar)

        # deduce the discretization tag from the grid geometry
        discretizationMethod = gridGeometry.discMethod
        discretizationMap = {
            "box": "BoxModel",
            "cctpfa": "CCTpfaModel",
        }
        if discretizationMethod in discretizationMap:
            inheritsFrom += [discretizationMap[discretizationMethod]]

        # call parent constructor
        super().__init__(name=genericName, inheritsFrom=inheritsFrom, gridGeometry=gridGeometry)

        # set the scalar type
        self.__setitem__("Scalar", Property.fromCppType(scalar))
