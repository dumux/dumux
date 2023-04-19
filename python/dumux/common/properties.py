# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The DuMux property system in Python consisting of
Property, TypeTag and Model
"""

import os
import re
from dataclasses import dataclass
from typing import List, Union
from dune.common.hashit import hashIt
import dumux


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
    propertyString += f"struct {propertyName}<TypeTag, TTag::{typeTagName}>\n{{\n"

    if isinstance(typeArg, (float)):
        propertyString += "    using type = double;\n"
    elif isinstance(typeArg, (int)):
        propertyString += "    using type = int;\n"
    elif isinstance(typeArg, Property):
        if typeArg.requiredPropertyTypes or typeArg.requiredPropertyValues:
            propertyString += "private:\n"
        if typeArg.requiredPropertyTypes is not None:
            for reqProp in typeArg.requiredPropertyTypes:
                propertyString += (
                    f"    using {reqProp} = " f"GetPropType<TypeTag, Properties::{reqProp}>;\n"
                )
        if typeArg.requiredPropertyValues is not None:
            for reqProp in typeArg.requiredPropertyValues:
                reqPropLowerCase = reqProp[0].lower() + reqProp[1:]
                propertyString += (
                    f"    static constexpr auto {reqPropLowerCase} = "
                    f"getPropValue<TypeTag, Properties::{reqProp}>();\n"
                )

        propertyString += "public:\n"
        propertyString += f"    using type = {typeArg.cppType};\n"

    propertyString += "};"

    return propertyString


def valuePropertyToString(propertyName, typeTagName, value):
    """Converts a Value Property to a string"""

    propertyString = "template<class TypeTag>\n"
    propertyString += f"struct {propertyName}<TypeTag, TTag::{typeTagName}>\n{{"

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
    print("The following TypeTags are available:")
    for key, value in _typeTags.items():
        print(key, ":", value["description"])
    print("\n**********************************")


def propertiesHeaderPath():
    """Find the path to the properties.hh C++ header"""

    path, _ = os.path.split(dumux.__file__)
    metaDataFile = os.path.join(path, "data/dumux.cmake")
    if os.path.exists(metaDataFile):
        data = {}
        with open(metaDataFile, "r") as metaData:
            for line in metaData:
                try:
                    key, value = line.split("=", 1)
                    data[key] = value.strip()
                except ValueError:  # no '=' in line
                    pass
        return os.path.abspath(
            os.path.join(
                data["DEPBUILDDIRS"].split(";")[0],
                "python",
                "properties.hh",
            )
        )

    # as fall-back try relative path
    propertiesHeader = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "../../../../dumux/common/properties.hh",
        )
    )

    if os.path.exists(propertiesHeader):
        return propertiesHeader

    raise RuntimeError("Could not find properties.hh header")


def predefinedProperties():
    """Create a list of properties defined in properties.hh"""

    propertiesHeader = propertiesHeaderPath()
    with open(propertiesHeader, encoding="utf-8") as header:
        properties = []
        pattern = re.compile(r"DUMUX_DEFINE_PROPERTY\((.*?)\)")
        for line in header:
            if line.startswith("DUMUX_DEFINE_PROPERTY("):
                properties.append(pattern.search(line).group(1))
            elif line.startswith("struct"):
                properties.append(line.split(" ")[1])
        return properties


def predefinedTypeTags(inheritsFrom=None):
    """Create a list of predefined TypeTags"""
    result = _typeTags
    if inheritsFrom is not None:
        for parentTypeTag in inheritsFrom:
            if isinstance(parentTypeTag, CppTypeTag):
                result[parentTypeTag.name] = {
                    "include": parentTypeTag.includes[0],
                    "description": parentTypeTag.description,
                }
    return result


class CppTypeTag:
    """Creates a TypeTag from an existing c++ TypeTag"""

    def __init__(self, *, name, include, description=None):
        self.name = name
        self.includes = [include]
        self.description = description
        self.isExistingTypeTag = True
        self.properties = {}
        self.newPropertyDefinitions = []


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

        existingTypeTags = predefinedTypeTags(self.inheritsFrom)

        if self.name in existingTypeTags:
            if inheritsFrom is not None:
                raise ValueError(
                    f"Existing TypeTag {name} cannot inherit from other TypeTags."
                    f" Use TypeTag({name}) only."
                )
            self.isExistingTypeTag = True
            self.includes = [existingTypeTags[self.name]["include"]]
        else:
            self.isExistingTypeTag = False

        if self.inheritsFrom is not None:
            # treat existing TypeTags by converting the given string to a real TypeTag instance
            for idx, parentTypeTag in enumerate(self.inheritsFrom):
                if not isinstance(parentTypeTag, (CppTypeTag, TypeTag)):
                    if not isinstance(parentTypeTag, str):
                        raise ValueError(
                            f"Unknown parent TypeTag {parentTypeTag}. "
                            "Use either argument of type TypeTag "
                            "or a string for an existing TypeTag. "
                            f"List of existing TypeTags: {existingTypeTags.keys()}"
                        )
                    if parentTypeTag not in existingTypeTags:
                        raise ValueError(
                            f"Unknown TypeTag {parentTypeTag}. "
                            f"List of existing TypeTags: {existingTypeTags.keys()}"
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

        file = f"#ifndef DUMUX_{self.name.upper()}_PROPERTIES_HH\n"
        file += f"#define DUMUX_{self.name.upper()}_PROPERTIES_HH\n\n"

        file += "#include <dumux/common/properties.hh>\n"

        for include in self.includes:
            assert "<" not in include
            file += f"#include <{include}>\n"

        file += "\n"

        file += "namespace Dumux::Properties {\n\n"
        file += "namespace TTag {\n"

        if self.inheritsFrom is not None:
            for otherTypeTag in self.inheritsFrom:
                if not otherTypeTag.isExistingTypeTag:
                    file += f"struct {otherTypeTag.name} {{}}; \n"

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
            file += typePropertyToString(
                "Grid", self.name, Property.fromCppType("typename TypeTag::GridGeometry::Grid")
            )
            file += "\n\n"

        file += "} // end namespace Dumux::Properties \n\n"
        file += "#endif"

        return file


class Model(TypeTag):
    """A DuMux model specifies all properties necessary to build a DuMux simulator"""

    def __init__(self, *, inheritsFrom, gridGeometry, scalar: str = "double"):

        inheritsFromNames = []
        for parentModel in inheritsFrom:
            if isinstance(parentModel, str):
                inheritsFromNames.append(parentModel)
            elif isinstance(parentModel, (CppTypeTag, TypeTag)):
                inheritsFromNames.append(parentModel.name)
            else:
                raise ValueError(
                    "Arguments passed to inheritsFrom must be either names (strings) of "
                    "known TypeTags, TypeTag objects or CppTypeTag objects."
                )

        # generate a generic name
        genericName = "TypeTag" + hashIt(
            "".join(inheritsFromNames) + gridGeometry._typeName + scalar
        )

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
