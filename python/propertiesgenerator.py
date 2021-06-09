#!/usr/bin/env python3

"""
Properties generator
"""

def typePropertyToString(propertyName, typeTagName, typeArg):
    propertyString = 'template<class TypeTag>\n'
    propertyString += 'struct {}<TypeTag, TTag::{}>\n{{\n'.format(propertyName, typeTagName)

    if isinstance(typeArg, (float)):
        propertyString += '    using type = {};\n'.format('double')
    elif isinstance(typeArg, (int)):
        propertyString += '    using type = {};\n'.format('int')
    elif not isinstance(typeArg, (str)):
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
    if isinstance(value, bool):
        value = str(value).lower()
    propertyString = 'template<class TypeTag>\n'
    propertyString += 'struct {}<TypeTag, TTag::{}>\n{{'.format(propertyName, typeTagName)
    propertyString += '\n    using value = {};\n'.format(value)
    propertyString += '};'

    # print(propertyString)
    return propertyString

class TypeTag:
    def __init__(self, name, inheritsFrom = None):
        self.name = name
        self.inheritsFrom = inheritsFrom
        self.includes = []
        self.properties = {}

        if inheritsFrom is not None:
            for parentTypeTag in reversed(inheritsFrom): # TODO check order
                for key in parentTypeTag.properties:
                    self.properties[key] = parentTypeTag.properties[key]

            for include in parentTypeTag.includes:
                self.includes.append(include)

        self._typeName = 'Dumux::Properties::TTag::' + name

    def __setitem__(self, key, value):
        self.properties[key] = value
        if hasattr(value, '_includes'):
            for include in value._includes:
                self.includes.append(include)

    def __getitem__(self, key):
        return self.properties[key]

    def getTypeTag():
        return 'Dumux::Properties::TTag::' + self.name

    def getProperties(self):
        file = '#ifndef DUMUX_{}_PROPERTIES_HH\n'.format(self.name)
        file += '#define DUMUX_{}_PROPERTIES_HH\n\n'.format(self.name)

        for include in self.includes:
            file += '#include {}\n'.format(include)

        file += '\n'

        file += 'namespace Dumux::Properties {\n\n'
        file += 'namespace TTag {\n'

        args = ",".join(x.name for x in self.inheritsFrom)

        file += 'struct {} {{ using inheritsFrom = std::tuple<{}>; }} \n'.format(self.name, args)
        file += '} // end namespace TTag\n\n'

        for prop in self.properties:
            if isinstance(self[prop], (bool, int, float)):
                file += valuePropertyToString(prop, self.name, self[prop]) +'\n\n'
            else:
                file += typePropertyToString(prop, self.name, self[prop]) +'\n\n'

        file += '} // end namespace Dumux::Properties \n\n'
        file += '#endif'

        return file


### Test stuff TODO remove

class Scalar:
    def __init__(self):
        self._typeName = 'double'
        self._includes = ['<dumux/common/math.hh>']

class CppInt:
    def __init__(self):
        self._typeName = 'int'


base = TypeTag('Base')
base['Scalar'] = CppInt()
base['UseMoles'] = True

second = TypeTag('Second')
second['Scalar'] = Scalar()
second['NumPhases'] = 4

class MyProblem:
    def __init__(self):
        self._typeName = 'Dumux::MyProblem<TypeTag, Scalar, numPhases>'
        self._includes = ['<test/myproblem.hh>']
        self._requiredPropertyTypes = ['Scalar']
        self._requiredPropertyValues = ['NumPhases']

model = TypeTag('MyModel', [second, base])
model['Problem'] = MyProblem()


file = model.getProperties()
print(file)



# # construct a GridGeometry from a gridView
# # the grid geometry is JIT compiled
# def FVAssembler(model, ...):
#     includes = model._includes
#     typeName = "Dumux::FVAssembler<" + model._typeName + ">"
#
#
#     moduleName = "fvassembler_" + hashIt(typeName)
#     holderType = "std::shared_ptr<{}>".format(typeName)
#     generator = SimpleGenerator("FVAssembler", "Dumux::Python")
#     module = generator.load(includes, typeName, moduleName, options=[holderType])
#     return module.FVAssembler(...)
