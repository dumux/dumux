#!/usr/bin/env python3

"""
Automatically updates newparameters.json by searching all *.hh files
for usage of getParam or getParamFromGroup.
"""

import os
import argparse
import json


# find the content of the given string between the first matching pair of opening/closing keys
def getEnclosedContent(string, openKey, closeKey):

    # cut off everything before the first occurence of openKey
    string = openKey + string.partition(openKey)[2]

    # get content between mathing pair
    rest = string.partition(closeKey)
    result, rest = rest[0] + closeKey, rest[2]
    while result.count(openKey) != result.count(closeKey):
        rest = rest.partition(closeKey)
        if rest[1] == "":
            raise IOError(
                'Could not get content between "{}" and "{}" in given string "{}"'.format(
                    openKey, closeKey, string
                )
            )
        result, rest = result + rest[0] + closeKey, rest[2]

    return result.partition(openKey)[2].rpartition(closeKey)[0]


# extract a parameter from a given line
def extractParamName(line):

    # split occurrences of getParam<T>(CALLARGS) or getParamFromGroup<T>(CALLARGS)
    # into the template arguments T and the function arguments CALLARGS
    if "getParamFromGroup<" in line:
        line = line.split("getParamFromGroup")[1]
        hasGroup = True
    elif "getParam<" in line:
        line = line.split("getParam")[1]
        hasGroup = False
    else:
        return {}

    # TODO: Support this also
    if line.count("getParam") > 1:
        raise IOError('Cannot process multiple occurrences of "getParam" in one line')

    # remove trailing spaces and cut off everything behind semicolon
    line = line.strip("\n").strip(" ").split(";")[0]

    # extract template arg between '<' and '>'
    paramType = getEnclosedContent(line, "<", ">")

    # extract function arguments
    functionArgs = line.partition("<" + paramType + ">")[2]
    functionArgs = getEnclosedContent(functionArgs, "(", ")")

    if hasGroup:
        functionArgs = functionArgs.partition(",")[2]
    functionArgs = functionArgs.partition(",")
    paramName = functionArgs[0]
    defaultValue = None if not functionArgs[2] else functionArgs[2]

    paramType = paramType.strip(" ")
    paramName = paramName.strip(" ")
    if defaultValue:
        defaultValue = defaultValue.strip(" ")

    # if interior spaces occur in the parameter name, we can't identify it
    if paramName[0] != '"' or paramName[-1] != '"' or " " in paramName:
        raise IOError("Could not correctly process parameter name")

    return {
        "paramType": paramType,
        "paramName": paramName.strip('"'),
        "defaultValue": defaultValue,
    }


# extract all parameters from a given file
def getParamsFromFile(file, errorMsg):
    parameters = []
    errors = {}
    with open(file) as f:
        for lineIdx, line in enumerate(f):
            try:
                param = extractParamName(line)
                if param:
                    parameters.append(param)
            except IOError as e:
                errors[lineIdx] = {"line": line.strip(), "message": e}

    # print encountered errors
    if errors:
        errorMsg.append(
            "\n\n{} parameter{} in file {} could not be retrieved automatically. Please check them yourself:".format(
                len(errors), "s" if len(errors) > 1 else "", file
            )
        )
        for lineIdx in errors:
            errorMsg.append("\n\t-> line {}: {}".format(lineIdx, errors[lineIdx]["line"]))
            errorMsg.append("\n\t\t-> error message: {}".format(errors[lineIdx]["message"]))
        errorMsg.append(" ")

    return parameters

class CheckExistAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        if os.path.isfile(values):
            setattr(namespace, self.dest, values)
            setattr(namespace, "hasInput", True)
        else:
            parser.error("File {} does not exist!".format(values))


parser = argparse.ArgumentParser(
    description="""
    Generate a json file of parameters list from header files.
    ----------------------------------------------------------------
    If input file is given, the descriptions will be copied from the input.
    If no input file is given, the descriptions will be empty.
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "--root",
    help="root path of Dumux",
    metavar="rootpath",
    default=os.path.abspath(os.path.join(os.path.abspath(__file__), "../../../")),
)
parser.add_argument(
    "--input",
    help="json file of given paratemers",
    action=CheckExistAction,
    metavar="input",
    dest="inputFile",
)
"""
parser.add_argument(
    "--output",
    help="relative path (to the root path) of the output file",
    metavar="output",
    default="doc/doxygen/extradoc/newparameterlist.json",
)"""

args = vars(parser.parse_args())

maxExplanationWidth = 0

# get the explanations from current parameter json file
oldDataDict = {}
if args["inputFile"]:
    with open(args["inputFile"], "r") as f:
        data = f.read()
        oldDataDict = json.loads(data)

# search all *.hh files for parameters
# TODO: allow runtime args with extensions and folder(s) to be checked
parameters = []
errorLog = []
for root, _, files in os.walk(args["root"]):
    for file in files:
        if (
            os.path.splitext(file)[1] == ".hh"
            and os.path.splitext(file)[0] != "parameters"
        ):
            parameters.extend(getParamsFromFile(os.path.join(root, file),errorLog))

with open("parameterlist_error.log","w") as f:
    f.writelines(errorLog)

# make sorted dictionary of the entries
# treat duplicates (could have differing default values or type names - e.g. via aliases)
parameterDict = {}
for params in parameters:
    key = params["paramName"]
    if key in parameterDict:
        parameterDict[key]["defaultValue"].append(params["defaultValue"])
        parameterDict[key]["paramType"].append(params["paramType"])
    else:
        parameterDict[key] = params
        parameterDict[key]["defaultValue"] = [params["defaultValue"]]
        parameterDict[key]["paramType"] = [params["paramType"]]
parameterDict = {key: value for key, value in sorted(parameterDict.items())}

# determine actual entries (from duplicates)
tableEntryData = {}
multiValueLog = []
for key in parameterDict:

    entry = parameterDict[key]
    hasGroup = True if entry["paramName"].count(".") != 0 else False
    groupEntry = "-" if not hasGroup else entry["paramName"].split(".")[0]
    paramName = (
        entry["paramName"] if not hasGroup else entry["paramName"].partition(".")[2]
    )

    # In case of multiple occurrences, we use the first entry that is not None and print the others for possible manual editing
    paramType = entry["paramType"][0]
    defaultValue = next((e for e in entry["defaultValue"] if e), "-")

    hasMultiplePT = (
        True if not all(pt == paramType for pt in entry["paramType"]) else False
    )
    hasMultipleDV = (
        True
        if not all(
            dv == (defaultValue if defaultValue != "-" else None)
            for dv in entry["defaultValue"]
        )
        else False
    )
    if hasMultiplePT or hasMultipleDV:
        multiValueLog.append(
            "\nFound multiple occurrences of parameter "
            + paramName
            + " with differing specifications: " + '\n'
        )
    if hasMultiplePT:
        multiValueLog.append(" -> Specified type names:" + '\n')
        for typeName in set(entry["paramType"]):
            multiValueLog.append(" " * 8 + typeName + '\n')
        multiValueLog.append(
            " ---> For the parameters list, " + '\n'
            + " " * 8 + paramType + '\n'
            + " has been chosen. Please adapt manually if desired." + '\n'
        )

    if hasMultipleDV:
        multiValueLog.append(" -> Specified default values:" + '\n' )
        for default in set(entry["defaultValue"]):
            multiValueLog.append(" " * 8 + (default if default else "- (none given)") + '\n' )
        multiValueLog.append(
            " ---> For the parameters list," + '\n'
            + " " * 8 + defaultValue  + '\n'
            + " has been chosen. Please adapt manually if desired."  + '\n'
        )

    paramKey = groupEntry + "." + paramName
    explanation = (
        "" if paramKey not in oldDataDict else oldDataDict[paramKey]["Explanation"]
    )

    keysList = ["Group", "Parameter", "Type", "Default Value", "Explanation"]
    valuesList = [groupEntry,
                  paramName,
                  paramType,
                  defaultValue,
                  explanation]
    tableEntryData.update({paramKey:dict(zip(keysList, valuesList))})

#save the output of multiple values
with open("parameterlist_multivalue.log","w") as f:
    f.writelines(multiValueLog)

#make the json file readable
entryDataJson = json.dumps(tableEntryData, indent = 4, ensure_ascii= False, sort_keys = True)
with open('newparameters.json', 'w') as fp:
    fp.writelines(entryDataJson)

print("File newparameters.json is generated.")
