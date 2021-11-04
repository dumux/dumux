#!/usr/bin/env python3

"""
Automatically updates parameterlist.txt by searching all *.hh files
for usage of getParam or getParamFromGroup.
"""

import os
import argparse
import json


def getEnclosedContent(string, openKey, closeKey):
    """find the content of the given string between
    the first matching pair of opening/closing keys"""

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


def extractParamName(line):
    """extract a parameter from a given line"""

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

    # if we get an empty name
    if len(paramName) == 0:
        raise IOError("Could not correctly process parameter name")
    # if interior spaces occur in the parameter name, we can't identify it
    if paramName[0] != '"' or paramName[-1] != '"' or " " in paramName:
        raise IOError("Could not correctly process parameter name")

    return {"paramType": paramType, "paramName": paramName.strip('"'), "defaultValue": defaultValue}


def getParamsFromFile(file, log):
    """extract all parameters from a given file"""
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
        log.append(
            (
                "\n\n{} parameter{} in file {} could not be retrieved automatically. "
                "Please check them yourself:"
            ).format(len(errors), "s" if len(errors) > 1 else "", file)
        )
        for lineIdx in errors:
            log.append("\n\t-> line {}: {}".format(lineIdx, errors[lineIdx]["line"]))
            log.append("\n\t\t-> error message: {}".format(errors[lineIdx]["message"]))
        log.append(" ")

    return parameters


class CheckExistAction(argparse.Action):
    """check if the input file exists"""

    def __call__(self, parser, namespace, values, option_string):
        if os.path.isfile(values):
            setattr(namespace, self.dest, values)
            setattr(namespace, "hasInput", True)
        else:
            parser.error("File {} does not exist!".format(values))


argumentParser = argparse.ArgumentParser(
    description="""
    This script generates parameters list from header files.
    The header files "test" and "examples" folders are not included.
    ----------------------------------------------------------------
    If input file is given, the descriptions will be copied from the input.
    Multientry of parameters are allowed in input files.
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
argumentParser.add_argument(
    "--root",
    help="root path of Dumux",
    metavar="rootpath",
    default=os.path.abspath(os.path.join(os.path.abspath(__file__), "../../../")),
)
argumentParser.add_argument(
    "--input",
    help="json file of given paratemers",
    action=CheckExistAction,
    metavar="input",
    dest="inputFile",
    default="../../doc/doxygen/extradoc/parameters.json",
)
argumentParser.add_argument(
    "--output",
    help="relative path (to the root path) of the output file",
    metavar="output",
    default="doc/doxygen/extradoc/parameterlist.txt",
)
args = vars(argumentParser.parse_args())

# search all *.hh files for parameters
# TODO: allow runtime args with extensions and folder(s) to be checked
parameters = []
log = []
rootDir = args["root"]
for root, dirs, files in os.walk(rootDir):
    dirs[:] = [d for d in dirs if d not in ("test", "examples")]
    for file in files:
        if os.path.splitext(file)[1] == ".hh" and os.path.splitext(file)[0] != "parameters":
            parameters.extend(getParamsFromFile(os.path.join(root, file), log))

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

# get the explanations from current parameter json file
inputDict = {}
if args["inputFile"]:
    with open(args["inputFile"], "r") as f:
        data = f.read()
        inputDict = json.loads(data)

# add the missing parameters from input
missingParameters = [key for key in inputDict if key.replace("-.", "") not in parameterDict]

for missingkey in missingParameters:

    paramName = inputDict[missingkey]["group"] + "." + inputDict[missingkey]["parameter"]

    hasMode = "mode" in inputDict[missingkey].keys()
    mode = ""
    if hasMode: mode = inputDict[missingkey]["mode"]
    if  mode == "manual":
        key = missingkey.replace("-.", "")
        parameterDict[key] = inputDict[missingkey]
        parameterDict[key]["defaultValue"] = inputDict[missingkey]["default"]
        parameterDict[key]["paramType"] = inputDict[missingkey]["type"]
        parameterDict[key]["paramName"] = paramName

        log.append("\nAdd to parameter list: '"
        + paramName
        + "' which could not be extracted from code but is explicitly added in "
        + args["inputFile"] +  "\n")

    else:
        log.append("\nFound parameter in "
        + args["inputFile"] + ": '"
        + paramName
        + "', set mode to manuel if it is to be kept.\n")
parameterDict = dict(sorted(parameterDict.items(), key=lambda kv: kv[0]))
# determine actual entries (from duplicates)
# and determine maximum occurring column widths

tableEntryData = []
for key in parameterDict:

    entry = parameterDict[key]
    hasGroup = entry["paramName"].count(".") > 0
    groupEntry = "-" if not hasGroup else entry["paramName"].split(".")[0]
    paramName = entry["paramName"] if not hasGroup else entry["paramName"].partition(".")[2]

    # In case of multiple occurrences,
    # we prefer the default value from input
    # otherwise use the first entry that is not None
    # and write the others in log for possible manual editing
    # determin multiple entries in input
    keyInput = groupEntry + "." + paramName
    numOfEntries = 0
    if keyInput in inputDict:
        numOfEntries = max(
            len(inputDict[keyInput]["default"]),
            len(inputDict[keyInput]["type"]),
            len(inputDict[keyInput]["explanation"]),
        )

    paramType = entry["paramType"][0]
    defaultValue = next((e for e in entry["defaultValue"] if e), "-")

    hasMultiplePT = (not all(pt == paramType for pt in entry["paramType"]))
    hasMultipleDV = (
        not all(
            dv == (defaultValue if defaultValue != "-" else None) for dv in entry["defaultValue"]
        )
    )
    if hasMultiplePT or hasMultipleDV:
        log.append(
            "\n\nFound multiple occurrences of parameter "
            + paramName
            + " with differing specifications: "
        )
    if hasMultiplePT:
        log.append("\n -> Specified type names:")
        for typeName in entry["paramType"]:
            log.append("\n" + " " * 8 + typeName)
        if numOfEntries == 0:
            log.append(
                "\n ---> For the parameters list, "
                + paramType
                + " has been chosen. Please adapt manually if desired."
            )
        else:
            paramType = inputDict[keyInput]["type"]
            log.append(
                "\n ---> For the parameters list, "
                + str(paramType)
                + " has been chosen. Type is from input file."
            )
    if hasMultipleDV:
        log.append("\n -> Specified default values:")
        for default in entry["defaultValue"]:
            log.append("\n" + " " * 8 + (default if default else "- (none given)"))
        if numOfEntries == 0:
            log.append(
                "\n ---> For the parameters list, "
                + defaultValue
                + " has been chosen. Please adapt manually if desired."
            )
        else:
            defaultValue = inputDict[keyInput]["default"]
            log.append(
                "\n ---> For the parameters list, "
                + str(defaultValue)
                + " has been chosen. Default value is from input file."
            )

    # get explanation
    explanation = ""
    if numOfEntries > 0:
        explanation = inputDict[keyInput]["explanation"]
        paramType = inputDict[keyInput]["type"]
        defaultValue = inputDict[keyInput]["default"]

    if numOfEntries == 0:
        tableEntryData.append(
            {
                "group": groupEntry,
                "name": paramName,
                "type": paramType,
                "default": defaultValue,
                "explanation": explanation,
            }
        )
    else:
        for i in range(numOfEntries):
            # maybe fewer entries for some keys
            if len(defaultValue) < i + 1:
                defaultValue.append(defaultValue[i - 1])
            if len(explanation) < i + 1:
                explanation.append(explanation[i - 1])
            if len(paramType) < i + 1:
                paramType.append(paramType[i - 1])

            tableEntryData.append(
                {
                    "group": groupEntry,
                    "name": paramName,
                    "type": paramType[i],
                    "default": defaultValue[i],
                    "explanation": explanation[i],
                }
            )

# generate actual table entries
tableEntriesWithGroup = []
tableEntriesWithoutGroup = []
previousGroupEntry = None

groupEntryLength = 20
paramNameLength = 45
paramTypeLength = 24
defaultValueLength = 15
explanationLength = 80
for data in tableEntryData:

    groupEntry = data["group"]
    paramName = data["name"]
    paramType = data["type"]
    defaultValue = data["default"]
    explanation = data["explanation"]

    if groupEntry != previousGroupEntry:
        previousGroupEntry = groupEntry
        if groupEntry != "-":
            groupEntry = "\\b " + groupEntry

    if len(explanation.strip()) == 0:
        log.append("\n parameter " + groupEntry + "." + paramName + " has no explanation.")

    tableEntry = " * | {} | {} | {} | {} | {} |".format(
        groupEntry.ljust(groupEntryLength),
        paramName.ljust(paramNameLength),
        paramType.ljust(paramTypeLength),
        defaultValue.ljust(defaultValueLength),
        explanation.ljust(explanationLength)
    )

    if groupEntry != "-":
        tableEntriesWithGroup.append(tableEntry)
    else:
        tableEntriesWithoutGroup.append(tableEntry)

# combine entries
tableEntries = tableEntriesWithoutGroup + tableEntriesWithGroup

header = """/*!
 *\\internal
   _    _                  _
  | |  | |                (_)
  | |  | | __ _ _ __ _ __  _ _ __   __ _
  | |/\| |/ _` | '__| '_ \| | '_ \ / _` |
  \  /\  / (_| | |  | | | | | | | | (_| |
   \/  \/ \__,_|_|  |_| |_|_|_| |_|\__, |
                                    __/ |
                                   |___/

 * WARNING: This file is auto-generated. Do not manually edit.
 * Run bin/doc/generate_parameterlist.py
 *\\endinternal
 *
 *\\file
 *\\ingroup Parameter
 *
 *\\brief List of currently useable run-time parameters
 *
 * The listed run-time parameters are available in general,
 * but we point out that a certain model might not be able
 * to use every parameter!
 *\n"""
header += " * | " + "Group".ljust(groupEntryLength)
header += " | " + "Parameter".ljust(paramNameLength)
header += " | " + "Type".ljust(paramTypeLength)
header += " | " + "Default Value".ljust(defaultValueLength)
header += " | Explanation ".ljust(explanationLength) + "|\n"

header += " * | " + ":-".ljust(groupEntryLength)
header += " | " + ":-".ljust(paramNameLength)
header += " | " + ":-".ljust(paramTypeLength)
header += " | " + ":-".ljust(defaultValueLength)
header += " | :- ".ljust(explanationLength) + "|\n"

header += " * | " + "-".ljust(groupEntryLength)
header += " | " + "ParameterFile".ljust(paramNameLength)
header += " | " + "std::string".ljust(paramTypeLength)
header += " | " + "executable.input".ljust(defaultValueLength)
header += " | :- ".ljust(explanationLength) + "|\n"

# overwrite the old parameterlist.txt file
with open(os.path.join(rootDir, args["output"]), "w") as outputfile:
    outputfile.write(header)
    for e in tableEntries:
        outputfile.write(e + "\n")
    outputfile.write(" */\n")
    print("Finished, check the output file and log.")

with open("generate_parameterlist.log", "w") as f:
    f.writelines(log)
