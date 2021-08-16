#!/usr/bin/env python3

"""
Automatically updates parameterlist.txt by searching all *.hh files
for usage of getParam or getParamFromGroup.
"""

import os
import argparse
import json

parser = argparse.ArgumentParser(
    description="""
    Generate a json file of parameters list from header files.
    ----------------------------------------------------------------
    If input file is given, the descriptions will be copied from the input.
    If no input file is given, the descriptions will be empty.
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

class CheckExistAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        if os.path.isfile(values):
            setattr(namespace, self.dest, values)
            setattr(namespace, "hasInput", True)
        else:
            parser.error("File {} does not exist!".format(values))

parser.add_argument(
    "--input",
    help="json file of given paratemers",
    action=CheckExistAction,
    metavar="input",
    dest="inputFile",
    required=True,
)
parser.add_argument(
    "--root",
    help="root path of Dumux",
    metavar="rootpath",
    default=os.path.abspath(os.path.join(os.path.abspath(__file__), "../../../")),
)
parser.add_argument(
    "--output",
    help="relative path (to the root path) of the output file",
    metavar="output",
    default="doc/doxygen/extradoc/newparameterlist.txt",
)

args = vars(parser.parse_args())

# get the data from parameter json file
parameterDict = {}
with open(args["inputFile"], "r") as file:
    data = file.read()
    parameterDict = json.loads(data)


# determine the column widths
maxGroupWidth = 0
maxParamWidth = 0
maxTypeWidth = 0
maxDefaultWidth = 0
maxExplanationWidth = 0
tableEntryData = {}

for key in parameterDict:

    groupEntry = parameterDict[key]["Group"]
    paramName = parameterDict[key]["Parameter"]
    paramType = parameterDict[key]["Type"]
    defaultValue = parameterDict[key]["Default Value"]
    explanation = parameterDict[key]["Explanation"]

    # +3 because \b will be added later
    maxGroupWidth = max(maxGroupWidth, len(groupEntry) + 3)
    maxParamWidth = max(maxParamWidth, len(paramName))
    maxTypeWidth = max(maxTypeWidth, len(paramType))
    maxDefaultWidth = max(maxDefaultWidth, len(defaultValue))
    maxExplanationWidth = max(maxExplanationWidth, len(explanation))

# generate actual table entries
tableEntriesWithGroup = []
tableEntriesWithoutGroup = []
previousGroupEntry = None

for key in parameterDict:

    groupEntry = parameterDict[key]["Group"]
    paramName = parameterDict[key]["Parameter"]
    paramType = parameterDict[key]["Type"]
    defaultValue = parameterDict[key]["Default Value"]
    explanation = parameterDict[key]["Explanation"]

    if groupEntry != previousGroupEntry:
        previousGroupEntry = groupEntry
        if groupEntry != "-":
            groupEntry = "\\b " + groupEntry

    tableEntry = " * | {} | {} | {} | {} | {} |".format(
        groupEntry.ljust(maxGroupWidth),
        paramName.ljust(maxParamWidth),
        paramType.ljust(maxTypeWidth),
        defaultValue.ljust(maxDefaultWidth),
        explanation.ljust(maxExplanationWidth),
    )

    if groupEntry != "-":
        tableEntriesWithGroup.append(tableEntry)
    else:
        tableEntriesWithoutGroup.append(tableEntry)

# combine entries
tableEntries = tableEntriesWithoutGroup + tableEntriesWithGroup

header = """/*!
 *\\file
 *\ingroup Parameter
 *
 *\\brief List of currently useable run-time parameters
 *
 * The listed run-time parameters are available in general,
 * but we point out that a certain model might not be able
 * to use every parameter!
 *\n"""
header += " * | " + "Group".ljust(maxGroupWidth)
header +=   " | " + "Parameter".ljust(maxParamWidth)
header +=   " | " + "Type".ljust(maxTypeWidth)
header +=   " | " + "Default Value".ljust(maxDefaultWidth)
header +=   " | " + "Explanation".ljust(maxExplanationWidth) + "|\n"

header += " * | " + ":-".ljust(maxGroupWidth)
header +=   " | " + ":-".ljust(maxParamWidth)
header +=   " | " + ":-".ljust(maxTypeWidth)
header +=   " | " + ":-".ljust(maxDefaultWidth)
header +=   " | " + ":-".ljust(maxExplanationWidth) + "|\n"

header += " * | " + "".ljust(maxGroupWidth)
header +=   " | " + "ParameterFile".ljust(maxParamWidth)
header +=   " | " + "std::string".ljust(maxTypeWidth)
header +=   " | " + "executable.input".ljust(maxDefaultWidth)
header +=   " | " + ":-".ljust(maxExplanationWidth) + "|\n"

# overwrite the old parameterlist.txt file
with open(os.path.join(args["root"], args["output"]), "w") as outputfile:
    outputfile.write(header)
    for e in tableEntries:
        outputfile.write(e + "\n")
    outputfile.write(" */\n")
    print(
        "The parameter list is written to {}".format(
            os.path.join(args["root"], args["output"])
        )
    )
