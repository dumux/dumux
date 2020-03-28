#!/usr/bin/env python3

import os
import argparse

from convert_code_to_doc import *

def convertToMarkdownAndMerge(dir, includeList):
    with open(os.path.join(dir, "README.md"), "w") as readme:
        for include, description in includeList:
            fileExtension = os.path.splitext(include)[1]
            if fileExtension == ".md":
                with open(os.path.join(dir, include), "r") as markdown:
                    readme.write(markdown.read())
            elif fileExtension == ".hh" or fileExtension == ".cc":
                title = description + " (`{}`)\n\n\n" if description else "The file `{}`\n\n\n"
                title = title.format(os.path.split(include)[1])
                readme.write("\n\n## " + title)
                with open(os.path.join(dir, include), "r") as cppCode:
                    readme.write(transformCode(cppCode.read(), cppRules()) + "\n")
            else:
                raise IOError("Unsupported or unknown file extension *{}".format(fileExtension))

def generateReadme(dir):
    includeList = None
    # look for .doc_config, if not found we pass
    try:
        configname = os.path.join(dir, ".doc_config")
        with open(configname, 'r') as config:
            def readline(line):
                line = line.split(' ', 1)
                if len(line) == 1:
                    return [line[0], ""]
                else:
                    return [line[0], line[1]]
            includeList = [readline(line) for line in config.read().splitlines()]
    except FileNotFoundError:
        pass
    if includeList is not None:
        convertToMarkdownAndMerge(dir, includeList)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="The folder to look for examples", default=".")
    args = vars(parser.parse_args())

    for path in os.listdir(args["directory"]):
        abspath = os.path.join(os.path.abspath(args["directory"]), path)
        if os.path.isdir(abspath):
            generateReadme(abspath)
