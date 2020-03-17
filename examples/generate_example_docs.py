#!/usr/bin/env python3
import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", help="The folder to look for examples", default=".")
parser.add_argument("-m", "--markdowngenerator", help="The markdown generator script taking a list of files", default="../bin/doc/cpp_to_md.sh")
args = vars(parser.parse_args())

def convertToMarkdownAndMerge(dir, includeList):
    script = os.path.join(os.path.abspath(args["markdowngenerator"]))
    with open(os.path.join(dir, "README.md"), "w") as readme:
        for include in includeList:
            if os.path.splitext(include)[1] == ".md":
                with open(include, "r") as markdown:
                    readme.write(markdown.read())
            else:
                readme.write("\n\n## The file `{}`\n\n".format(os.path.split(include)[1]))
                markdown = subprocess.check_output(["bash", script, include], encoding="utf-8")
                readme.write(markdown + "\n")

def generateReadme(dir):
    includeList = None
    try:
        configname = os.path.join(dir, ".doc_config")
        with open(configname, 'r') as config:
            includeList = [os.path.join(dir, include) for include in config.read().splitlines()]
    except FileNotFoundError:
        print("Error: The example directory {} does not contain a .doc_config file! Could not generate README.md!".format(dir))
        raise
    if includeList is not None:
        convertToMarkdownAndMerge(dir, includeList)

for path in os.listdir(args["directory"]):
    abspath = os.path.join(os.path.abspath(args["directory"]), path)
    if os.path.isdir(abspath):
        generateReadme(abspath)
