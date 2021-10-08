#!/usr/bin/env python3

"""
Create a cff file (citation file format) for dumux.
"""

import sys
import yaml
from pathlib import Path

try:
    sys.path.append(str(Path(__file__).absolute().parents[1]))
    from util.common import runCommand, indent
    from util.moduleinfo import getModuleInfo
except Exception:
    sys.stderr.write("Failed to import utility modules\n")
    sys.exit(1)

# TODO: detect duplicates?
# TODO: filter with a minimum number of commits?
# TODO: Author ORCID?
# TODO: List only a few authors? How to select?
# TODO: If person has several first and/or last names, autodetection fails
# TODO: Umlaute are destryoed in the yml file
# TODO: Nicknames like 'NED' are in the yml file


# TODO: port getContributors.sh to python and provide these functions somewhere else?
def getFirstCommitSha():
    firstCommits = runCommand("git rev-list --max-parents=0 HEAD")
    return firstCommits.split('\n')[0]


def getCurrentCommitSha():
    return runCommand("git rev-parse HEAD")


def getCommitDate(commitSha):
    formatPlaceHolder = "%cs"
    return runCommand(f"git show --format={formatPlaceHolder}").strip("\n")


def getDumuxVersion():
    modPath = Path(__file__).absolute().parents[2]
    return getModuleInfo(modPath, "Version")


def getCFFBaseConfig():
    version = getDumuxVersion()
    date = getCommitDate(getCurrentCommitSha())
    return {
        "cff-version": "1.2.0",
        "title": "Dumux - an open-source simulator for "
                 "flow and transport in porous media",
        "version": version,
        "doi": "https://doi.org/10.5281/zenodo.2479594",
        "date-released": date
    }


def printStatus(message, file=sys.stdout):
    indentationWidth = 4
    indentation = " " * indentationWidth
    message = indent(message, indentation)
    message = " -- " + message[indentationWidth:]
    print(f"{message}", file=file)


def getAllAuthorsAndNumberOfCommits():
    """Return the list of authors that contributed to dumux so far"""
    filePath = Path(__file__).absolute().parents[0]
    scriptPath = filePath / Path("getcontributors.sh")
    authorList = runCommand(
        "{} -from {} -to {}"
        .format(str(scriptPath), getFirstCommitSha(), getCurrentCommitSha())
    )
    authorList = authorList.strip('\n')
    authorList = authorList.split('\n')[1:]
    return [
        {"name": a.split(maxsplit=1)[1], "numCommits": a.split(maxsplit=1)[0]}
        for a in authorList
    ]


def authorHasLastName(author):
    return len(author["name"].split(" ", maxsplit=1)) > 1


def removeAuthors(authors, toRemove):
    for author in toRemove:
        authors.remove(author)


class CFFWriter:
    def __init__(self, baseConfig):
        self.baseConfig = baseConfig
        self.authors = []

    def addAuthors(self, authors):
        cffAuthors = []
        for author in authors:
            firstName, lastName = author["name"].split(maxsplit=1)
            cffAuthors.append({
                "family-names": lastName,
                "given-names": firstName
            })
        self.authors.extend(cffAuthors)

    def write(self, fileName):
        with open(fileName, "w") as ymlFile:
            yaml.dump(self.baseConfig, ymlFile)
            yaml.dump({"authors": self.authors}, ymlFile)


if __name__ == "__main__":
    authors = getAllAuthorsAndNumberOfCommits()
    printStatus("number of all authors found: {}".format(len(authors)))

    singleNameAuthors = filter(lambda a: not authorHasLastName(a), authors)
    singleNameAuthors = list(singleNameAuthors)
    printStatus("Single-name authors to be exlcuded:\n" + '\n'.join(
        a["name"] for a in singleNameAuthors
    ))
    removeAuthors(authors, singleNameAuthors)
    printStatus("number of remaining authors: {}".format(len(authors)))

    writer = CFFWriter(getCFFBaseConfig())
    writer.addAuthors(authors)
    writer.write("cff.yml")
