#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


"""
One click install script for dumux
"""
import os
import sys
import argparse
import itertools
import logging
import subprocess
import time
from shutil import which

logger = logging.getLogger("installdumux")
logger.setLevel(logging.DEBUG)

streamHandler = logging.StreamHandler(stream=sys.stdout)
streamHandler.setLevel(logging.INFO)
streamHandler.setFormatter(logging.Formatter("-- %(levelname)s: %(message)s"))
logger.addHandler(streamHandler)

fileHandler = logging.FileHandler("installdumux.log", encoding="utf-8")
fileHandler.setLevel(logging.DEBUG)
logger.addHandler(fileHandler)

_SPINNER = itertools.cycle(range(4))


class _Version:
    def __init__(self, version: str) -> None:
        self._version = [int(v) for v in version.strip(" ").strip("\n").split(".")]

    def __lt__(self, other) -> bool:
        for versionSelf, versionOther in zip(self._version, other._version):
            if versionSelf < versionOther:
                return True
            if versionSelf > versionOther:
                return False
        return False


parser = argparse.ArgumentParser(
    prog="installdumux",
    usage="./installdumux.py [OPTIONS]",
    description="This script downloads and compiles the latest release of Dumux.",
)
# Optional arguments
parser.add_argument("--dune-version", default="2.9", help="Dune version to be checked out.")
parser.add_argument("--dumux-version", default="3.8", help="Dumux version to be checked out.")
args = vars(parser.parse_args())

duneBranch = (
    args["dune_version"] if args["dune_version"] == "master" else "releases/" + args["dune_version"]
)
dumuxBranch = (
    args["dumux_version"]
    if args["dumux_version"] == "master"
    else "releases/" + args["dumux_version"]
)


def showMessage(message):
    """Pretty print message"""
    logger.info("*" * 120)
    logger.info(message)
    logger.info("*" * 120)


def spin(condition, *arguments, **keyWordArguments):
    """Display a progress symbol while condition is true"""
    while condition(*arguments, **keyWordArguments):
        sys.stdout.write("༄ " * next(_SPINNER))
        sys.stdout.flush()
        time.sleep(0.2)
        sys.stdout.write("\r\033[K")


def checkCppVersion():
    """Check compiler version"""
    requiredversion = "9.3"
    result = subprocess.check_output(["g++", "-dumpversion"]).decode().strip()
    if _Version(result) < _Version(requiredversion):
        logger.error("An error occurred while checking for prerequisites.")
        logger.error(f"g++ greater than or equal to {requiredversion} is required.")
        sys.exit(1)


def runCommand(command, workdir="."):
    """Run command with error checking"""
    logger.debug(f"Running {' '.join(command)}")
    with open("../installdumux.log", "a") as log:
        with subprocess.Popen(
            command,
            stdout=log,
            stderr=log,
            universal_newlines=True,
            cwd=workdir,
        ) as popen:
            spin(lambda p: p.poll() is None, popen)
            returnCode = popen.wait()
            if returnCode:
                logger.error(f"The command {' '.join(command)} returned with non-zero exit code.")
                logger.error("You find the error message in the file 'installdumux.log'.")
                logger.error("If you can't fix the problem yourself consider reporting your issue")
                logger.error("on the mailing list (dumux@listserv.uni-stuttgart.de) and")
                logger.error("attach the file 'installdumux.log'.")
                sys.exit(1)


def gitClone(url, branch=None):
    """Clone git repo"""
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    runCommand(command=[*clone, url])


def gitSetBranch(folder, branch):
    """Checkout specific git branch"""
    checkout = ["git", "checkout", branch]
    runCommand(command=checkout, workdir=folder)


# clear the log file
with open("installdumux.log", "w") as _:
    pass

#################################################################
#################################################################
# (1/3) Check some prerequisites
#################################################################
#################################################################
programs = ["git", "gcc", "g++", "cmake", "pkg-config"]
showMessage(f"(1/3) Checking all prerequisites: {' '.join(programs)}...")

# check some prerequisites
for program in programs:
    if which(program) is None:
        logger.error("An error occurred while checking for prerequisites.")
        logger.error(f"The required program '{program}' has not been found.")
        sys.exit(1)

if which("paraview") is None:
    logger.warning("ParaView could not be found. (You might have it but we can't find it.)")
    logger.warning("We recommend installing ParaView to view simulation results.")

checkCppVersion()

showMessage("(1/3) Step completed. All prerequisites found.")

#################################################################
#################################################################
# (2/3) Clone modules
#################################################################
#################################################################
# make a new folder containing everything
os.makedirs("./dumux", exist_ok=True)
os.chdir("dumux")

showMessage(
    "(2/3) Cloning repositories. This may take a while. "
    "Make sure to be connected to the internet..."
)

# the core modules
for module in ["common", "geometry", "grid", "localfunctions", "istl"]:
    if not os.path.exists(f"dune-{module}"):
        gitClone(f"https://gitlab.dune-project.org/core/dune-{module}.git", duneBranch)
    else:
        logger.info(f"Skip cloning dune-{module} because the folder already exists.")
        gitSetBranch(f"dune-{module}", duneBranch)

# dumux
if not os.path.exists("dumux"):
    gitClone("https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git", dumuxBranch)
else:
    logger.info("Skip cloning dumux because the folder already exists.")
    gitSetBranch("dumux", dumuxBranch)


showMessage("(2/3) Step completed. All repositories have been cloned into a containing folder.")

#################################################################
#################################################################
# (3/3) Configure and build
#################################################################
#################################################################
showMessage(
    "(3/3) Configure and build dune modules and dumux using dunecontrol. "
    "This may take several minutes..."
)

# run dunecontrol
runCommand(command=["./dune-common/bin/dunecontrol", "--opts=dumux/cmake.opts", "all"])

showMessage("(3/3) Step completed. Successfully configured and built dune and dumux.")

#################################################################
#################################################################
# Show message how to check that everything works
#################################################################
#################################################################
TEST_PATH = "dumux/dumux/build-cmake/test/porousmediumflow/1p"
if dumuxBranch == "master" or _Version(args["dumux_version"]) > _Version("3.3"):
    TEST_PATH += "/isothermal"
else:
    TEST_PATH += "/implicit/isothermal"

showMessage("༄  DuMuˣ installation complete 🎉")
print(
    "To test if everything works, "
    "please run the following commands (can be copied to command line):\n\n"
    f"  cd {TEST_PATH}\n"
    "  make test_1p_tpfa\n"
    "  ./test_1p_tpfa\n"
    "  paraview *pvd\n"
)
