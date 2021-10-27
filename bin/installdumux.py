#!/usr/bin/env python3

"""
One click install script for dumux
"""
import os
import sys
import argparse
import subprocess
import textwrap
from distutils.spawn import find_executable
from distutils.version import LooseVersion

parser = argparse.ArgumentParser(
    prog="installdumux",
    usage="./installdumux.py [OPTIONS]",
    description="This script downloads and compiles the latest release of Dumux.",
)
# Optional arguments
parser.add_argument("--dune-version", default="2.8", help="Dune version to be checked out.")
parser.add_argument("--dumux-version", default="3.4", help="Dumux version to be checked out.")
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
    print("*" * 120)
    print(message)
    print("*" * 120)


def checkCppVersion():
    """Check compiler version"""
    requiredversion = "7"
    result = subprocess.check_output(["g++", "-dumpversion"]).decode().strip()
    if LooseVersion(result) < LooseVersion(requiredversion):
        print("-- An error occured while checking for prerequistes.")
        raise Exception(
            f"g++ greater than or equal to {requiredversion} "
            "is required for dumux releases >=3.2!"
        )


def runCommand(command, workdir="."):
    """Run command with error checking"""
    with open("../installdumux.log", "a") as log:
        with subprocess.Popen(
            command,
            stdout=log,
            stderr=log,
            universal_newlines=True,
            cwd=workdir,
        ) as popen:
            returnCode = popen.wait()
            if returnCode:
                message = textwrap.dedent(
                    f"""\

                    (Error) The command {command} returned with non-zero exit code
                    If you can't fix the problem yourself consider reporting your issue
                    on the mailing list (dumux@listserv.uni-stuttgart.de) and
                    attach the file 'installdumux.log'
                    """
                )
                showMessage(message)
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
# (1/3) Check some prerequistes
#################################################################
#################################################################
programs = ["git", "gcc", "g++", "cmake", "pkg-config"]
showMessage("(1/3) Checking all prerequistes: " + " ".join(programs) + "...")

# check some prerequistes
for program in programs:
    if find_executable(program) is None:
        print("-- An error occured while checking for prerequistes.")
        raise Exception(f"Program {program} has not been found.")

if find_executable("paraview") is None:
    print(
        "-- Warning: paraview seems to be missing. You may not be able to view simulation results!"
    )

checkCppVersion()

showMessage("(1/3) Step completed. All prerequistes found.")

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
        print(f"-- Skip cloning dune-{module} because the folder already exists.")
        gitSetBranch(f"dune-{module}", duneBranch)

# dumux
if not os.path.exists("dumux"):
    gitClone("https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git", dumuxBranch)
else:
    print("-- Skip cloning dumux because the folder already exists.")
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

showMessage("(3/3) Step completed. Succesfully configured and built dune and dumux.")

#################################################################
#################################################################
# Show message how to check that everything works
#################################################################
#################################################################
TEST_PATH = "dumux/dumux/build-cmake/test/porousmediumflow/1p"
if dumuxBranch == "master" or LooseVersion(args["dumux_version"]) > LooseVersion("3.3"):
    TEST_PATH += "/isothermal"
else:
    TEST_PATH += "/implicit/isothermal"

showMessage(
    "(Installation complete) To test if everything works, "
    "please run the following commands (can be copied to command line):\n\n"
    f"  cd {TEST_PATH}\n"
    "  make test_1p_tpfa\n"
    "  ./test_1p_tpfa\n"
    "  paraview *pvd\n"
)
