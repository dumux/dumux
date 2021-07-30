#!/usr/bin/env python3

"""
One click install script for dumux
"""
import os
import sys
import argparse
import subprocess
from distutils.spawn import find_executable
from distutils.version import LooseVersion

parser = argparse.ArgumentParser(
    prog="installdumux",
    usage="./installdumux.py [OPTIONS]",
    description="This script downloads and compiles the latest release of Dumux.",
)
# Optional arguments
parser.add_argument("--dune-version", default="2.7", help="Dune version to be checked out.")
parser.add_argument("--dumux-version", default="3.4", help="Dumux version to be checked out.")
args = vars(parser.parse_args())

dune_branch = (
    args["dune_version"] if args["dune_version"] == "master" else "releases/" + args["dune_version"]
)
dumux_branch = (
    args["dumux_version"]
    if args["dumux_version"] == "master"
    else "releases/" + args["dumux_version"]
)


def show_message(message):
    print("*" * 120)
    print(message)
    print("*" * 120)


def check_cpp_version():
    requiredversion = "7"
    result = subprocess.check_output(["g++", "-dumpversion"]).decode().strip()
    if LooseVersion(result) < LooseVersion(requiredversion):
        print("-- An error occured while checking for prerequistes.")
        raise Exception(
            "g++ greater than or equal to {} is required for dumux releases >=3.2!".format(
                requiredversion
            )
        )


def run_command(command, workdir="."):
    with open("../installdumux.log", "a") as log:
        popen = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            cwd=workdir,
        )
        for line in popen.stdout:
            log.write(line)
            print(line, end="")
        for line in popen.stderr:
            log.write(line)
            print(line, end="")
        popen.stdout.close()
        popen.stderr.close()
        return_code = popen.wait()
        if return_code:
            print("\n")
            message = "\n    (Error) The command {} returned with non-zero exit code\n".format(
                command
            )
            message += "\n    If you can't fix the problem yourself consider reporting your issue\n"
            message += "    on the mailing list (dumux@listserv.uni-stuttgart.de) and attach the file 'installdumux.log'\n"
            show_message(message)
            sys.exit(1)


def git_clone(url, branch=None):
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    run_command(command=[*clone, url])


def git_setbranch(folder, branch):
    checkout = ["git", "checkout", branch]
    run_command(command=checkout, workdir=folder)


# clear the log file
open("installdumux.log", "w").close()

#################################################################
#################################################################
# (1/3) Check some prerequistes
#################################################################
#################################################################
programs = ["git", "gcc", "g++", "cmake", "pkg-config"]
show_message("(1/3) Checking all prerequistes: " + " ".join(programs) + "...")

# check some prerequistes
for program in programs:
    if find_executable(program) is None:
        print("-- An error occured while checking for prerequistes.")
        raise Exception("Program {} has not been found.".format(program))

if find_executable("paraview") is None:
    print(
        "-- Warning: paraview seems to be missing. You may not be able to view simulation results!"
    )

check_cpp_version()

show_message("(1/3) Step completed. All prerequistes found.")

#################################################################
#################################################################
# (2/3) Clone modules
#################################################################
#################################################################
# make a new folder containing everything
os.makedirs("./dumux", exist_ok=True)
os.chdir("dumux")

show_message(
    "(2/3) Cloning repositories. This may take a while. Make sure to be connected to the internet..."
)

# the core modules
for module in ["common", "geometry", "grid", "localfunctions", "istl"]:
    if not os.path.exists("dune-{}".format(module)):
        git_clone("https://gitlab.dune-project.org/core/dune-{}.git".format(module), dune_branch)
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))
        git_setbranch("dune-{}".format(module), dune_branch)

# dumux
if not os.path.exists("dumux"):
    git_clone("https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git", dumux_branch)
else:
    print("-- Skip cloning dumux because the folder already exists.")
    git_setbranch("dumux", dumux_branch)


show_message("(2/3) Step completed. All repositories have been cloned into a containing folder.")

#################################################################
#################################################################
# (3/3) Configure and build
#################################################################
#################################################################
show_message(
    "(3/3) Configure and build dune modules and dumux using dunecontrol. This may take several minutes..."
)

# run dunecontrol
run_command(command=["./dune-common/bin/dunecontrol", "--opts=dumux/cmake.opts", "all"])

show_message("(3/3) Step completed. Succesfully configured and built dune and dumux.")

#################################################################
#################################################################
# Show message how to check that everything works
#################################################################
#################################################################
test_path = "dumux/dumux/build-cmake/test/porousmediumflow/1p"
if dumux_branch == "master" or LooseVersion(args["dumux_version"]) > LooseVersion("3.3"):
    test_path += "/isothermal"
else:
    test_path += "/implicit/isothermal"

show_message(
    "(Installation complete) To test if everything works, please run the following commands (can be copied to command line):\n\n"
    "  cd {}\n"
    "  make test_1p_tpfa\n"
    "  ./test_1p_tpfa\n"
    "  paraview *pvd\n".format(test_path)
)
