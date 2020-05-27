#!/usr/bin/env python3

"""
One click install script for dumux
"""
import os
import subprocess
from distutils.spawn import find_executable
from pkg_resources import parse_version

def show_message(message):
    print("*" * 120)
    print(message)
    print("*" * 120)


def check_cpp_version():
    requiredversion = "7"
    result = subprocess.check_output(["g++", "-dumpversion"]).decode().strip()
    if parse_version(result) < parse_version(requiredversion):
        print("-- An error occured while checking for prerequistes.")
        raise Exception("g++ greater than or equal to {} is required for dumux releases >=3.2!".format(requiredversion))


def git_clone(url, branch=None):
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    result = subprocess.run([*clone, url])
    if result.returncode != 0:
        print("-- Failed to clone the repositories. Look for repository-specific errors.")
        raise Exception(result.stderr)


#################################################################
#################################################################
## (1/3) Check some prerequistes
#################################################################
#################################################################
programs = ['wget', 'git', 'gcc', 'g++', 'cmake', 'pkg-config']
show_message("(1/3) Checking all prerequistes: " + " ".join(programs) + "...")

# check some prerequistes
for program in programs:
    if find_executable(program) is None:
        print("-- An error occured while checking for prerequistes.")
        raise Exception("Program {} has not been found.".format(program))

if find_executable('paraview') is None:
    print("-- Warning: paraview seems to be missing. You may not be able to view simulation results!")

check_cpp_version()

show_message("(1/3) Step completed. All prerequistes found.")

#################################################################
#################################################################
## (2/3) Clone modules
#################################################################
#################################################################
# make a new folder containing everything
os.makedirs("./DUMUX", exist_ok=True)
os.chdir("DUMUX")

show_message("(2/3) Cloning repositories. This may take a while. Make sure to be connected to the internet...")

dune_version=2.7
dumux_version=3.2
# the core modules
for module in ['common', 'geometry', 'grid', 'localfunctions', 'istl']:
    if not os.path.exists("dune-{}".format(module)):
        git_clone('https://gitlab.dune-project.org/core/dune-{}.git'.format(module), "releases/{}".format(dune_version))
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))

# dumux
if not os.path.exists("dumux"):
    git_clone('https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git', "releases/{}".format(dumux_version))
else:
    print("-- Skip cloning dumux because the folder already exists.")


show_message("(2/3) Step completed. All repositories have been cloned into a containing folder.")

#################################################################
#################################################################
## (3/3) Configure and build
#################################################################
#################################################################
show_message("(3/3) Configure and build dune modules and dumux using dunecontrol. This may take several minutes...")

# run dunecontrol
if not os.path.isfile("cmake.opts"):
    subprocess.run(["wget","https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/releases/3.2/cmake.opts"])
else:
    print("-- The file cmake.opts already exists. The existing file will be used to configure dumux.")


result = subprocess.run(["./dune-common/bin/dunecontrol", "--opts=cmake.opts", "all"])
if result.returncode != 0:
    print("-- Failed to build the dune libaries.")
    raise Exception(result.stderr)

show_message("(3/3) Step completed. Succesfully configured and built dune and dumux.")

#################################################################
#################################################################
## Show message how to check that everything works
#################################################################
#################################################################
show_message("(Installation complete) To test if everything works, please run the following commands (can be copied to command line):\n\n"
             "  cd DUMUX/dumux/build-cmake/test/porousmediumflow/1p/implicit/isothermal\n"
             "  make test_1p_tpfa\n"
             "  ./test_1p_tpfa\n"
             "  paraview *pvd\n")
