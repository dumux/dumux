#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightInfo: Copyright (c) 2020 Andrew Scott
# SPDX-License-Identifier: MIT
# The class TarSafe is adapted from https://github.com/beatsbears/tarsafe, licensed under MIT. All
# Changes are licensed under GPL-3.0-or-later.

"""
install external stuff for dumux
"""
import os
import shutil
import re
import urllib.request
import sys
import subprocess
import argparse
import textwrap
import pathlib
import tarfile as _unsafe_tarfile


# pylint: disable=C0103,W0212,W0622,C0116
class ChoicesAction(argparse._StoreAction):
    """Action to show choices in argparse"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.choices is None:
            self.choices = []
        self._choices_actions = []

    def add_choice(self, choice, help=""):
        self.choices.append(choice)
        choice_action = argparse.Action(option_strings=[], dest=choice, help=help)
        self._choices_actions.append(choice_action)

    def _get_subactions(self):
        return self._choices_actions


class TarSafe:
    """
    A safe subclass of the TarFile class for interacting with tar files.
    Runs all necessary checks for the safety of a tarfile (tar).

    Adapted from the code https://github.com/beatsbears/tarsafe
    which has the following license:

    MIT License

    Copyright (c) 2020 Andrew Scott

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    """

    def __init__(self, *args, **kwargs):
        # pylint: disable=consider-using-with
        self._tarfile = _unsafe_tarfile.TarFile(*args, **kwargs)
        self.directory = os.getcwd()

    @classmethod
    def open(cls, name=None, mode="r", fileobj=None, bufsize=_unsafe_tarfile.RECORDSIZE, **kwargs):
        """
        Open a tar archive
        """
        return _unsafe_tarfile.TarFile.open(name, mode, fileobj, bufsize, **kwargs)

    def extract(self, member, path="", set_attrs=True, *, numeric_owner=False):
        """
        Override the parent extract method and add safety checks.
        """
        self._safetar_check()
        self._tarfile.extract(member, path, set_attrs=set_attrs, numeric_owner=numeric_owner)

    def extractall(self, path=".", members=None, *, numeric_owner=False):
        """
        Override the parent extractall method and add safety checks.
        """
        self._safetar_check()
        self._tarfile.extractall(path, members, numeric_owner=numeric_owner)

    def _safetar_check(self):
        """
        Runs all necessary checks for the safety of a tarfile.
        """
        for tarinfo in iter(self._tarfile):
            if self._is_traversal_attempt(tarinfo=tarinfo):
                raise IOError(f"TarSafe: Attempted directory traversal for member: {tarinfo.name}")
            if self._is_unsafe_symlink(tarinfo=tarinfo):
                raise IOError(
                    f"TarSafe: Attempted directory traversal"
                    f" via symlink for member: {tarinfo.linkname}"
                )
            if self._is_unsafe_link(tarinfo=tarinfo):
                raise IOError(
                    f"TarSafe: Attempted directory traversal"
                    f" via link for member: {tarinfo.linkname}"
                )
            if tarinfo.ischr() or tarinfo.isblk():
                raise IOError("TarSafe: tarfile returns true for isblk() or ischr()")

    def _is_traversal_attempt(self, tarinfo):
        if not os.path.abspath(os.path.join(self.directory, tarinfo.name)).startswith(
            self.directory
        ):
            return True
        return False

    def _is_unsafe_symlink(self, tarinfo):
        if tarinfo.issym():
            symlink_file = pathlib.Path(
                os.path.normpath(os.path.join(self.directory, tarinfo.linkname))
            )
            if not os.path.abspath(os.path.join(self.directory, symlink_file)).startswith(
                self.directory
            ):
                return True
        return False

    def _is_unsafe_link(self, tarinfo):
        if tarinfo.islnk():
            link_file = pathlib.Path(
                os.path.normpath(os.path.join(self.directory, tarinfo.linkname))
            )
            if not os.path.abspath(os.path.join(self.directory, link_file)).startswith(
                self.directory
            ):
                return True
        return False


# pylint: enable=C0103,W0212,W0622,C0116


def showMessage(message):
    """Pretty print for info MESSAGES"""
    print("*" * 120)
    print(message)
    print("")
    print("*" * 120)


if len(sys.argv) == 1:
    showMessage(
        "No options given. For more information "
        "run the following command: \n ./installexternal.py --help"
    )
    sys.exit()

parser = argparse.ArgumentParser(
    prog="installexternal",
    usage="./installexternal.py [OPTIONS] PACKAGES",
    description="This script downloads extensions for dumux and dune \
                                     and installs some External Libraries and Modules.",
)
parser.register("action", "store_choice", ChoicesAction)
# Positional arguments
group = parser.add_argument_group(title="your choice of packages")
packages = group.add_argument("packages", nargs="+", metavar="PACKAGES", action="store_choice")
packages.add_choice(
    "dumux-extensions",
    help="Download dumux-course, dumux-lecture, dune-alugrid, dune-foamgrid and dune-subgrid.",
)
packages.add_choice(
    "dumux-course", help="Download dumux-course, dune-alugrid, dune-foamgrid and dune-subgrid."
)
packages.add_choice(
    "dune-extensions",
    help="Download dune-uggrid, dune-alugrid, dune-foamgrid, \
                    dune-subgrid, dune-spgrid, dune-mmesh and dune-functions.",
)
packages.add_choice("optimization", help="Download and install glpk and nlopt.")
packages.add_choice("others", help="Download and install opm and gstat.")
packages.add_choice("lecture", help="Download dumux-lecture.")
packages.add_choice("course", help="Download dumux-course.")
packages.add_choice("ug", help="Download dune-uggrid.")
packages.add_choice("alugrid", help="Download dune-alugrid.")
packages.add_choice("foamgrid", help="Download dune-foamgrid.")
packages.add_choice("subgrid", help="Download dune-subgrid.")
packages.add_choice("spgrid", help="Download dune-spgrid.")
packages.add_choice("mmesh", help="Download dune-mmesh.")
packages.add_choice("functions", help="Download dune-functions.")
packages.add_choice("glpk", help="Download and install glpk.")
packages.add_choice("nlopt", help="Download and install nlopt.")
packages.add_choice("opm", help="Download opm modules required for cornerpoint grids.")
packages.add_choice("gstat", help="Download and install gstat.")


# Optional arguments
options = parser.add_mutually_exclusive_group(required=False)
options.add_argument(
    "--clean", action="store_true", default=False, help="Delete all files for the given packages."
)
options.add_argument(
    "--download", action="store_true", default=False, help="Only download the packages."
)

parser.add_argument("--dune-branch", default="releases/2.9", help="Dune branch to be checked out.")
parser.add_argument(
    "--dumux-branch", default="releases/3.8", help="Dumux branch to be checked out."
)
parser.add_argument("--opm-branch", default="release/2022.10", help="Opm branch to be checked out.")
parser.add_argument("--mmesh-branch", default="release/1.4", help="Mmesh branch to be checked out.")

cmdArgs = vars(parser.parse_args())


def runCommand(command, currentDir="."):
    """Helper function to run commands with error checking and reporting"""
    with open(currentDir + "/installexternal.log", "a") as log:
        with subprocess.Popen(command, stdout=log, stderr=log, universal_newlines=True) as popen:
            returnCode = popen.wait()
            if returnCode:
                message = textwrap.dedent(
                    f"""
                    (Error) The command {command} returned with non-zero exit code
                    If you can't fix the problem yourself consider reporting your issue
                    on the mailing list (dumux@listserv.uni-stuttgart.de) and
                    attach the file 'installexternal.log'
                """
                )
                showMessage(message)
                sys.exit(1)


def gitClone(url, branch=None):
    """Clone a repository from a given URL"""
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    runCommand(command=[*clone, url])


def branchName(package, parameters):
    """Get the correct branch name"""
    # Set the branch
    if "dumux" in package:
        return parameters["dumux_branch"]
    if "mmesh" in package:
        return parameters["mmesh_branch"]
    if "dune" in package:
        return parameters["dune_branch"]
    if "opm" in package:
        return parameters["opm_branch"]
    return ""


def cleanPackage(package, finalMessage):
    """Clean up after a package"""
    if os.path.isfile(package + ".tar.gz"):
        os.remove(package + ".tar.gz")
    if os.path.exists(package):
        shutil.rmtree(package)
        finalMessage.append(f"{package} has been removed.")
    else:
        # Save message to be shown at the end
        finalMessage.append(f"The folder {package} does not exist.")


def filterPackageList(packageListOld):
    """Filter the package list and add possible dependencies"""
    packageList = []
    for pkg in packageListOld:
        if pkg in PACKAGE_NAMES:
            packageList.extend(PACKAGE_NAMES[pkg])
        else:
            packageList.extend([key for key in EXTERNAL_URLS if pkg in key])
    return packageList


def checkLocation():
    """check call location of this script"""
    if not os.path.isdir("dune-common"):
        showMessage(
            "You have to call " + sys.argv[0] + " for " + sys.argv[1] + " from\n"
            "the same directory in which dune-common is located.\n"
            "You cannot install it in this folder."
        )
        raise Exception("Script called in wrong location. Aborting.")


def installFromTarball(package, parameters, externalDir, finalMessage):
    """Install a package that uses a tarball as source code archive"""
    # Download the tarfile
    with urllib.request.urlopen(EXTERNAL_URLS[package]) as fileData:
        dataToWrite = fileData.read()
        with open(externalDir + "/" + package + ".tar.gz", "wb") as file:
            file.write(dataToWrite)

    # Save message to be shown at the end
    finalMessage.append(f"{package} has been successfully downloaded.")

    # Start Installation if the flag download is set to false.
    if not parameters["download"]:
        # Extract
        with TarSafe.open(package + ".tar.gz") as tarArchive:
            tarArchive.extractall()
            shutil.move(os.path.commonprefix(tarArchive.getnames()), package)  # rename

        # Start the configuration
        os.chdir(externalDir + "/" + package)
        if package == "gstat":
            with open("configure", "r+") as file:
                content = file.read()
                file.seek(0)
                file.truncate()
                file.write(content.replace("doc/tex/makefile", ""))

        # Run Configuration command
        configCmd = "./configure"
        runCommand(configCmd, currentDir=externalDir)
        try:
            runCommand("make", currentDir=externalDir)
        except subprocess.CalledProcessError as exc:
            raise Exception(f"{package} installation has failed.") from exc
        # Save message to be shown at the end
        if os.path.exists(externalDir + "/" + package):
            finalMessage.append(f"{package} has been successfully installed.")


def installExternal(parameters):
    """Main driver: install external packages"""

    topDir = os.getcwd()
    externalDir = topDir + "/external"
    parameters["packages"] = filterPackageList(parameters["packages"])

    # print the list of packages to be downloaded/installed/removed
    action = "removed" if parameters["clean"] else "downloaded"
    print(
        f"The following package(s) will be {action}:\n",
        ", ".join(parameters["packages"]),
        "\n",
    )

    checkLocation()

    # clear the log file
    logDir = externalDir if os.path.exists(externalDir) else topDir
    with open(logDir + "/installexternal.log", "w") as _:
        pass

    finalMessage = []
    for package in parameters["packages"]:
        os.chdir(topDir)
        # Package name for final message
        finalMessage.append("[---" + package + "---]")

        # Set the directory: create externalDir for external packages
        if not any(re.compile(p).match(package) for p in ["dumux", "dune", "opm"]):
            os.makedirs(externalDir, exist_ok=True)
            os.chdir(externalDir)

        branch = branchName(package, parameters)

        # Run the requested command
        if parameters["clean"]:
            cleanPackage(package, finalMessage)
            continue

        # Check if tarball
        tarball = EXTERNAL_URLS[package].endswith("tar.gz")

        if not os.path.exists(package):
            if tarball:
                installFromTarball(package, parameters, externalDir, finalMessage)
            else:
                # Clone from repo
                gitClone(EXTERNAL_URLS[package], branch)
                # Save message to be shown at the end
                finalMessage.append(f"{package} has been successfully cloned.")
        else:
            if tarball:
                finalMessage.append(f"{package} has been already installed.")
            else:
                # Checkout to the requested branch
                os.chdir(topDir + "/" + package)
                with subprocess.Popen(["git", "checkout", branch]) as _:
                    # Save message to be shown at the end
                    finalMessage.append(
                        f"-- Skip cloning {package}, because the folder already exists."
                    )
                    finalMessage.append(f"-- Checking out {package} " + branch)
                    continue

        # Save post installation message if there is any.
        if package in MESSAGES:
            finalMessage.extend(MESSAGES[package])

        # Change to topDir
        os.chdir(topDir)

    # Save post installation message about dunecontrol if need be.
    if not parameters["clean"] and any(
        x in pkg for pkg in parameters["packages"] for x in ["dumux", "dune", "opm"]
    ):
        finalMessage.append(
            "\n\nPlease run the following commands "
            "to clear the CMake cache, reconfigure and rebuild "
            "(can be copied to command line):\n\n  "
            "./dune-common/bin/dunecontrol bexec rm -r CMakeFiles CMakeCache.txt\n  "
            "./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts all"
        )

    # If cleanup and only logfile in the external directory, remove the directory
    if os.path.isdir(externalDir):
        _, _, files = next(os.walk(externalDir))
        if parameters["clean"] and len(files) == 1 and "installexternal.log" in files:
            shutil.rmtree(externalDir)

    return "\n".join(finalMessage)


#################################################################
#################################################################
# (1/2) Define the necessary packages and their URLS
#################################################################
#################################################################
DUNE_GIT_BASEURL = "https://gitlab.dune-project.org/"
DUMUX_GIT_BASEURL = "https://git.iws.uni-stuttgart.de/dumux-repositories/"
EXTERNAL_URLS = {
    "dumux-lecture": DUMUX_GIT_BASEURL + "dumux-lecture.git",
    "dumux-course": DUMUX_GIT_BASEURL + "dumux-course.git",
    "dune-uggrid": DUNE_GIT_BASEURL + "/staging/dune-uggrid.git",
    "dune-alugrid": DUNE_GIT_BASEURL + "extensions/dune-alugrid.git",
    "dune-foamgrid": DUNE_GIT_BASEURL + "extensions/dune-foamgrid.git",
    "dune-subgrid": DUNE_GIT_BASEURL + "extensions/dune-subgrid.git",
    "dune-spgrid": DUNE_GIT_BASEURL + "extensions/dune-spgrid.git",
    "dune-mmesh": DUNE_GIT_BASEURL + "samuel.burbulla/dune-mmesh.git",
    "dune-functions": DUNE_GIT_BASEURL + "staging/dune-functions.git",
    "dune-typetree": DUNE_GIT_BASEURL + "staging/dune-typetree.git",
    "glpk": "http://ftp.gnu.org/gnu/glpk/glpk-4.60.tar.gz",
    "nlopt": "http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz",
    "opm-common": "https://github.com/OPM/opm-common",
    "opm-grid": "https://github.com/OPM/opm-grid",
    "gstat": "http://gstat.org/gstat.tar.gz",
}

PACKAGE_NAMES = {
    "dumux-extensions": [
        "dumux-lecture",
        "dumux-course",
        "dune-alugrid",
        "dune-foamgrid",
        "dune-subgrid",
    ],
    "dumux-course": ["dumux-course", "dune-alugrid", "dune-foamgrid", "dune-subgrid"],
    "dune-extensions": [
        "dune-uggrid",
        "dune-alugrid",
        "dune-foamgrid",
        "dune-subgrid",
        "dune-spgrid",
        "dune-mmesh",
        "dune-functions",
        "dune-typetree",
    ],
    "functions": ["dune-functions", "dune-typetree"],
    "optimization": ["glpk", "nlopt"],
    "others": ["opm-common", "opm-grid", "gstat"],
}

MESSAGES = {
    "glpk": [
        "In addition, it might be necessary to set manually",
        "the glpk path in the CMAKE_FLAGS section of the .opts-file:",
        "  -DGLPK_ROOT=/path/to/glpk \\",
    ],
    "dune-mmesh": [
        "Maybe you also have to install CGAL",
        "(see cgal.org/download.html)",
        "Maybe you also need to change your core dune modules' branches to their newest versions.",
    ],
    "opm-common": [
        "In addition, it might be necessary to set manually some",
        "CMake variables in the CMAKE_FLAGS section of the .opts-file:",
        "  -DUSE_MPI=ON",
        "It might also be required to apply a patch,",
        "have a look at patches/README.md in the dumux folder.",
        "Currently, compiling opm with clang is not possible.",
        "",
        "Maybe you also have to install the following packages (see the",
        " opm prerequisites at opm-project.org):",
        "  BLAS, LAPACK, Boost, SuiteSparse, Zoltan",
    ],
}


#################################################################
#################################################################
# (2/2) Download/config/clean the requested packages
#################################################################
#################################################################
# Start download/configuration/cleaning tasks
showMessage(installExternal(cmdArgs))
