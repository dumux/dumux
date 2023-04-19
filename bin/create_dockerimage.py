#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


"""
Script to create a Docker image from a Dune module
"""

import os
import sys
import string
import shutil
import argparse
import subprocess
from util.moduleinfo import getModuleFile
from util.moduleinfo import extractModuleInfos
from util.common import queryYesNo

# require python 3
if sys.version_info[0] < 3:
    sys.exit("\nERROR: Python3 required")


def substituteAndWrite(template, target, mapping):
    """substitute content from template and write to target"""
    if not os.path.exists(template):
        sys.exit("Template file '" + template + "' could not be found")
    with open(target, "w") as targetFile:
        with open(template) as tmp:
            raw = string.Template(tmp.read())
            targetFile.write(raw.substitute(**mapping))


if __name__ == "__main__":

    # input argument parser
    parser = argparse.ArgumentParser(
        description="Create a docker image for a given module and install script."
    )

    parser.add_argument("-m", "--modulepath", required=True, help="the path to the your module")
    parser.add_argument(
        "-i",
        "--installScript",
        required=True,
        help="Specify the installation script. In order to create portable containers, make sure"
        " that your installation script uses cmake options such as those in dumux/docker.opts (to"
        " avoid hardware-specific compiler flags).",
    )
    parser.add_argument(
        "-t", "--templateFolder", required=False, help="Specify the folder with the template files"
    )

    args = vars(parser.parse_args())

    # get information on the module
    modulePath = os.path.abspath(args["modulepath"])
    modInfo = extractModuleInfos(getModuleFile(modulePath), ["Module", "Maintainer"])
    moduleName = modInfo["Module"]
    moduleMaintainer = modInfo["Maintainer"]
    dockerTag = moduleName.lower()  # docker only supports lower case

    # get folder with the template files
    templateFolder = args["templateFolder"]
    if not templateFolder:
        templateFolder = os.path.join(os.path.dirname(__file__), "../docker")
        print(f"Folder with template files not given. Defaulting to '{templateFolder}'.")
    if not os.path.exists(templateFolder):
        sys.exit(f"Template folder {templateFolder} could not be found")

    print(
        "--> Ready to create the Dockerfile. In order to create portable containers, make sure that"
        " your installation script uses cmake options that do not set any hardware-specific"
        " compiler flags. For instance, you may use 'dumux/docker.opts'."
    )
    useCurrentCmakeOpts = queryYesNo("Do you understand the implications and want to continue?")
    if not useCurrentCmakeOpts:
        sys.exit("Abort.")

    print("*" * 54)
    print("\n-- Creating a Docker image for module " + moduleName + " --\n")
    print("*" * 54)

    if os.path.exists("docker"):
        delete = queryYesNo(
            "\nA docker folder already exists. Continue anyway? - will be overwritten"
        )
        if delete:
            shutil.rmtree("docker")
            print("--> Deleted old docker folder.")
        else:
            sys.exit("Abort.")

    os.mkdir("docker")
    print("--> Created the folder 'docker'.")

    # copy install script into docker folder and make it executable
    installScriptPath = args["installScript"]
    installScriptName = os.path.split(installScriptPath)[1]
    installScript = os.path.join(os.path.join(os.getcwd(), "docker"), installScriptName)
    shutil.copy(installScriptPath, installScript)
    os.system(f"chmod +x {installScript}")
    print(
        f"--> Using install script: {installScript} "
        f"to install dependencies for module {moduleName}."
    )

    # write setpermissions helper script
    substituteAndWrite(
        template=os.path.join(templateFolder, "setpermissions.sh.template"),
        target=os.path.join(os.getcwd(), "docker/setpermissions.sh"),
        mapping={},
    )
    print("--> Created permission helper script for easier container setup.")

    # write welcome message file
    substituteAndWrite(
        template=os.path.join(templateFolder, "WELCOME.template"),
        target=os.path.join(os.getcwd(), "docker/WELCOME"),
        mapping={"modName": moduleName, "modFolder": moduleName},
    )
    print("--> Created welcome message displayed on Docker container startup.")

    # write readme file
    substituteAndWrite(
        template=os.path.join(templateFolder, "README.md.template"),
        target=os.path.join(os.getcwd(), "docker/README.md"),
        mapping={"modName": moduleName, "dockerTag": dockerTag},
    )
    print("--> Created README.md on how to use the docker image.")

    # write helper file for container spin-up (make it executable after creation)
    dockerScript = os.path.join(os.getcwd(), f"docker/docker_{dockerTag}.sh")
    substituteAndWrite(
        template=os.path.join(templateFolder, "docker.sh.template"),
        target=dockerScript,
        mapping={"dockerTag": dockerTag},
    )
    os.system("chmod +x " + dockerScript)
    print("--> Created helper script to spin up the docker container.")

    # write the docker file
    substituteAndWrite(
        template=os.path.join(templateFolder, "Dockerfile.template"),
        target=os.path.join(os.getcwd(), "docker/Dockerfile"),
        mapping={
            "modName": moduleName,
            "modMaintainer": moduleMaintainer,
            "dockerTag": dockerTag,
            "instScript": installScriptName,
        },
    )
    print("--> Created Dockerfile. You can adapt it to your needs.")

    build = queryYesNo("Do you want to directly build the Docker image?")
    if build:
        print("Building Docker image... this may take several minutes.")
        try:
            os.chdir("docker")
            subprocess.run(
                ["docker", "build", "-f", "Dockerfile", "-t", dockerTag, "."], check=True
            )
            os.chdir("../")
        except subprocess.CalledProcessError:
            os.chdir("../")
            sys.exit("ERROR: docker image build failed")

        print(
            "",
            f"Successfully built image: {dockerTag}. ",
            "Have a look at docker/README.md.",
            "Check the container by running "
            f"'docker run -it {dockerTag} /bin/bash' in the same "
            "directory as the Dockerfile, and try using the convenience script "
            f"docker_{dockerTag}.sh",
            "See docker/README.md for more information.",
        )
    else:
        print(
            "You can build your Docker image later by running "
            f"'docker build -f Dockerfile -t {dockerTag}' ."
            "from within the folder 'docker' that was created by this script, "
            "and in which you should find the 'Dockerfile'."
        )
