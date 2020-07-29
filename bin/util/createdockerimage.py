#!/usr/bin/env python3
import os
import sys
import shutil
import argparse
from argparse import RawTextHelpFormatter

# require python 3
if sys.version_info[0] < 3:
    sys.exit("\nERROR: Python3 required")

# input argument parser
parser = argparse.ArgumentParser(
        description="This is to be run inside an extracted module, e.g. created by the extractmodulepart script.\n\n" + \
                    "USAGE: python3 " + sys.argv[0] + " -i installDepScript\n\n" + \
                    "installDepScript defaults to the file install{moduleName}.sh." + \
                    "You can also manually specify a script installing the DUNE dependencies for your module\n",
        formatter_class=RawTextHelpFormatter
)

parser.add_argument('-i', '--installScript', help="Specify the script that installs dependencies")
args = vars(parser.parse_args());

# prints a status message with lines of '=' characters above and below
def showMessage(message):
    print("=" * 54)
    print(message)
    print("=" * 54)

# get the module name
if not os.path.exists("dune.module"):
    sys.exit("\n"
             "ERROR: Could not find dune.module.\n"
             "Make sure that you are inside your module!")

for line in open("dune.module", 'r').readlines():
    if "Module:" in line:
        moduleName = line.split()[1]
    if "Maintainer:" in line:
        moduleMaintainer = line.split()[1]

# docker only supports lower case so we convert the name
dockerTag = moduleName.lower()
showMessage("Creating a Docker image from module {}".format(moduleName))

# if installScript was not specified use the default
installDepScript = args['installScript']
if installDepScript == None: installDepScript = "install" + moduleName + ".sh"
print("--> Using install script: {} to install dependencies for module {}.".format(installDepScript, moduleName))

# make it executable
os.system("chmod +x {}".format(installDepScript))

# remove docker folder if exist
if os.path.exists("docker"):
    print("\nA docker folder already exists. Continue anyway? - will be overwritten - [y/N]\n")
    delete = input()
    if delete == "y" or delete == "Y":
        shutil.rmtree("docker")
        print("--> Deleted old docker folder.")
    else:
        sys.exit("Abort.")

# dockerIgnore file
dockerIgnore = open(".dockerIgnore", "w+")
dockerIgnore.writelines(["build*\n", "CMakeFiles"])
dockerIgnore.close()
print("--> Created dockerIgnore file.")

# make docker folder
os.mkdir("docker")

# also copy it here in case we build here
shutil.copyfile(".dockerIgnore", "docker/.dockerIgnore")

# setpermissions helper script
setPermissionsContent = [
    '#!/bin/bash',
    '# The user can pass the user and group id by passing',
    '# --env HOST_UID=$(id -u $USER) --env HOST_GID=$(id -g $USER)',
    '# with the UID on the host to the container. Should work for Linux, Mac and Windows.',
    '# Allows to manage the writes for shared volumes.',
    'if [ "$HOST_UID" ]; then',
    '    usermod -u $HOST_UID dumux',
    'fi',
    'if [ "$HOST_GID" ]; then',
    '    groupmod -g $HOST_GID dumux',
    'fi',
    '# Make sure that everything in /dumux is accessible by the dumux user',
    '# sed "1d" removes forst line which is the folder itself',
    '# exclude the shared folder using grep. This folder should be accessible by setting the UID/GID.',
    '# chown transfers ownership to group dumux user dumux',
    '# the OR true trick avoids the script exiting due to use of set -e',
    'find /dumux -maxdepth 1 | sed "1d" | grep -v "/dumux/shared" | xargs chown -R dumux:dumux 2> /dev/null || true',
    '']
setpermissions = open("docker/setpermissions.sh", "x")
setpermissions.write('\n'.join(setPermissionsContent) + "\n")
setpermissions.close()
print("--> Created permission helper script for easier container setup.")

welcomeContent = [
    "Welcome to the dumux-pub Docker container {}. You can run pre-compiled examples".format(moduleName),
    "in the folder {}/build-cmake/appl, or compile them there using make <example>.".format(moduleName),
    "The container doesn't have graphics support.",
    "Copying the output like VTK files to the folder /dumux/shared will make them available",
    "outside this container on the host for display."]
welcome = open("docker/WELCOME", "x")
welcome.write('\n'.join(welcomeContent) + "\n")
welcome.close()
print("--> Created welcome message displayed on Docker container startup.")

readmeContent = [
    '# readme for the dumux pub table {}'.format(moduleName),
    '',
    'You created a Docker image {}. Next steps:'.format(dockerTag),
    '',
    '* Try your container by running pubtable_{} open'.format(dockerTag),
    '  See below for instructions how to share files with the host system.',
    '',
    '* Push the docker image to DockerHub or the GitLab Docker registry of your dumux-pub module.',
    '  Look at the Registry tab of your dumux-pub module for help.',
    '',
    '* Replace the image name in pubtable_{} with the actual image name'.format(dockerTag),
    '  e.g. git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a.',
    '',
    '* [Optional] Add the Dockerfile to the git repository.',
    '',
    '* Add the pubtable_{} script to the git repository (this is for the user)'.format(dockerTag),
    '  and add the following lines to your README.md:',
    '',
    'Using the pub table {} with docker'.format(moduleName),
    '=============================================',
    '',
    'In order to run simulations of this pub table look',
    'at the convenience script pubtable_{}.'.format(dockerTag),
    'First download the script from the git repository.',
    '',
    'The simplest way is to spin up a container',
    'is creating a new folder "dumux"',
    '$ mkdir dumux',
    'change to the new folder',
    '$ cd dumux',
    'and open the pub table by running',
    '$ pubtable_{} open'.format(dockerTag),
    '',
    'The container will spin up. It will mount the "dumux"',
    'directory into the container at /dumux/shared. Put files',
    'in this folder to share them with the host machine. This',
    'could be e.g. VTK files produced by the simulation that',
    'you want to visualize on the host machine.',
    '']
readme = open("docker/README.md", "x")
readme.write('\n'.join(readmeContent) + "\n")
readme.close()
print("--> Created README.md on how to use the docker image.")

pubTableContent = [
    '#!/usr/bin/env bash',
    '',
    '# TODO set the image name here to a global address',
    "# e.g. 'git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a'",
    'PUB_IMAGE_NAME={}'.format(dockerTag),
    '',
    '# the host directory that is mounted into...',
    'SHARED_DIR_HOST="$(pwd)"',
    '# ...this container directory.',
    'SHARED_DIR_CONTAINER="/dumux/shared"',
    '',
    'executeandlog ()',
    '{',
    '    echo "[$@]"',
    '    eval $@',
    '}',
    '',
    'help ()',
    '{',
    '    echo ""',
    '    echo "Usage: pubtable_{} <command> [options]"'.format(dockerTag),
    '    echo ""',
    '    echo "  pubtable_{} open [image]      - open a pub table."'.format(dockerTag),
    '    echo "  pubtable_{} help              - display this message."'.format(dockerTag),
    '    echo ""',
    '    echo "Optionally supply an Docker image name to the open command."',
    '    echo ""',
    '}',
    '',
    '# open a pub table. Only argument is the Docker image.',
    'open()',
    '{',
    '    IMAGE="$1"',
    '    COMMAND="docker create -ti \\',
    '             -e HOST_UID=$(id -u \$USER) \\',
    '             -e HOST_GID=$(id -g \$USER) \\',
    '             -v $SHARED_DIR_HOST:\$SHARED_DIR_CONTAINER \\',
    '             --name dumuxpub_ \\'.format(dockerTag),
    '             $IMAGE /bin/bash"',
    ' ',
    '    CONTAINER_NAME=\$(executeandlog \$COMMAND | tail -n 1)',
    '    executeandlog docker start -a -i \$CONTAINER_NAME',
    '    executeandlog docker rm \$CONTAINER_NAME',
    '}',
    '',
    '# Check if user specified valid command otherwise print help message',
    'if [ "$1" == "open" ]; then',
    '    IMAGE="$2" : \${IMAGE:="$PUB_IMAGE_NAME"}',
    '    open $IMAGE',
    'else',
    '    help',
    '    exit 1',
    'fi',
    '']
pubTable = open("docker/pubTable_{}".format(dockerTag), "x")
pubTable.write('\n'.join(pubTableContent) + "\n")
pubTable.close()
os.system("chmod +x docker/pubTable_{}".format(dockerTag))
print("--> Created pubtable_{} script to spin up the docker container.".format(dockerTag))

dockerFileContent = [
    '# {} docker container'.format(moduleName),
    '# see https://github.com/phusion/baseimage-docker for information on the base image',
    '# It is Ubuntu LTS customized for better Docker compatibility',
    'FROM phusion/baseimage:0.9.22',
    'MAINTAINER {}'.format(moduleMaintainer),
    '',
    '# run Ubuntu update as advised on https://github.com/phusion/baseimage-docker',
    'RUN apt-get update \\',
    '    && apt-get upgrade -y -o Dpkg::Options::="--force-confold" \\',
    '    && apt-get clean \\',
    '    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*',
    '',
    '# install the basic dependencies',
    'RUN apt-get update \\',
    '    && apt-get install --no-install-recommends --yes \\',
    '    ca-certificates \\',
    '    vim \\',
    '    python-dev \\',
    '    python-pip \\',
    '    git \\',
    '    pkg-config \\',
    '    build-essential \\',
    '    gfortran \\',
    '    cmake \\',
    '    mpi-default-bin \\',
    '    mpi-default-dev \\',
    '    libsuitesparse-dev \\',
    '    libsuperlu-dev \\',
    '    libeigen3-dev \\',
    '    doxygen \\',
    '    && apt-get clean \\',
    '    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*',
    '',
    '# add the permission helper script to the my_init service',
    'COPY ./docker/setpermissions.sh /etc/my_init.d/setpermissions.sh',
    '',
    '# create a dumux user',
    '# add the welcome message (copied further down) output to bashrc',
    '# make the set permission helper script executable',
    '# add user to video group which enables graphics if desired',
    'RUN useradd -m --home-dir /dumux dumux \\',
    '    && echo "cat /dumux/WELCOME" >> /dumux/.bashrc \\',
    '    && chmod +x /etc/my_init.d/setpermissions.sh \\',
    '    && usermod -a -G video dumux',
    '',
    '# Turn off verbose syslog messages as described here:',
    '# https://github.com/phusion/baseimage-docker/issues/186',
    'RUN touch /etc/service/syslog-forwarder/down',
    '',
    '# copy the extracted dumux-pub module and make dumux own it',
    'COPY . /dumux/{}'.format(moduleName),
    'RUN chown -R dumux:dumux /dumux/{}'.format(moduleName),
    '',
    '# switch to the dumux user and set the working directory',
    'USER dumux',
    'WORKDIR /dumux',
    '',
    '# create a shared volume communicating with the host',
    'RUN mkdir /dumux/shared',
    'VOLUME /dumux/shared',
    '',
    '# This is the message printed on entry',
    'COPY ./docker/WELCOME /dumux/WELCOME',
    '',
    '# install dumux-pub module dependencies',
    'COPY {} /dumux/{}'.format(installDepScript, installDepScript),
    'RUN ./{} && rm -f /dumux/{}'.format(installDepScript, installDepScript),
    '',
    '# configure module',
    'RUN /dumux/dune-common/bin/dunecontrol --opts=/dumux/dumux/optim.opts all',
    '',
    '# build doxygen documentation and tests',
    '# all applications that use dune_add_test will be built like this',
    'RUN cd {}/build-cmake && make doc && make -j4 build_tests'.format(moduleName),
    '',
    '# switch back to root',
    'USER root',
    '',
    '# set entry point like advised https://github.com/phusion/baseimage-docker',
    '# this sets the permissions right, see above',
    'ENTRYPOINT ["/sbin/my_init","--quiet","--","/sbin/setuser","dumux","/bin/bash","-l","-c"]',
    '',
    '# start interactive shell',
    'CMD ["/bin/bash","-i"]',
    '']
dockerfile = open("docker/Dockerfile", "x")
dockerfile.write('\n'.join(dockerFileContent) + "\n")
dockerfile.close()

print("--> Created Dockerfile. You can adapt it to your needs, or\n"
      "build the Docker image directly? [y/N]")
build = input()
if build == "y" or build == "Y":
    print("Building Docker image... this may take several minutes.")
    os.system("docker build -f docker/Dockerfile -t {} .".format(dockerTag))
    print("\n"
          + "Successfully built docker image: {}. Have a look at docker/README.md.\n".format(dockerTag)
          + "Check the container running docker run -it {} /bin/bash \n".format(dockerTag)
          + "in the same directory as the Dockerfile.\n"
          + "And try using the convenience script pubtable_{}, see docker/README.md.".format(dockerTag))
else:
    print("You can build your Docker image later by running docker build -f docker/Dockerfile -t {} .\n".format(dockerTag)
          + "in your module directory, i.e. above the docker folder containing Dockerfile.")
