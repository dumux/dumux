#!/usr/bin/env python3
import os
import sys
import shutil
import argparse
import subprocess
from getmoduleinfo import getModuleFile
from getmoduleinfo import extractModuleInfos

# require python 3
if sys.version_info[0] < 3:
    sys.exit("\nERROR: Python3 required")

# input argument parser
parser = argparse.ArgumentParser(description="This script creates a docker image for a given module and module installation script.")
parser.add_argument('-m', '--modulepath', required=True, help='the path to the your module')
parser.add_argument('-i', '--installScript', required=True, help="Specify the installation script")
args = vars(parser.parse_args());

# get the module name
modulePath = args['modulepath']
moduleFolder = os.path.relpath(modulePath, os.path.join(modulePath, '../'))
modInfo = extractModuleInfos(getModuleFile(moduleFolder), ['Module', 'Maintainer'])
moduleName = modInfo['Module']
moduleMaintainer = modInfo['Maintainer']

print("*"*54)
print("\n-- Creating a Docker image for module " + moduleName + " --\n" + "*"*54)

# docker only supports lower case so we convert the name
dockerTag = moduleName.lower()

if os.path.exists("docker"):
    print("\nA docker folder already exists. Continue anyway? - will be overwritten - [y/N]\n")
    delete = input()
    if delete == "y" or delete == "Y":
        shutil.rmtree("docker")
        print("--> Deleted old docker folder.")
    else:
        sys.exit("Abort.")

os.mkdir("docker")
print("--> Created the folder 'docker'.")

# copy install script into docker folder and make it executable
installScriptArg = args['installScript']
installScriptName = os.path.split(installScriptArg)[1]
installScript = os.path.join(os.path.join(os.getcwd(), 'docker'), installScriptName)
shutil.copy(installScriptArg, installScript)
os.system("chmod +x {}".format(installScript))
print("--> Using install script: {} to install dependencies for module {}.".format(installScript, moduleName))

# setpermissions helper script
setPermissionsFile = os.path.join(os.getcwd(), 'docker/setpermissions.sh')
with open(setPermissionsFile, 'w') as permissionFile:
    permissionFile.write('\n'.join(['#!/usr/bin/env bash',
                                    '',
                                    'echo "Setting permissions for shared folder"',
                                    '',
                                    '# if HOST_UID or HOST_GID are passed to the container',
                                    '# as environment variables, e.g. by calling',
                                    '# docker run -e HOST_UID=$(id -u $USER) -e HOST_GID=$(id -g $USER),',
                                    '# then we set the permissions of the files in the shared folder',
                                    'if [ "$HOST_UID" ]; then',
                                    '    echo "Changing user id to the provided one"',
                                    '    usermod -u $HOST_UID dumux',
                                    'fi',
                                    'if [ "$HOST_GID" ]; then',
                                    '    echo "Changing group id to the provided one"',
                                    '    groupmod -g $HOST_GID dumux',
                                    'fi',
                                    '',
                                    '# Change permissions only if both user and group id were passed.',
                                    '# Otherwise, this would change ownership to the default id of dumux,',
                                    '# which could lead to permission issues with the host user.',
                                    'if [ "${HOST_UID}" -a "${HOST_GID}" ]; then',
                                    '    # find all data in /dumux/shared/ and transfer ownership.',
                                    '    # sed "1d" removes the /dumux/shared folder itself (first line) that',
                                    '    # is within the results of the find command. If no files are present,',
                                    '    # chown returns an error because arguments are missing. Therefore, errors',
                                    '    # are redirected into /dev/null. Still, the script might return with an error',
                                    '    # in this case, and we guarantee successful execution with the || true trick at the end',
                                    '    find /dumux/shared/ | sed "1d" | xargs chown -R dumux:dumux 2> /dev/null || true',
                                    'else',
                                    '    echo "Skipping ownership transfer as host user and/or group id were not provided"',
                                    'fi']))
print("--> Created permission helper script for easier container setup.")

# welcome message
welcomeFile = os.path.join(os.getcwd(), 'docker/WELCOME')
with open(welcomeFile, 'w') as welcome:
    welcome.write('\n'.join(['',
                             '*'*50,
                             'Welcome to the dumux-pub Docker container "{}"!'.format(moduleName),
                             '',
                             'You can run pre-compiled examples in the sub-folders of PATH_TO_MODULES/{}/build-cmake/,'.format(moduleFolder),
                             'or compile them there using make <example>.',
                             '',
                             'Note that the container does not have graphics support, but copying the output (e.g. VTK files) to the folder',
                             '/dumux/shared will make them available outside this container on the host for display.\n']))
print("--> Created welcome message displayed on Docker container startup.")

# readme
readmeFile = os.path.join(os.getcwd(), 'docker/README.md')
with open(readmeFile, 'w') as readme:
    readme.write('\n'.join(['# Readme for the dumux pub table {}'.format(moduleName),
                            '',
                            'You created the Docker image {}. Next steps:'.format(dockerTag),
                            '',
                            '* Try your container by running pubtable_{}.sh open'.format(dockerTag),
                            '  See below for instructions how to share files with the host system.',
                            '* Push the docker image to DockerHub or the GitLab Docker registry of your dumux-pub module.',
                            '  Look at the Registry tab of your dumux-pub module for help.',
                            '',
                            '* Replace the image name in pubtable_{}.sh with the actual image name'.format(dockerTag),
                            '  e.g. git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a.',
                            '',
                            '* [Optional] Add the Dockerfile to the git repository.',
                            '',
                            '* Add the pubtable_{}.sh script to the git repository (this is for the user)'.format(dockerTag),
                            '  and add the following lines to your README.md:',
                            '',
                            'Using the pub table {} with docker'.format(moduleName),
                            '=============================================',
                            '',
                            'In order to run simulations of this pub table look',
                            'at the convenience script pubtable_{}.sh.'.format(dockerTag),
                            'First download the script from the git repository.',
                            '',
                            'The simplest way is to spin up a container',
                            'is creating a new folder "dumux"',
                            '$ mkdir dumux',
                            'change to the new folder',
                            '$ cd dumux',
                            'and open the pub table by running',
                            '$ pubtable_{}.sh open'.format(dockerTag),
                            '',
                            'The container will spin up. It will mount the "dumux"',
                            'directory into the container at /dumux/shared. Put files',
                            'in this folder to share them with the host machine. This',
                            'could be e.g. VTK files produced by the simulation that',
                            'you want to visualize on the host machine.']))
print("--> Created README.md on how to use the docker image.")

# pub table file (make it executable)
pubTableFile = os.path.join(os.getcwd(), 'docker/pubtable_{}.sh'.format(dockerTag))
with open(pubTableFile, "w") as pubTable:
    pubTable.write('\n'.join(['#!/usr/bin/env bash',
                              '',
                              '# TODO set the image name here to a global address',
                              '#      after adding your container to your modules container registry',
                              '# e.g. git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a',
                              'PUB_IMAGE_NAME={}'.format(dockerTag),
                              '',
                              '# the host directory ...',
                              'SHARED_DIR_HOST="$(pwd)"',
                              '# ... that is mounted into this container directory:',
                              'SHARED_DIR_CONTAINER="/dumux/shared"',
                              '',
                              'help ()',
                              '{',
                              '    echo ""',
                              '    echo "Usage: pubtable_{}.sh <command> [options]"'.format(dockerTag),
                              '    echo ""',
                              '    echo "  pubtable_{}.sh open [image]      - open a pub table."'.format(dockerTag),
                              '    echo "  pubtable_{}.sh help              - display this message."'.format(dockerTag),
                              '    echo ""',
                              '    echo "Optionally supply a Docker image name to the open command."',
                              '    echo ""',
                              '}',
                              '',
                              '# open a pub table. Only argument is the Docker image.',
                              'open()',
                              '{',
                              '    IMAGE="$1"',
                              '    docker run -it \\',
                              '                -e HOST_UID=$(id -u $USER) \\',
                              '                -e HOST_GID=$(id -g $USER) \\',
                              '                -v $SHARED_DIR_HOST:$SHARED_DIR_CONTAINER \\',
                              '                --name dumuxpub_{} \\'.format(dockerTag),
                              '                $IMAGE /bin/bash',
                              '}',
                              '',
                              '# Check if user specified valid command otherwise print help message',
                              'if [ "$1" == "open" ]; then',
                              '    IMAGE="$2" : ${IMAGE:="$PUB_IMAGE_NAME"}',
                              '    open $IMAGE',
                              'else',
                              '    help',
                              '    exit 1',
                              'fi']))
os.system("chmod +x " + pubTableFile)
print("--> Created " + pubTableFile + " script to spin up the docker container.")

# write the docker file
with open("docker/Dockerfile", "w") as dockerFile:
    dockerFile.write('\n'.join(['# ' + moduleName + ' docker container',
                                '# see https://github.com/phusion/baseimage-docker for information on the base image',
                                '# It is Ubuntu LTS customized for better Docker compatibility',
                                'FROM phusion/baseimage:18.04-1.0.0',
                                'MAINTAINER ' + moduleMaintainer + '',
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
                                '    python3-dev \\',
                                '    python3-pip \\',
                                '    git \\',
                                '    pkg-config \\',
                                '    build-essential \\',
                                '    gfortran \\',
                                '    mpi-default-bin \\',
                                '    mpi-default-dev \\',
                                '    libsuitesparse-dev \\',
                                '    libsuperlu-dev \\',
                                '    libeigen3-dev \\',
                                '    doxygen \\',
                                '    wget \\',
                                '    && apt-get clean \\',
                                '    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*',
                                '',
                                '# get cmake 3.18 and export it to path',
                                'RUN mkdir -p cmake-3.18',
                                'RUN wget -qO- "https://cmake.org/files/v3.18/cmake-3.18.0-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C cmake-3.18',
                                'ENV PATH="/cmake-3.18/bin:${PATH}"',
                                '',
                                '# add the permission helper script to the my_init service',
                                'COPY setpermissions.sh /etc/my_init.d/setpermissions.sh',
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
                                '# switch to the dumux user and set the working directory',
                                'USER dumux',
                                'WORKDIR /dumux',
                                '',
                                '# create a shared volume communicating with the host',
                                'RUN mkdir /dumux/shared',
                                'VOLUME /dumux/shared',
                                '',
                                '# This is the message printed on entry',
                                'COPY WELCOME /dumux/WELCOME',
                                '',
                                '# Install the dumux-pub module and its dependencies',
                                '# This expects the install script to do everything from clone to configure',
                                'COPY ' + installScriptName + ' /dumux/' + installScriptName + '',
                                'RUN ./' + installScriptName + ' && rm -f /dumux/' + installScriptName + '',
                                '',
                                '# build doxygen documentation',
                                'WORKDIR /dumux',
                                'RUN cd ' + moduleFolder + '/build-cmake && make doc',
                                '',
                                '# switch back to root',
                                'USER root',
                                '',
                                '# set entry point like advised https://github.com/phusion/baseimage-docker',
                                '# this sets the permissions right, see above',
                                'ENTRYPOINT ["/sbin/my_init","--quiet","--","/sbin/setuser","dumux","/bin/bash","-l","-c"]',
                                '',
                                '# start interactive shell',
                                'CMD ["/bin/bash","-i"]']))
print("--> Created Dockerfile. You can adapt it to your needs.")
print()
print("Do you want to directly build the Docker image? [y/N]")

build = input()
if build == "y" or build == "Y":
    print("Building Docker image... this may take several minutes.")
    try:
        os.chdir('docker')
        subprocess.run(['docker', 'build', '-f', 'Dockerfile', '-t', dockerTag, '.'], check=True)
        os.chdir('../')
    except:
        os.chdir('../')
        sys.exit("ERROR: docker image build failed")
    print()
    print("Successfully built docker image: " + dockerTag + ". Have a look at docker/README.md.")
    print("Check the container running docker run -it " + dockerTag + " /bin/bash in the same")
    print("directory as the Dockerfile and try using the convenience script pubtable_" + dockerTag + ".sh")
    print("See docker/README.md for more information.")
else:
    print("You can build your Docker image later by running docker build -f Dockerfile -t " + dockerTag)
    print("from within the folder 'docker' that was created by this script, and where the Dockerfile is.")
