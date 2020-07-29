#!/usr/bin/env python3

import os
import sys
import shutil

def show_message(message):
    print("=" * 54)
    print(message)
    print("=" * 54)

def print_help():
    pass

# get the module name
dune_module_file = "dune.module"
#dune_module_file = "dune.module"
if not os.path.exists(dune_module_file):
    sys.exit("\n"
             "ERROR: Could not find dune.module.\n"
             "Be sure that you are inside your module!")

file_object = open(dune_module_file, "r")
for line in file_object:
    if "Module:" in line:
        module_name = line.split()[1]
    if "Maintainer:" in line:
        module_maintainer = line.split()[1]
# docker only support lower case so we convert the name
docker_tag = module_name.lower()

show_message("Creating a Docker image from module {}".format(module_name))

# if install_dep_script was not specified use the default
#install_dep_script = "install{}.sh".format(module_name)
install_dep_script = "install.sh"
print("--> Using install script: {} to install dependencies for module {}."
      .format(install_dep_script, module_name))

# make it executable
os.system("chmod +x {}".format(install_dep_script))

# remove docker folder if exist
if os.path.exists("docker"):
    print("\nA docker folder already exists. Continue anyway?"
          "- will be overwritten - [y/N]\n")
    delete = input()
    if delete == "y" or delete == "Y":
        shutil.rmtree("docker")
        print("--> Deleted old docker folder.")
    else:
        sys.exit("Abort.")

# dockerignore file
dockerignore = open(".dockerignore", "w+")
dockerignore.writelines(["build*\n", "CMakeFiles"])
dockerignore.close()
print("--> Created dockerignore file.")

# make docker folder
os.mkdir("docker")

# also copy it here in case we build here
shutil.copyfile(".dockerignore","docker/.dockerignore")

# setpermissions helper script
setpermissions_content = [
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
setpermissions.write('\n'.join(setpermissions_content) + "\n")
setpermissions.close()
print("--> Created permission helper script for easier container setup.")

welcome_content = [
    "Welcome to the dumux-pub Docker container {}. You can run pre-compiled examples".format(module_name),
    "in the folder {}/build-cmake/appl, or compile them there using make <example>.".format(module_name),
    "The container doesn't have graphics support.",
    "Copying the output like VTK files to the folder /dumux/shared will make them available",
    "outside this container on the host for display."]
welcome = open("docker/WELCOME", "x")
welcome.write('\n'.join(welcome_content) + "\n")
welcome.close()
print("--> Created welcome message displayed on Docker container startup.")

readme_content = [
    '# readme for the dumux pub table {}'.format(module_name),
    '',
    'You created a Docker image {}. Next steps:'.format(docker_tag),
    '',
    '* Try your container by running pubtable_{} open'.format(docker_tag),
    '  See below for instructions how to share files with the host system.',
    '',
    '* Push the docker image to DockerHub or the GitLab Docker registry of your dumux-pub module.',
    '  Look at the Registry tab of your dumux-pub module for help.',
    '',
    '* Replace the image name in pubtable_{} with the actual image name'.format(docker_tag),
    '  e.g. git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a.',
    '',
    '* [Optional] Add the Dockerfile to the git repository.',
    '',
    '* Add the pubtable_{} script to the git repository (this is for the user)'.format(docker_tag),
    '  and add the following lines to your README.md:',
    '',
    'Using the pub table {} with docker'.format(module_name),
    '=============================================',
    '',
    'In order to run simulations of this pub table look',
    'at the convenience script pubtable_{}.'.format(docker_tag),
    'First download the script from the git repository.',
    '',
    'The simplest way is to spin up a container',
    'is creating a new folder "dumux"',
    '$ mkdir dumux',
    'change to the new folder',
    '$ cd dumux',
    'and open the pub table by running',
    '$ pubtable_{} open'.format(docker_tag),
    '',
    'The container will spin up. It will mount the "dumux"',
    'directory into the container at /dumux/shared. Put files',
    'in this folder to share them with the host machine. This',
    'could be e.g. VTK files produced by the simulation that',
    'you want to visualize on the host machine.',
    '']
readme = open("docker/README.md", "x")
readme.write('\n'.join(readme_content) + "\n")
readme.close()
print("--> Created README.md on how to use the docker image.")

putable_content = [
    '#!/usr/bin/env bash',
    '',
    '# TODO set the image name here to a global address',
    "# e.g. 'git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a'",
    'PUB_IMAGE_NAME={}'.format(docker_tag),
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
    '    echo "Usage: pubtable_{} <command> [options]"'.format(docker_tag),
    '    echo ""',
    '    echo "  pubtable_{} open [image]      - open a pub table."'.format(docker_tag),
    '    echo "  pubtable_{} help              - display this message."'.format(docker_tag),
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
    '             --name dumuxpub_ \\'.format(docker_tag),
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
putable = open("docker/putable_{}".format(docker_tag), "x")
putable.write('\n'.join(putable_content) + "\n")
putable.close()
os.system("chmod +x docker/putable_{}".format(docker_tag))
print("--> Created pubtable_{} script to spin up the docker container.".format(docker_tag))

dockerfile_content = [
    '# {} docker container'.format(module_name),
    '# see https://github.com/phusion/baseimage-docker for information on the base image',
    '# It is Ubuntu LTS customized for better Docker compatibility',
    'FROM phusion/baseimage:0.9.22',
    'MAINTAINER {}'.format(module_maintainer),
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
    'COPY . /dumux/{}'.format(module_name),
    'RUN chown -R dumux:dumux /dumux/{}'.format(module_name),
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
    'COPY {} /dumux/{}'.format(install_dep_script, install_dep_script),
    'RUN ./{} && rm -f /dumux/{}'.format(install_dep_script, install_dep_script),
    '',
    '# configure module',
    'RUN /dumux/dune-common/bin/dunecontrol --opts=/dumux/dumux/optim.opts all',
    '',
    '# build doxygen documentation and tests',
    '# all applications that use dune_add_test will be built like this',
    'RUN cd {}/build-cmake && make doc && make -j4 build_tests'.format(module_name),
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
dockerfile.write('\n'.join(dockerfile_content) + "\n")
dockerfile.close()

print("--> Created Dockerfile. You can adapt it to your needs, or\n"
      "build the Docker image directly? [y/N]")
build = input()
if build == "y" or build == "Y":
    print("Building Docker image... this may take several minutes.")
    os.system("docker build -f docker/Dockerfile -t {} .".format(docker_tag))
    print("\n"
          + "Successfully built docker image: {}. Have a look at docker/README.md.\n".format(docker_tag)
          + "Check the container running docker run -it {} /bin/bash \n".format(docker_tag)
          + "in the same directory as the Dockerfile.\n"
          + "And try using the convenience script pubtable_{}, see docker/README.md.".format(docker_tag))
else:
    print("You can build your Docker image later by running docker build -f docker/Dockerfile -t {} .\n".format(docker_tag)
          + "in your module directory, i.e. above the docker folder containing Dockerfile.")
