#! /bin/bash

# check if help is needed
if test "$1" = "--help" || test "$1" = "-help" \
   || test "$1" = "help" || test "$1" = "-h"; then
  echo ""
  echo "This is to be run inside an extracted module, e.g. created by the extractmodulepart script."
  echo ""
  echo "USAGE: $0 [INSTALL_DEP_SCRIPT ...]"
  echo ""
  echo "INSTALL_DEP_SCRIPT defaults to the file install{MODULE_NAME}.sh. You can also manually "
  echo "specify a script installing the DUNE dependencies for your module"
  echo ""
  exit 0
fi

# get the module name
DUNE_MODULE_FILE=$(find . -name 'dune.module')
if test "$DUNE_MODULE_FILE" = ""; then
  echo "ERROR: Could not find dune.module."
  echo "Be sure that you are inside your module!"
  exit 1;
fi

MODULE_NAME=$(awk '/Module/{print $NF}' $DUNE_MODULE_FILE)
MODULE_MAINTAINER=$(awk '/Maintainer/{print $NF}' $DUNE_MODULE_FILE)
# docker only support lower case so we convert the name
DOCKER_TAG=$(echo $MODULE_NAME | tr '[:upper:]' '[:lower:]')

echo "======================================================"
echo "Creating a Docker image from module $MODULE_NAME"
echo "======================================================"

# if INSTALL_DEP_SCRIPT was not specified use the default
INSTALL_DEP_SCRIPT=${1:-install$MODULE_NAME.sh}
echo "--> Using install script: $INSTALL_DEP_SCRIPT to install dependencies for module $MODULE_NAME."

# make it executable
chmod +x $INSTALL_DEP_SCRIPT

# make a docker folder
if [ -d docker ]; then
    echo ""
    read -p "   A docker folder already exists. Continue anyway? - will be overwritten - [y/N] " DELETE
    echo ""
    if test x$DELETE = xy -o x$DELETE = xY; then
        rm -r docker
        echo "--> Deleted old docker folder."
    else
        echo "Abort."
        exit 1
    fi
fi

# dockerignore file
if [ -e .dockerignore ]; then
    rm .dockerignore
fi
touch .dockerignore
echo "build*" >> .dockerignore
echo "CMakeFiles" >> .dockerignore
echo "--> Created dockerignore file."

mkdir docker
cd docker

# also copy it here in case we build here
cp ../.dockerignore .

# setpermissions helper script
touch setpermissions.sh
cat <<EOT >> setpermissions.sh
#!/bin/bash
# The user can pass the user and group id by passing
# --env HOST_UID=\$(id -u \$USER) --env HOST_GID=\$(id -g \$USER)
# with the UID on the host to the container. Should work for Linux, Mac and Windows.
# Allows to manage the writes for shared volumes.
if [ "\$HOST_UID" ]; then
    usermod -u \$HOST_UID dumux
fi
if [ "\$HOST_GID" ]; then
    groupmod -g \$HOST_GID dumux
fi
# Make sure that everything in /dumux is accessible by the dumux user
# sed "1d" removes forst line which is the folder itself
# exclude the shared folder using grep. This folder should be accessible by setting the UID/GID.
# chown transfers ownership to group dumux user dumux
# the OR true trick avoids the script exiting due to use of set -e
find /dumux -maxdepth 1 | sed "1d" | grep -v "/dumux/shared" | xargs chown -R dumux:dumux 2> /dev/null || true

EOT

echo "--> Created permission helper script for easier container setup."

touch WELCOME
cat <<EOT >> WELCOME
Welcome to the dumux-pub Docker container $MODULE_NAME. You can run pre-compiled examples
in the folder $MODULE_NAME/build-cmake/appl, or compile them there using make <example>.
The container doesn't have graphics support.
Copying the output like VTK files to the folder /dumux/shared will make them available
outside this container on the host for display.
EOT

echo "--> Created welcome message displayed on Docker container startup."

touch README.md
cat <<EOT >> README.md
# readme for the dumux pub table $MODULE_NAME

You created a Docker image $DOCKER_TAG. Next steps:

* Try your container by running pubtable_$DOCKER_TAG open
  See below for instructions how to share files with the host system.

* Push the docker image to DockerHub or the GitLab Docker registry of your dumux-pub module.
  Look at the Registry tab of your dumux-pub module for help.

* Replace the image name in pubtable_$DOCKER_TAG with the actual image name
  e.g. git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a.

* [Optional] Add the Dockerfile to the git repository.

* Add the pubtable_$DOCKER_TAG script to the git repository (this is for the user)
  and add the following lines to your README.md:

Using the pub table $MODULE_NAME with docker
=============================================

In order to run simulations of this pub table look
at the convenience script pubtable_$DOCKER_TAG.
First download the script from the git repository.

The simplest way is to spin up a container
is creating a new folder "dumux"
$ mkdir dumux
change to the new folder
$ cd dumux
and open the pub table by running
$ pubtable_$DOCKER_TAG open

The container will spin up. It will mount the "dumux"
directory into the container at /dumux/shared. Put files
in this folder to share them with the host machine. This
could be e.g. VTK files produced by the simulation that
you want to visualize on the host machine.

EOT

echo "--> Created README.md on how to use the docker image."

touch pubtable_$DOCKER_TAG
chmod +x pubtable_$DOCKER_TAG
cat <<EOT >> pubtable_$DOCKER_TAG
#!/usr/bin/env bash

# TODO set the image name here to a global address
# e.g. 'git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a'
PUB_IMAGE_NAME=$DOCKER_TAG

# the host directory that is mounted into...
SHARED_DIR_HOST="\$(pwd)"
# ...this container directory.
SHARED_DIR_CONTAINER="/dumux/shared"

executeandlog ()
{
    echo "[\$@]"
    eval \$@
}

help ()
{
    echo ""
    echo "Usage: pubtable_$DOCKER_TAG <command> [options]"
    echo ""
    echo "  pubtable_$DOCKER_TAG open [image]      - open a pub table."
    echo "  pubtable_$DOCKER_TAG help              - display this message."
    echo ""
    echo "Optionally supply an Docker image name to the open command."
    echo ""
}

# open a pub table. Only argument is the Docker image.
open()
{
    IMAGE="\$1"
    COMMAND="docker create -ti \\
             -e HOST_UID=\$(id -u \$USER) \\
             -e HOST_GID=\$(id -g \$USER) \\
             -v \$SHARED_DIR_HOST:\$SHARED_DIR_CONTAINER \\
             --name dumuxpub_$DOCKER_TAG \\
             \$IMAGE /bin/bash"

    CONTAINER_NAME=\$(executeandlog \$COMMAND | tail -n 1)
    executeandlog docker start -a -i \$CONTAINER_NAME
    executeandlog docker rm \$CONTAINER_NAME
}

# Check if user specified valid command otherwise print help message
if [ "\$1" == "open" ]; then
    IMAGE="\$2" : \${IMAGE:="\$PUB_IMAGE_NAME"}
    open \$IMAGE
else
    help
    exit 1
fi

EOT

echo "--> Created pubtable_$DOCKER_TAG script to spin up the docker container."


touch Dockerfile
cat <<EOT >> Dockerfile
# $MODULE_NAME docker container
# see https://github.com/phusion/baseimage-docker for information on the base image
# It is Ubuntu LTS customized for better Docker compatibility
FROM phusion/baseimage:18.04-1.0.0
MAINTAINER $MODULE_MAINTAINER

# run Ubuntu update as advised on https://github.com/phusion/baseimage-docker
RUN apt-get update \\
    && apt-get upgrade -y -o Dpkg::Options::="--force-confold" \\
    && apt-get clean \\
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install the basic dependencies
RUN apt-get update \\
    && apt-get install --no-install-recommends --yes \\
    ca-certificates \\
    vim \\
    python-dev \\
    python-pip \\
    git \\
    pkg-config \\
    build-essential \\
    gfortran \\
    cmake \\
    mpi-default-bin \\
    mpi-default-dev \\
    libsuitesparse-dev \\
    libsuperlu-dev \\
    libeigen3-dev \\
    doxygen \\
    && apt-get clean \\
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# add the permission helper script to the my_init service
COPY ./docker/setpermissions.sh /etc/my_init.d/setpermissions.sh

# create a dumux user
# add the welcome message (copied further down) output to bashrc
# make the set permission helper script executable
# add user to video group which enables graphics if desired
RUN useradd -m --home-dir /dumux dumux \\
    && echo "cat /dumux/WELCOME" >> /dumux/.bashrc \\
    && chmod +x /etc/my_init.d/setpermissions.sh \\
    && usermod -a -G video dumux

# copy the extracted dumux-pub module and make dumux own it
COPY . /dumux/$MODULE_NAME
RUN chown -R dumux:dumux /dumux/$MODULE_NAME

# switch to the dumux user and set the working directory
USER dumux
WORKDIR /dumux

# create a shared volume communicating with the host
RUN mkdir /dumux/shared
VOLUME /dumux/shared

# This is the message printed on entry
COPY ./docker/WELCOME /dumux/WELCOME

# install dumux-pub module dependencies
COPY $INSTALL_DEP_SCRIPT /dumux/$INSTALL_DEP_SCRIPT
RUN ./$INSTALL_DEP_SCRIPT && rm -f /dumux/$INSTALL_DEP_SCRIPT

# configure module
RUN /dumux/dune-common/bin/dunecontrol --opts=/dumux/dumux/cmake.opts all

# build doxygen documentation and tests
# all applications that use dune_add_test will be built like this
RUN cd $MODULE_NAME/build-cmake && make doc && make -j4 build_tests

# switch back to root
USER root

# set entry point like advised https://github.com/phusion/baseimage-docker
# this sets the permissions right, see above
ENTRYPOINT ["/sbin/my_init","--quiet","--","/sbin/setuser","dumux","/bin/bash","-l","-c"]

# start interactive shell
CMD ["/bin/bash","-i"]

EOT

echo "--> Created Dockerfile. You can adapt it to your needs, or"

echo ""
read -p "   build the Docker image directly? [y/N]" BUILD
echo ""
if test x$BUILD = xy -o x$BUILD = xY; then
    echo "Building Docker image... this may take several minutes."
    cd ..
    docker build -f docker/Dockerfile -t $DOCKER_TAG .
    echo ""
    echo "Successfully built docker image: $DOCKER_TAG. Have a look at docker/README.md."
    echo "Check the container running docker run -it $DOCKER_TAG /bin/bash "
    echo "in the same directory as the Dockerfile."
    echo "And try using the convenience script pubtable_$DOCKER_TAG, see docker/README.md."
else
    cd ..
    echo "You can build your Docker image later by running docker build -f docker/Dockerfile -t $DOCKER_TAG ."
    echo "in your module directory, i.e. above the docker folder containing Dockerfile."
fi
