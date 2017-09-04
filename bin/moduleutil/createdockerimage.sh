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

echo "Creating a Dockerfile from module $MODULE_NAME..."

# if INSTALL_DEP_SCRIPT was not specified use the default
INSTALL_DEP_SCRIPT=${1:-install$MODULE_NAME.sh}
echo "Using install script: $INSTALL_DEP_SCRIPT to install dependencies for module $MODULE_NAME."

# make it executable
chmod +x $INSTALL_DEP_SCRIPT

# check if the Dockerfile already exists
if [ -e Dockerfile ]; then
    read -p "   A Dockerfile already exists. Continue anyway? [y/N] " DELETE
    if test x$DELETE = xy -o x$DELETE = xY; then
        rm Dockerfile
        if [ -e .dockerignore ]; then
            rm .dockerignore
        fi
        echo "Deleted old Dockerfile."
    else
        echo "Abort..."
        exit 1
    fi
fi

touch .dockerignore
echo "build*" >> .dockerignore
echo "CMakeFiles" >> .dockerignore

touch Dockerfile
cat <<EOT >> Dockerfile
# $MODULE_NAME docker container
FROM ubuntu:latest
MAINTAINER $MODULE_MAINTAINER

# get the package list
RUN apt-get update && apt-get dist-upgrade --no-install-recommends --yes && \\
    apt-get install --no-install-recommends --yes \\
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
    && rm -rf /var/lib/apt/lists/*

# create a dumux user
RUN useradd -m --home-dir /dumux dumux

# copy the module itself
ADD . /dumux/$MODULE_NAME
RUN chown -R dumux /dumux/$MODULE_NAME

# add user to video group for graphics
RUN usermod -a -G video dumux
USER dumux
WORKDIR /dumux

# install dependencies
COPY $INSTALL_DEP_SCRIPT /dumux/$INSTALL_DEP_SCRIPT
RUN ./$INSTALL_DEP_SCRIPT && rm -f /dumux/$INSTALL_DEP_SCRIPT

# configure module
RUN ./dune-common/bin/dunecontrol --opts=./dumux/optim.opts all

# build doxygen documentation
RUN cd $MODULE_NAME/build-cmake && make doc

# run bash shell
CMD ["/bin/bash"]

EOT

echo "Created Dockerfile. You can adapt it to your needs, or"

read -p "   Build Docker image directly? [y/N]" BUILD
if test x$BUILD = xy -o x$BUILD = xY; then
    echo "Building Docker image... this may take several minutes."
    DOCKER_TAG=$(echo $MODULE_NAME | tr '[:upper:]' '[:lower:]')
    docker build -t $DOCKER_TAG .
    echo "Successfully built docker image: $DOCKER_TAG."
    echo "Check the container running docker run -it $DOCKER_TAG /bin/bash"
else
    echo "Abort. You can build your Docker image by running docker build -t $DOCKER_TAG ."
    echo "in the same directory as the Dockerfile (your module directory)."
fi
