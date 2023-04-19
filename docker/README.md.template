# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# Readme for the dumux module ${modName}

You created the Docker image ${dockerTag}. Next steps:

* Try your container by running docker_${dockerTag}.sh open
  See below for instructions how to share files with the host system.
* Push the docker image to DockerHub or the GitLab Docker registry of your dumux module.
  Look at the Registry tab of your dumux module for help.

* Replace the image name in docker_${dockerTag}.sh with the actual image name
  e.g. git.iws.uni-stuttgart.de:4567/dumux-pub/koch2017a.

* [Optional] Add the Dockerfile to the git repository.

* Add the docker_${dockerTag}.sh script to the git repository (this is for the user)
  and add the following lines to your README.md:

Using dumux module ${modName} with docker
=============================================

In order to run simulations of this module look
at the convenience script docker_${dockerTag}.sh.
First download the script from the git repository.

The simplest way is to spin up a container
is creating a new folder "dumux"
$$ mkdir dumux
change to the new folder
$$ cd dumux
and open the module image by running
$$ docker_${dockerTag}.sh open

The container will spin up. It will mount the "dumux"
directory into the container at /dumux/shared. Put files
in this folder to share them with the host machine. This
could be e.g. VTK files produced by the simulation that
you want to visualize on the host machine.
