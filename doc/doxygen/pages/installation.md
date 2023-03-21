# Installation

[TOC]

## Prerequisites

DuMu<sup>x</sup> builds and runs on **Linux** and **Mac** operating systems.
Installation on Windows is possible but it is definitely not something to try if you are a beginner.
If you use Windows, we recommend the [Ubuntu bash on Windows](https://msdn.microsoft.com/commandline/wsl/).
Alternatively, you can try to employ MinGW, Cygwin or a Linux Virtual Machine.

In order to build DuMu<sup>x</sup> you need at least the following software:

* C++17 compiler (GCC 9.3 or newer, Clang 10 or newer)
* CMake (CMake 3.13 or newer)
* pkg-config

Detailed information on supported compiler and CMake versions can be found in the [DuMux handbook](/docs/#handbook).
The following software is recommended but optional:

* MPI (either OpenMPI, lam, or mpich)
* ParaView (to visualize the results)
* a browser (to access the GitLab instance and README files)
* python3 with numpy (to execute a number of different scripts used for installation, testing, post-processing, etc.)


## 1. Installation via script

We provide you with a Python script that facilitates setting up a Dune/DuMux directory
tree and configures all modules using CMake. Download [installdumux.py](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/bin/installdumux.py) and place it in a directory where you want to install Dune/DuMux. The script will create a folder "DUMUX" in which all dependencies will be downloaded and built. Note that this installation method requires Python 3.X, `wget` (for automatically downloading dependencies) and `git` (for cloning the code repositories) to be installed on your system.

Run the script by typing:

    python3 installdumux.py

Follow the instructions printed after the script has completed successfully to verify everything works as expected.

For more detailed explanation on installation and building of DuMux, please have a look below or at the [handbook](/docs/#handbook) of the latest release.


## 2. Manual installation from source

DuMux depends on [Dune](https://dune-project.org/).
You can obtain the required Dune modules in form of binary packages for Debian, Ubuntu and openSUSE, see the [Dune binary packages](http://www.dune-project.org/binary/). The Dune releases can also be obtained as [tarballs](https://www.dune-project.org/releases/). Alternatively, you can use [Git](https://www.dune-project.org/dev/downloadgit/) and download the modules (recommended and described here). To clone the Dune core modules (minimum requirement), you can run:

    for module in common geometry grid localfunctions istl; do
      git clone -b releases/2.8 https://gitlab.dune-project.org/core/dune-$module.git
    done

Our Git repository (just like Dune's) offers anonymous read access. To clone the 3.6 release version, you can type:

    git clone -b releases/3.6 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git

You can clone the master branch (developer version) by typing:

    git clone https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git

Alternatively, it is also possible to download release tarballs from [Zenodo](https://doi.org/10.5281/zenodo.2479594) or from [GitLab](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/releases).

After obtaining all modules (at least the Dune core modules `dune-common`,`dune-geometry`,`dune-grid`,`dune-istl`,`dune-localfunctions` and `dumux`), DuMux is built together with other Dune modules. Assuming that the Dune core modules and DuMux are folders in the same directory

    installation folder
    |- dune-common
    |- dune-geometry
    |- dune-grid
    |- dune-istl
    |- dune-localfunctions
    |- dumux

you can configure and build the module stack with the `dunecontrol` helper script:

    ./dune-common/bin/dunecontrol --opts=dumux/cmake.opts all

This will create a separate build folder called `build-cmake` individually in each module. In case you want to build the module stack in a separate build folder use

    ./dune-common/bin/dunecontrol --opts=dumux/cmake.opts --builddir=$(pwd)/build all

More details on the Dune build system can be found in the [Dune installation notes](http://www.dune-project.org/doc/installation/). Dune and DuMux rely heavily on compiler optimization. The speed difference between running a compiler-optimized versus a non-optimized DuMux executable can easily exceed a factor of 10.
The default `cmake.opts` in `dumux/cmake.opts` already enable compiler optimisations.
To use debug options use `-DCMAKE_BUILD_TYPE=Debug` instead of `-DCMAKE_BUILD_TYPE=Release` in the `cmake.opts`, or
add `set(CMAKE_BUILD_TYPE Debug)` to any `CMakeLists.txt` containing a test that you want to compile non-optimised and with debug symbols enabled.


## 3. Install external dependencies via script

There are various external libraries and modules that provide additional functionality but are
not generally required to run DuMux. A list of external libraries and modules can be found in the [DuMux handbook](/docs/#handbook).
DuMux contains the script `installexternal.py` which allows you to install extension from your DuMux installation directory.

If you run the script with the option \-\-help

    python3 dumux/bin/installexternal.py --help

it will show a list of installable packages. Dune modules can be installed by adding the name of the package.
For example, running the script by typing:

    python3 dumux/bin/installexternal.py alugrid

downloads the ALUGrid module `dune-alugrid` in a separate folder in your installation directory. After the download
you can run following command to clean the CMake cache:

    ./dune-common/bin/dunecontrol bexec rm -r CMakeFiles CMakeCache.txt

Afterwards you can reconfigure and build DuMux with the `dunecontrol` script:

    ./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts all

If you install an external library which is not a Dune module, the script will
install the library in the directory `external`. If you install external libraries in an non-standard location you
must usually either set the root path of the package (e.g. `CGAL_ROOT=<path/to/folder>`) or add the path to the `CMAKE_PREFIX_PATH`
variable. These variables can be added in the option file `cmake.opts` before running `dunecontrol`. Compiler definitions are
passed to CMake with the `-D` prefix, e.g. `-DCGAL_ROOT=<path/to/folder>`, so variables passed via `cmake.opts` need `-D` at the beginning.
After you have installed an external library you should also clean the CMake cache and
reconfigure and build DuMux as mentioned above.


## Compatible versions of Dune and DuMux

Only the following Dune and DuMux versions are compatible:

<table class="table table-hover">
<tr>
    <th>DuMux version</th>
    <th>Dune version</th>
</tr>
<tr align="center">
    <td align="left">3.5, 3.6, 3.7-git</td>
    <td align="left">2.8, 2.9-git<sup>**</sup></td>
</tr>
<tr align="center">
    <td align="left">3.3, 3.4</td>
    <td align="left">2.7, 2.8</td>
</tr>
<tr align="center">
    <td align="left">3.1, 3.2</td>
    <td align="left">2.6<sup>*</sup>, 2.7</td>
</tr>
<tr align="center">
    <td align="left">3.0</td>
    <td align="left">2.6<sup>*</sup>, 2.7</td>
</tr>
<tr align="center">
    <td align="left">2.9, 2.10, 2.11, 2.12</td>
    <td align="left">2.4, 2.5, 2.6<sup>*</sup></td>
</tr>
<tr align="center">
    <td align="left">2.6, 2.7, 2.8</td>
    <td align="left">2.3, 2.4</td>
</tr>
<tr align="center">
    <td align="left">2.5</td>
    <td align="left">2.2, 2.3</td>
</tr>
</table>

<sup>*</sup> Please use the most recent corresponding Git branches `releases/2.6` instead of the `2.6.0` tarballs.

<sup>**</sup> Compatibility with Dune `2.9-git` is tested automatically only for DuMux `3.7-git`.
