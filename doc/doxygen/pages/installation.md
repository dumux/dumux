# Installation

[TOC]

## Prerequisites

DuMu<sup>x</sup> builds and runs on **Linux** and **Mac** operating systems.
Installation on Windows is possible but it is definitely not something to try if you are a beginner.
If you use Windows, we recommend the [Ubuntu bash on Windows (WSL)](https://msdn.microsoft.com/commandline/wsl/)
and a [guide to using Dumux with the WSL can be found in the wiki](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/wikis/Installation-of-DuMux-inside-the-Windows-Subsystem-for-Linux-(WSL)).
Alternatively, you can try to employ MinGW, Cygwin or a Linux Virtual Machine.

In order to build DuMu<sup>x</sup> you need at least the following software:

* Standard-compliant C++17 compiler supporting the common feature set supported by GCC 9.3 and Clang 10
* CMake 3.16 or newer
* pkg-config
* The DUNE core modules (>= 2.9), see installation instructions below

The following software is recommended but optional:

* MPI (either OpenMPI, lam, or mpich; only OpenMPI support is automatically tested)
* SuiteSparse (for the direct solver UMFPack)
* ParaView (to visualize the results)
* A web browser (to access the GitLab instance and README files)
* Python >= 3.7 with `numpy` (to execute a number of different scripts used for installation, testing, post-processing, etc.)


## 1. Installing DuMux and DUNE via script

This installation method requires

* `wget` (for automatically downloading dependencies)
* `git` (for cloning the code repositories)
* `python3` (>= 3.7)

We provide you with a Python script that facilitates setting up a Dune/DuMux directory
tree and configures all modules using CMake. Download [installdumux.py](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/bin/installdumux.py) and place it in a directory where you want to install Dune/DuMux, for example by running

    wget https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/bin/installdumux.py

Executing the script creates a folder "dumux" in which all dependencies will be downloaded and built.
Run the script by typing:

    python3 installdumux.py

Follow the instructions printed after the script has completed successfully to verify everything works as expected.

In brief, the script performs the following steps:

* create a folder "dumux"
* clone the DUNE core modules in the current release version
* clone the DuMux module
* configure DUNE and DuMux with CMake by using the script `dune-common/bin/dunecontrol` and the options in `dumux/cmake.opts`
* build the DUNE and DuMux libraries using CMake and generated GNU Makefiles

Note that this process can take several minutes. The next section will guide
you through the same process providing the necessary commands for executing each step.

## 2. Manually installing DUNE and Dumux from source

### 2.1 Obtaining the DUNE core modules

DuMux depends on [Dune](https://dune-project.org/).
Required are the libraries/modules `dune-common`, `dune-geometry`,
`dune-grid`, `dune-localfunctions`, and `dune-istl`.
You can obtain the required Dune modules in form of binary packages
for Debian, Ubuntu and openSUSE, see the [Dune binary packages](http://www.dune-project.org/binary/).
The Dune releases can also be obtained as [tarballs](https://www.dune-project.org/releases/).
Alternatively and recommended and described below, you can use [git](https://www.dune-project.org/dev/downloadgit/)
to download the source code and configure and build with the `dunecontrol` script.
To clone the Dune core modules, run:

    for module in common geometry grid localfunctions istl; do
      git clone -b releases/2.9 https://gitlab.dune-project.org/core/dune-$module.git
    done

### 2.2 Obtaining the Dumux source code

To clone the latest 3.8 release version, run

    git clone -b releases/3.8 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git

The master branch (developer version) can be cloned with

    git clone https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git

Alternatively, it is also possible to download release tarballs
from [DaRUS](https://doi.org/10.18419/darus-3788) or
from [GitLab](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/releases).

### 2.3 Configure and build

After obtaining all modules (at least `dune-common`,`dune-geometry`,`dune-grid`,`dune-istl`,`dune-localfunctions` and `dumux`),
DuMux is built together with other Dune modules. Assuming that the Dune core modules and DuMux are folders in the same directory

    DUMUX (installation folder)
    |- dune-common
    |- dune-geometry
    |- dune-grid
    |- dune-istl
    |- dune-localfunctions
    |- dumux

you can configure and build the module stack with the `dunecontrol` helper script:

    ./dune-common/bin/dunecontrol --opts=dumux/cmake.opts all

In case you have obtained the DUNE modules via a package manager,
`dunecontrol` should be an available program in your path environment
(replace `./dune-common/bin/dunecontrol` by `dunecontrol`).

Running `dunecontrol` will create a build folder called
`build-cmake` in each of the module folders.
In case you want to build the module stack in a separate build folder use

    ./dune-common/bin/dunecontrol --opts=dumux/cmake.opts --builddir=$(pwd)/build all

More details on the DUNE build system can be found in the [Dune installation notes](http://www.dune-project.org/doc/installation/).

### Compiler options

Dune and DuMux rely heavily on compiler optimization.
The speed difference between running a compiler-optimized versus a non-optimized DuMux executable can easily exceed a factor of $10$.
The default `cmake.opts` in `dumux/cmake.opts` already enable compiler optimisations.
To use debug options use `-DCMAKE_BUILD_TYPE=Debug` instead of `-DCMAKE_BUILD_TYPE=Release` in the `cmake.opts`, or
add `set(CMAKE_BUILD_TYPE Debug)` to any `CMakeLists.txt` containing a test that you want
to compile non-optimised and with debug symbols enabled.

In case you are running macOS on a recent arm chip (M1, M2) your compiler might not support the flag `-march=native` yet.
You can replace this flag by a `mcpu=<cpu>` flag, where `<cpu>` is the best match of available option shown by
`clang --print-supported-cpus`.

Often it makes sense to create and keep around a custom option file for `dunecontrol` tailored
to your local setup. You can use `dumux/cmake.opts` as a starting point.

## 3. Install external dependencies via script

There are various external libraries and modules that provide additional functionality but are
not generally required to run DuMux.
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

| DuMux version     | Dune version              |
|-------------------|---------------------------|
| master            | 2.9, master               |
| 3.7, 3.8          | 2.9                       |
| 3.5, 3.6          | 2.8, 2.9                  |
| 3.3, 3.4          | 2.7, 2.8                  |
| 3.1, 3.2          | 2.6<sup>*</sup>, 2.7      |
| 3.0               | 2.6<sup>*</sup>, 2.7      |
| 2.10, 2.11, 2.12  | 2.4, 2.5, 2.6<sup>*</sup> |

<sup>*</sup> Use the most recent version of the git branch `releases/2.6` instead of the `2.6.0` tarballs.
