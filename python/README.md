# Python bindings for DuMu<sup>x</sup>

The Python bindings have two objectives
* making it possible to use DuMu<sup>x</sup> from Python
* use Python code in DuMu<sup>x</sup>

The binding are experimental until further notice which means
the API might undergo unannounced changes breaking backwards
compatibility. Track changes by regularly checking for current merge requests
and newly merged Dumux commits on GitLab. Python bindings require Python 3.

##  Installation

Python bindings are part of the dune core modules >= version 2.7.
We recommend to use the Dune modules on the master branch to try out
Python bindings.

### Example development setup

This section is rather new and experimental. Please help improving
this documentation in the future.

Checkout the `master` branch of the Dune core modules and DuMu<sup>x</sup>

```
git clone https://gitlab.dune-project.org/core/dune-common
git clone https://gitlab.dune-project.org/core/dune-geometry
git clone https://gitlab.dune-project.org/core/dune-grid
git clone https://gitlab.dune-project.org/core/dune-localfunctions
git clone https://gitlab.dune-project.org/core/dune-istl
git clone https://git.iws.uni-stuttgart.de/dumux-repositories/dumux

cp dumux/cmake.opts .
```

Enable the Python bindings by adding the flag in `cmake.opts` (see comments inside the `.opts` file).
Then run dunecontrol which also builds the Dune Python bindings.

```
./dune-common/bin/dunecontrol --opts=cmake.opts all
```

Add the Python binding modules to your Python path like this and install them:

```
export PYTHONPATH=$(pwd)/dune-common/build-cmake/python:$(pwd)/dune-grid/build-cmake/python:$(pwd)/dune-geometry/build-cmake/python:$(pwd)/dune-istl/build-cmake/python:$(pwd)/dune-localfunctions/build-cmake/python:$(pwd)/dumux/build-cmake/python
python3 dune-common/bin/setup-dunepy.py --opts=cmake.opts install
```

If you are getting error with loading MPI in Python you might need to preload the MPI library

```
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
```

Replace the path with the path to the MPI library on your system.
Run your first DuMu<sup>x</sup> Python test

```
python3 dumux/test/python/test_py_gridgeometry.py
```

The Python bindings are based on just-in-time compilation of C++ code,
so the first execution might take considerably longer than subsequent executions.
In particular compiling the Dune grid interface may a few minutes.

You can run all currently existing DuMu<sup>x</sup> Python tests with
```
cd dumux/build-cmake
ctest -L python
```
