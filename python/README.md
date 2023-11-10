# Python bindings

The DuMu<sup>x</sup> Python bindings have two objectives
* making it possible to use DuMu<sup>x</sup> from Python
* use Python code in DuMu<sup>x</sup>

The current Python bindings are far from complete and only cover a
a small subset of DuMu<sup>x</sup> functionality
(see [test/python](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/test/python)
for some example Python scripts).
The binding are experimental until further notice which means
the API might undergo unannounced changes breaking backwards
compatibility. Track changes by regularly checking for current merge requests
and newly merged Dumux commits on GitLab. Python bindings require Python 3.
Feedback over the mailing list or the issue tracker is highly welcome.

[TOC]

## Installation of Python bindings

Using Python bindings require at least Dune core modules version 2.9
and at least DuMu<sup>x</sup> version 3.7. Python bindings are
configured **automatically** with `dunecontrol` if Python is found on your system
and DuMux and upstream modules are built using the default options in `cmake.opts`.
Nothing special needs to be done if you follow the installation instructions.

### Running a test

After configuring and installing DuMu<sup>x</sup> with `dunecontrol`
(or the installscript) you are ready to run the Python tests.

From the top-level dumux folder,
run your first DuMu<sup>x</sup> Python test using the helper script

```
./build-cmake/run-in-dune-env dumux/test/python/test_gridgeometry.py
```

See below what this script does and how you can get better control over
what is exactly happening.

The Python bindings are based on just-in-time compilation of C++ code,
so the first execution might take considerably longer than subsequent executions.
In particular, compiling the Dune grid interface may take a few minutes.

### Running all tests

To check that everything works correctly, you can
also run all currently existing DuMu<sup>x</sup> Python tests with
```
cd build-cmake
ctest -L python
```

## Setup with a Python virtual environment

When configuring `dune-common`, by default, an internal Python virtual environment setup is created
in the build folder of `dune-common` for the Python bindings. Hence, in the description above,
you were already running the script in a virtual environment. In this virtual environment,
all Dune module Python bindings are automatically installed
in editable mode (symlinked) when running `dunecontrol`.
After running `dunecontrol` the internal virtual environment can be
activated with

```
source ./dune-common/build-cmake/dune-env/bin/activate
```

Then you can run Python tests *without* the helper script like this

```
python3 dumux/test/python/test_gridgeometry.py
```

To have better control of the virtual environment, you can create and
activate a new virtual environment yourself before running `dunecontrol`.
The CMake build system of `dune-common` detects the activated virtual environment
and installs Python bindings into that virtual environment (again in editable mode (symlinked)).
This looks like this:

```
python3 -m venv venv
source ./venv/bin/activate
./dune-common/bin/bexec rm -r CMakeFiles CMakeCache.txt
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts all
```

This installs everything necessary into your virtual environment `venv`.
(You may choose any other name of course.) With activated virtual
environment, you can run Python script using Python bindings like this

```
python3 dumux/test/python/test_gridgeometry.py
```

## DuMux Python example

As an example of how to use the Python bindings have a look at this
[Python script solving a problem with one-phase flow in porous media](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/test/python/test_1p.py).

<!--DOXYGEN_ONLY @include test_1p.py-->


##  Development of Python bindings

You can install all development requirements with

```
python3 -m pip install -r requirements.txt
```

The `requirements.txt` can be found in the DuMux repository in the top level.
We recommend to use a virtual environment for development as described above.

All Python files are linted by the tool [`black`](https://pypi.org/project/black/).
You can install `black` with `pip install black` and run it from the dumux top-directory

```
black ./python
```

You can also run it on a specific file (replace `./python` by file name)
This will automatically format the Python files. Run black before every commit changing Python files.
The CI pipeline prevents merging code that would reformat when running `black`.

The `dumux` Python module should be get a score of `10` from
the tool [`pylint`](https://pypi.org/project/pylint/).
You can install `pylint` with `pip install pylint` and run it from the dumux top-directory

```
pylint build-cmake/python/dumux
```

The `pylint` configuration file `dumux/.pylintrc` can
be used to configure `pylint`. Some exceptions or other parameters than the default
might be sensible in the future but generally advice given by `pylint` leads to better code.
Different from `black`, `pylint` does no itself fix the code, you need to do this yourself.
Always run `black` before checking `pylint`.
We also run `flake8` on Python code. The CI pipeline automatically checks
that `pylint` and `flake8` return without warning.
