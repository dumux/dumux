# Python bindings for DuMu<sup>x</sup>

The Python bindings have two objectives
* making it possible to use DuMu<sup>x</sup> from Python
* use Python code in DuMu<sup>x</sup>

The binding are experimental until further notice which means
the API might undergo unannounced changes breaking backwards
compatibility. Track changes by regularly checking for current merge requests
and newly merged Dumux commits on GitLab. Python bindings require Python 3.

##  Installation

Python bindings are part of the dune core modules >= version 2.8 and enabled by default for
versions >= 2.9 and DuMu<sup>x</sup> >= version 3.7.

### Running a test

After configuring and installing DuMu<sup>x</sup> with `dunecontrol` (or the installscript) you are
ready to run the python tests.

Run your first DuMu<sup>x</sup> Python test using the helper script

```
./build-cmake/run-in-dune-env dumux/test/python/test_gridgeometry.py
```

The Python bindings are based on just-in-time compilation of C++ code,
so the first execution might take considerably longer than subsequent executions.
In particular compiling the Dune grid interface may a few minutes.

You can run all currently existing DuMu<sup>x</sup> Python tests with
```
cd dumux/build-cmake
ctest -L python
```

### Example development setup (Dune 2.9)

This section is rather new and experimental. Please help improving
this documentation in the future.

When configuring dune-common, by default, an internal Python virtual environment setup is created
for the Python bindings in which all following modules' Python bindings are automatically installed
in editable mode (symlinked). After running `dunecontrol` this internal virtual environment is
activated with

```
source ./dune-common/build-cmake/dune-env/bin/activate
```

Alternatively, to have better control of the virtual environment, you can create and activate a new virtual environment, before running `dunecontrol`. The CMake build system of `dune-common` detects the activated virtual environment and installs Python bindings into that virtual environment (again in editable mode (symlinked)).

```
python3 -m venv venv
source ./venv/bin/activate
./dune-common/bin/bexec rm -r CMakeFiles CMakeCache.txt
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts all
```

The helper script `run-in-dune-env` activates the virtual environment and executes python
scripts passed as arguments. Instead you can manually activate your chosen virtual
environment (rerunning dunecontrol is not required) and run test scripts directly.

```
source ./path-to-env/bin/activate
python3 dumux/test/python/test_gridgeometry.py
```

##  Development

All Python files should be linted by the tool [`black`](https://pypi.org/project/black/).
You can install `black` with `pip install black` and run it from the dumux top-directory

```
black ./python
```

You can also run it on a specific file (replace `./python` by file name)
This will automatically format the Python files. Run black before every commit changing Python files.

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
