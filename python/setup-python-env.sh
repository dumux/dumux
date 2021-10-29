#!/bin/bash

{ # try
    # Use internal venv of DUNE
    source ./dune-common/build-cmake/dune-env/bin/activate &&
} || { # catch
    # Install the python packages into custom venv (only needed for dune 2.8)
    python3 -m venv dune-venv
    source dune-venv/bin/activate
    ./dune-common/bin/dunecontrol --opts=cmake.opts make install_python
}
