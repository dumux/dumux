#!/bin/bash
# Adds all Python modules found in other Dune modules to the PYTHONPATH
#./dune-common/bin/dunecontrol bexec "echo -n :\$(pwd)/python >> $(pwd)/pythonpath.txt"
#export PYTHONPATH=$PYTHONPATH$(cat pythonpath.txt)
#rm pythonpath.txt
# Sets up the dune-py module for JIT compilation of Python binding codei (only needed for dune 2.8)
./dune-common/bin/dunecontrol --opts=cmake.opts make install_python
