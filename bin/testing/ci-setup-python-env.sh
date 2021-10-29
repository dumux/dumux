#!/bin/bash

if [ -d "./dune-common/build-cmake/dune-env" ]; then
    # Use internal venv of DUNE
    source ./dune-common/build-cmake/dune-env/bin/activate
else
    if [ -L /dune/bin/setup-python ] && [ -e /dune/bin/setup-python ] ; then
        dunecontrol bexec "echo -n :\$(pwd)/python >> $(pwd)/pythonpath.txt"
        export PYTHONPATH=$PYTHONPATH$(cat pythonpath.txt)
        rm pythonpath.txt
        setup-python --opts=$DUNE_OPTS_FILE install
    fi
fi
