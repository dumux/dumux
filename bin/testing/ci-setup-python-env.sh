#!/bin/bash

if [ -d "/dune/modules/dune-common/build-cmake/dune-env" ]; then
    # Use internal venv of DUNE
    echo "Activating the Python virtual environment of dune-common"
    source /dune/modules/dune-common/build-cmake/dune-env/bin/activate
else
    if [ -L /dune/bin/setup-python ] && [ -e /dune/bin/setup-python ] ; then
        dunecontrol bexec "echo -n :\$(pwd)/python >> $(pwd)/pythonpath.txt"
        export PYTHONPATH=$PYTHONPATH$(cat pythonpath.txt)
        rm pythonpath.txt
        dunecontrol --current --opts=$DUNE_OPTS_FILE make install_python
    fi
fi
