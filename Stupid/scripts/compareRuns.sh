#!/bin/bash

FILTER=".*" ~/Stupid/scripts/runsim.sh lenhard 1000 10 | head -n 500 > singularius.txt
~/Stupid/scripts/runonteleloch.sh lenhard 1000 10 | head -n 500 > teleloch.txt
kdiff3 singularius.txt teleloch.txt
