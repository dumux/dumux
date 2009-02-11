#!/usr/bin/python

# script to add an '#include "config.h"' as the first include all
# files which are given on as command line arguments, but only if the
# first is not config.h already.
#
# Typical usage: find -name "*.cc" | xargs ./addConfigH.py

import sys
import shutil

for filename in sys.argv[1:]:
    print "Processing file %s"%filename
    added = False
    tmpFile = open("/tmp/addConfigH.tmp", "w")
    for line in open(filename).xreadlines():
        if not added and line.find("#include") >= 0:
            if line.find("config.h") < 0:
                tmpFile.write("#include \"config.h\"\n")
            added = True
        tmpFile.write(line)
    tmpFile.close()
    shutil.move("/tmp/addConfigH.tmp", filename)
