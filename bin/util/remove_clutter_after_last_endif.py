#!/usr/bin/env python3
import os

# replace everything after last #endif with new line
def clearAfterLastEndIf(filename):
    with open(filename, 'r') as header:
        split = header.read().split("#endif")
        split[-1] = '\n'
    with open(filename, 'w') as header:
        header.write("#endif".join(split))

for root, _, files in os.walk(os.getcwd()):
    for file in files:
        if file.endswith(".hh"):
            clearAfterLastEndIf(os.path.join(root, file))
