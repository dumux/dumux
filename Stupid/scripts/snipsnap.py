#! /usr/bin/python

import sys

openFiles = {}

outFile = sys.stdout

n = 0
i = 0
for curLine in sys.stdin.xreadlines():
    if curLine.startswith("snip:"):
        i += 1
        key = curLine.split(":")[1].strip()
        if not openFiles.has_key(key):
            n+=1
            fileName = "%dstream_%s.csv"%(n, key)
            openFiles[key] = file(fileName, "w");
        outFile = openFiles[key]
        continue
    elif curLine.startswith("snap"):
        outFile = sys.stdout
        continue
    outFile.write(curLine)
