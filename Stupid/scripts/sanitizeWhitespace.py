#!/usr/bin/python
import sys

outbuf = ""
for curline in sys.stdin.xreadlines():
    curline = curline.rstrip().replace('\t', ' '*4)
    outbuf += curline + '\n'
    if curline:
        sys.stdout.write(outbuf)
        outbuf = ""

