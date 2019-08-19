#!/usr/bin/env python

"""Simple example of Galaxy Tool
"""

import sys

def hello(names, out=sys.stdout):
    print >>out, "Hello %s!"%(",".join(names),)

if __name__=="__main__":
    outpath = sys.argv[1]
    with open(outpath, "w") as outfile:
        hello(sys.argv[2:], out=outfile)

