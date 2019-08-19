"""
Given an output file path and a directory creates an output file that contains that directory.

"""

import sys

with open(sys.argv[1], "w") as outfile:
    outfile.write(sys.argv[2] + "\n")
