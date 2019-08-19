#!/usr/bin/env python
"""
"""
from optparse import OptionParser
import json


def convert_to_txt(inpath, indir, outpath):
    with open(outpath, "w") as of:
        of.write("Reading inpath...\n")
        of.write("inpath: %s\n"%(inpath,))
        of.write("indir: %s\n"%(indir,))
        with open(inpath, "r") as json_file:
            data = json.load(json_file)
            files = data["files"]
            of.write("files:\n")
            for f in files:
                of.write("%s\n"%(f["basename"]))
            parts = data["parts"]
            of.write("\n")
            of.write("parts:\n")
            for name in parts:
                part = parts[name]
                of.write("reads %s mapping %s\n"%(part["reads"], part["mapping"]))

                
if __name__ == "__main__":
    """
    """
    parser = OptionParser(usage=__doc__, version="%prog 0.01")
    parser.add_option("-o","--outpath",dest="outpath",
      help="output file path", default = 'fakeped')
    parser.add_option("-p","--indir",dest="indir",
      help="path for input files", default = './')
    parser.add_option("-i", "--inpath", dest="inpath",
      help="path for input primary file", default="./")
    (options,args) = parser.parse_args()
    convert_to_txt(options.inpath, options.indir, options.outpath)


        
