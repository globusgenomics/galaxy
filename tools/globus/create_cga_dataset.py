#!/usr/bin/env python
"""
"""

import ConfigParser
import glob
from optparse import OptionParser
import os
import sys
import json

# begin: imported for side affects
import galaxy
from galaxy import eggs
import galaxy.model
from galaxy.datatypes import sniff
# end: imported for side affects

from galaxy.datatypes.genetics import CgaData

# TODO: more reliable way of finding ini file?
CONFIG_FILE=os.path.join(
    os.path.dirname(
    os.path.dirname(
    os.path.dirname(
    os.path.abspath(galaxy.__file__)))),
    "universe_wsgi.ini")

def get_genomes():
    """Return option list of genome directories.
    """
    options = []

    #dir = "/Users/steder/T2DTest" # TODO pull this from galaxy.app.config
    genome_dir = None
    cp = ConfigParser.SafeConfigParser()
    with open(CONFIG_FILE, "r") as config_file:
        cp.readfp(config_file)
        genome_dir = cp.get("galaxy:tools", "complete_genomics_root")

    listing = os.listdir(genome_dir) 
    for path in listing:
        options.append((path, os.path.join(
            genome_dir, path), True))
        
    if len(options) >= 1:
        options.insert(0,('None','None',True))
    else:
        options = [('None','no genomes found',False),]
    return options


def write_dataset(outfile=None, input_directory=None,
                  output_directory=None):
    """Create a primary dataset file by reading input_directory.

    Also attempts to create the extra files directory expected
    by composite datasets by symbolicly linking input_directory
    to output_directory.
    """
    # instead of creating the output directory we'll
    # create a symbolic link to the input directory.
    # TODO: uncomment this when we know that galaxy won't follow the link
    # and accidentally delete our real dataset
    #os.symlink(input_directory, output_directory)

    # now create the primary dataset file:
    d = CgaData.dataset_dict_from_directory(input_directory)
    with open(outfile, "w") as out:
        json.dump(d, out)


if __name__ == "__main__":
    """
    """
    parser = OptionParser(usage=__doc__, version="%prog 0.01")
    parser.add_option("-o","--outf",dest="outf",
      help="Output file", default = 'fakeped')
    parser.add_option("-p","--outpath",dest="outpath",
      help="Path for output files", default = './')
    parser.add_option("-i", "--inpath", dest="inpath",
      help="Input data directory", default="./")
    (options,args) = parser.parse_args()
    write_dataset(outfile=options.outf,
                  input_directory=options.inpath,
                  output_directory=options.outpath)


        
