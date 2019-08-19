#!/usr/bin/python
import optparse
import sys
import os
from pprint import pprint
import shutil
import glob

parser = optparse.OptionParser()
parser.add_option( '-i', '--input-directory', dest='indir', help='Directory object path' )
parser.add_option( '-d', '--output-directory', dest='outdir', help='Output directory' )
parser.add_option( '-o', '--output', dest="output_file", help="output file" )
parser.add_option( '-l', '--list-file', dest="sample_list_file", help="list file of samples to keep" )
parser.add_option( '-t', '--list-string', dest="sample_list_string", help="string of samples to keep" )
(options, args) = parser.parse_args()

# get samples from list
samples_to_keep = []
if options.sample_list_file:
    fh = open(options.sample_list_file, "r")
    for line in fh:
        samples_to_keep.append(line.rstrip("\n"))
    fh.close()
elif options.sample_list_string:
    for sample in options.sample_list_string.split(","):
        samples_to_keep.append(sample.replace(" ", ""))


# get paths to copy to new directory
files = []
for sample in samples_to_keep:
    files.extend(glob.glob("%s/%s*" % (options.indir, sample)))

# copy files to new directory
if not os.path.exists(options.outdir):
    os.makedirs(options.outdir)

for file in files:
    shutil.copyfile(file, "%s/%s" % (options.outdir, os.path.basename(file)))

os.system("ls -l %s > %s" % (options.outdir, options.output_file))


