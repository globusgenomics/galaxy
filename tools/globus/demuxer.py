#!/usr/bin/env python
"""Creates a simple single file dataset from a composite CGA dataset

"""

import json
from optparse import OptionParser
import os
import re
import sys


def parse_lane_and_part_numbers(path):
    # parse lane number and part number out of:
    # the reads file name reads_GS21910-FS3-L04_004.tsv.bz2
    basename = os.path.basename(path)
    sans_extension = basename.split(".")[0]
    lane_part_segment = sans_extension.split("-")[-1]
    lane_string, part_string = lane_part_segment.split("_")
    lane_number = int(lane_string.replace("L", ""))
    part_number = int(part_string)
    return lane_number, part_number


def link(datapath, outpath):
    os.unlink(outpath)
    os.symlink(datapath, outpath)


def parse_chromosome_number(path):
    basename = os.path.basename(path)
    pattern = ".*chrm([0-9]+)_.*"
    m = re.search(pattern, basename)
    number = None
    if m:
        try:
            number = int(m.groups()[0])
        except ValueError:
            print "%s does not appear to have a chromosome number!"%(path,)
    return number


def extract_file(data, extension, lane_number, part_number, outpath, chromosome_number=1):
    datasets = data["files"]
    print  "lane_number:", lane_number, "part_number:", part_number
    for dataset in datasets:
        path = dataset["path"]
        print  "path:", path
        try:
            lane_no, part_no = parse_lane_and_part_numbers(path)
        except:
            lane_no, part_no = None, None
        print  "lane_no, part_no:", lane_no, part_no
        chrom_no = parse_chromosome_number(path)
        if ((lane_no == lane_number and part_no == part_number and path.endswith(extension)) or
            chromosome_number == chrom_no):
            datapath = path
            print  "Creating a symbolic link to the data %s at outpath %s"%(
                datapath, outpath)
            if not os.path.exists(datapath):
                print  "Some of these tools add another extension (e.g.: samtools sam2bam) .bam extensions..."
                datapath += ".%s"%(extension,)
            print  "os.path.exists(outpath):", os.path.exists(outpath)
            link(datapath, outpath)
            return
    raise Exception("File not found!")


def convert_to_dataset(inpath, indir, outpath,
                       output_type, lane_number, part_number):
    try:
        lane_number = int(lane_number)
    except ValueError:
        print  "Invalid lane number!"
    try:
        part_number = int(lane_number)
    except ValueError:
        print  "Invalid part number!"

    with open(inpath, "r") as json_file:
        data = json.load(json_file)
        files = data["files"]
        #for f in files:
        #    print  f
        if output_type == "dat":
            parts = data["parts"]
            found_file = False
            for name in parts:
                part = parts[name]
                print "part[%s] = %s"%(name, part)
                lane_no, part_no = parse_lane_and_part_numbers(part["reads"])
                if lane_no == lane_number and part_no == part_number:
                    datapath = part["reads"]
                    print  "Creating a symbolic link to the data %s at outpath %s"%(
                        datapath, outpath)
                    link(datapath, outpath)
                    found_file = True
            if not found_file:
                raise Exception("Unable to find read or mapping file for this lane or part!")
        elif output_type == "sam":
            extract_file(data, "sam", lane_number, part_number, outpath)
        elif output_type == "bam":
            extract_file(data, "bam", lane_number, part_number, outpath)
        elif output_type == "boxplot.png":
            extract_file(data, "boxplot.png", lane_number, part_number, outpath)
        else:
            raise Exception("Unknown output_type '%s' specified!"%(output_type,))

                
if __name__ == "__main__":
    parser = OptionParser(usage=__doc__, version="%prog 0.01")
    parser.add_option("-o","--outpath",dest="outpath",
      help="output file path", default = 'output.dat')
    parser.add_option("-t", "--output-type", dest="output_type",
                      help="output dataset type", default="dat")
    parser.add_option("-p","--indir",dest="indir",
      help="path for input files", default = './')
    parser.add_option("-i", "--inpath", dest="inpath",
      help="path for input primary file", default="./")
    parser.add_option("-l", "--lane", dest="lane_number",
                      help="lane number to extract", default="01")
    parser.add_option("-r", "--part", dest="part_number",
                      help="part number to extract", default="01")
    
    (options,args) = parser.parse_args()
    convert_to_dataset(options.inpath, options.indir, options.outpath,
                   options.output_type, options.lane_number,
                   options.part_number)


