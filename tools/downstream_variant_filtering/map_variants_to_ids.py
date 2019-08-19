#!/usr/bin/python
import os, sys, tempfile, optparse

parser = optparse.OptionParser()
parser.add_option( '-v', '--vcf', dest="vcf_input", help='vcf' )
parser.add_option( '-i', '--id', dest="id_input", help='id list' )
parser.add_option( '-o', '--output', dest="output_bed", help="bed_output" )
(options, args) = parser.parse_args()

id_dict = {}
fh = open(options.id_input, "r")
for line in fh:
    id_dict[line.rstrip("\n")] = 0
fh.close()

fh = open(options.vcf_input, "r")
for line in fh:
    values = line.split("\t")
    if values[4] in id_dict:
        print "IS_IN_FINAL\t%s" % line
        id_dict[values[4]] = 1
    else:
        print "NOT_IN_FINAL\t%s" % line
fh.close()

for i in id_dict:
    if id_dict[i] == 0:
        print "NOT_IN_FILTERED\t%s" % i

