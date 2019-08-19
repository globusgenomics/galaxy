#!/usr/bin/python
import os, sys
from optparse import OptionParser

args=None
parser = OptionParser()
parser.add_option("--input-minid", dest="input_minid", help="Input minid")
parser.add_option("--ultra_mem", dest="output_ultra", help="Output high")
parser.add_option("--high_mem", dest="output_high", help="Output high")
parser.add_option("--med_mem", dest="output_med", help="Output med")
parser.add_option("--low_mem", dest="output_low", help="Output low")
options, args = parser.parse_args(args)

minid_template = "ark:/57799/b989Nsg8j0Vv9I"
path = "/opt/galaxy/tools/wgsa/"
low_template = "low.minid.template.batch.submit.txt"
med_template = "medium.minid.template.batch.submit.txt"
high_template = "high.minid.template.batch.submit.txt"
ultra_template = "ultra.high.minid.template.batch.submit.txt"

batch_files = [options.output_low, options.output_med, options.output_high, options.output_ultra]
index = 0
for template in [low_template, med_template, high_template, ultra_template]:
    filer = "%s/%s" % (path, template)
    filew = batch_files[index]
    fhr = open(filer, "r")
    fhw = open(filew, "w")

    for line in fhr:
        fhw.write(line.replace(minid_template, options.input_minid))
    index += 1

