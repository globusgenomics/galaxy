#!/usr/bin/env python
"""
Converts binary output file to a tabulated text file
"""

import sys
import os
import subprocess
import argparse
import string
import struct

#Parse arguments
parser = argparse.ArgumentParser(description='Convert CRU3D global states file to text.')
parser.add_argument('--output_file')
parser.add_argument('--simulation')
args = parser.parse_args()

#Open states file
f_bin = open(args.simulation+'/states_0_0.out','rb')

#Send output to file
f_txt = open(args.output_file,'w')
#f_txt.write('time_ms\ti_ryr_pA\tryrs_open\ti_lcc_pA\tlccs_open\tca_jsr_uM\tca_ss_uM\n')
while 1:
    x = []
    for i in range(7):
        x.append(f_bin.read(4))
    if not x[0]:
        break
    for i in range(7):
        data = struct.unpack('f',x[i])
        f_txt.write('%f\t' % data[0])
    f_txt.write("\n")

f_txt.close()

f_bin.close();
