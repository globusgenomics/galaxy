#!/usr/bin/env python
"""
Converts binary output file to a tabulated text file
"""

import sys
import os
import argparse
import string
import struct
import numpy
import re

#Parse arguments
parser = argparse.ArgumentParser(description='Convert CRU3D global states file to text.')
parser.add_argument('--output_file')
parser.add_argument('--simulation')
parser.add_argument('--ryrs_per_cell')
parser.add_argument('--cyto_volume')
parser.add_argument('--ryr_thresh')
args = parser.parse_args()

#Get output files
dirlist = os.listdir(args.simulation)
state_files = []
p = re.compile("states_")
for f in dirlist:
    if p.match(f):
        state_files.append(f)

#Compute fidelity and leak
N_Sparks = 0;
Leak_Spark = 0;
Leak_Invis = 0;
for f in state_files:
    data = numpy.fromfile(args.simulation+'/'+f,numpy.dtype('<f'))
    M = numpy.reshape(data,(-1,7))
    n_open_max = M[:,2].max()
    t = M[:,0]
    t_steps = t[1:-1]-t[0:-2]
    leak = numpy.dot(t_steps,M[0:-2,1])
    if int(n_open_max) >= int(args.ryr_thresh):
        N_Sparks += 1
        Leak_Spark += leak
    else:
        Leak_Invis += leak
P_Spark = float(N_Sparks)/float(len(state_files))
Leak_Spark /= N_Sparks
Leak_Invis /= len(state_files)-N_Sparks

R_Open = (0.1107e-3)*pow(0.1,2.1) * 1000 #1/s
R_Quark = R_Open * float(args.ryrs_per_cell)
R_NonSpark = R_Quark*(1-P_Spark)
R_Spark = R_Quark * P_Spark

Spark_Leak_pA_ms = R_Spark * Leak_Spark;
Invis_Leak_pA_ms = R_NonSpark * Leak_Invis;
Spark_Leak_uM_s = Spark_Leak_pA_ms * (1e-3) * (1e-6) / (2*96485.3365*float(args.cyto_volume))
Invis_Leak_uM_s = Invis_Leak_pA_ms * (1e-3) * (1e-6) / (2*96485.3365*float(args.cyto_volume))

f_out = open(args.output_file,'w')
f_out.write('Number of trials: %d<br>' % len(state_files))
f_out.write('Number of sparks: %d<br>' % N_Sparks)
f_out.write('<b>Spark fidelity: %0.4f</b><br>' % P_Spark)
f_out.write('Non-spark rate: %0.2f (1/cell-s)<br>' % R_NonSpark)
f_out.write('<b>Spark rate: %0.2f (1/cell-s)</b><br>' % R_Spark)
f_out.write('Mean leak per non-spark (invisible): %0.4f (pA-ms)<br>' % Leak_Invis)
f_out.write('Mean leak per spark: %0.4f (pA-ms)<br>' % Leak_Spark)
f_out.write('Non-spark leak rate: %0.4f (uM/s)<br>' % Invis_Leak_uM_s)
f_out.write('<b>Spark leak rate: %0.4f (uM/s)</b><br>' % Spark_Leak_uM_s)
f_out.close()
