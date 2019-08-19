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
from scipy import interpolate

#Parse arguments
parser = argparse.ArgumentParser(description='Estimates ECC gain.')
parser.add_argument('--output_file')
parser.add_argument('--fluxfile')
parser.add_argument('--simulation')
parser.add_argument('--n_cru')
args = parser.parse_args()

#Get output files
dirlist = os.listdir(args.simulation)
state_files = []
p = re.compile("states_")
for f in dirlist:
    if p.match(f):
        state_files.append(f)

#Aggregate fluxes
Release_Total = []
Trigger_Total = []
Time = []
for f in state_files:
    data = numpy.fromfile(args.simulation+'/'+f,numpy.dtype('<f'))
    M = numpy.reshape(data,(-1,7))
    if len(Time)==0:
        Time = M[:,0]
        Release_Total = M[:,1]
        Trigger_Total = M[:,4]
    else:
        interp_release = interpolate.interp1d(M[:,0],M[:,1],bounds_error=False,fill_value=0)
        interp_trigger = interpolate.interp1d(M[:,0],M[:,4],bounds_error=False,fill_value=0)
        Release_Total += interp_release(Time)
        Trigger_Total += interp_trigger(Time)
Release_Total = Release_Total * (float(args.n_cru) / len(state_files))
Trigger_Total = Trigger_Total * (float(args.n_cru) / len(state_files))
Peak_Release = Release_Total.max()
Peak_Trigger = Trigger_Total.max()

if Peak_Trigger==0:
	Gain = 0
else:
	Gain = Peak_Release/Peak_Trigger


#Write computation results
f_out = open(args.output_file,'w')
f_out.write('Peak release : %0.2f (pA)<br>' % Peak_Release)
f_out.write('Peak trigger: %0.2f (pA)<br>' % Peak_Trigger)
f_out.write('<b>ECC Gain: %0.4f</b><br>' % Gain)
f_out.close()

#Write fluxes to tab-delimited file
data = numpy.column_stack((Time,Trigger_Total,Release_Total))
numpy.savetxt(args.fluxfile,data,delimiter="\t")
