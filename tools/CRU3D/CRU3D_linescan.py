#!/usr/bin/env python
"""
Generates CRU3D linescan.
"""

import sys
import os
import subprocess
import argparse
import string
import numpy
import scipy.misc
import numpy.random

#Parse arguments
parser = argparse.ArgumentParser(description='Generate CRU3D linescan.')
parser.add_argument('--root_dir')
parser.add_argument('--html_file')
parser.add_argument('--files_path')
parser.add_argument('--simulation')
parser.add_argument('--domain_x')
parser.add_argument('--domain_y')
parser.add_argument('--domain_z')
parser.add_argument('--c0')
parser.add_argument('--B_TOT_DYE')
parser.add_argument('--K_ON_DYE')
parser.add_argument('--K_OFF_DYE')
parser.add_argument('--num_scans')
parser.add_argument('--scan_interval')
parser.add_argument('--pixel_size')
parser.add_argument('--x_start')
parser.add_argument('--x_end')
parser.add_argument('--axis')
parser.add_argument('--offset_1')
parser.add_argument('--offset_2')
parser.add_argument('--var_x')
parser.add_argument('--var_y')
parser.add_argument('--var_z')
parser.add_argument('--noise_sigma')
args = parser.parse_args()

#Create output directory
#print "Creating directory: " + args.files_path
if not os.path.exists(args.files_path): 
    os.makedirs(args.files_path)

#Compute basal fluorescence
Dye_Kd = float(args.K_OFF_DYE)/float(args.K_ON_DYE);
F0 = float(args.B_TOT_DYE)*float(args.c0) / (float(args.c0) + Dye_Kd);

#Create XML
xml_content = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<analysis>
  <parameter>
    <symbol>scan_padding_res</symbol>
    <definition>Spacing between padding points [um]</definition>
    <label>Padding Spacing [um]</label>
    <value>0.30</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_padding_max</symbol>
    <definition>Upper-bound coordinate for domain padding</definition>
    <label>Upper-bound Domain Padding</label>
    <value>4.0</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_var_x</symbol>
    <definition>Variance of point spread function in X direction</definition>
    <label>PSF X Variance</label>
    <value>%(var_x)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_var_y</symbol>
    <definition>Variance of point spread function inY direction</definition>
    <label>PSF Y Variance</label>
    <value>%(var_y)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>fid_cyto_vol</symbol>
    <definition>Cytoplasmic volume in liters.</definition>
    <label>Cytosplasm volume [L]</label>
    <value>18e-12</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>outdir</symbol>
    <definition>Path to output files if needed</definition>
    <label>Output directory</label>
    <value>%(out_dir)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_offset_2</symbol>
    <definition>Second scan line offset parameter (Y if scan axis is Z, Z if scan axis is X or Y) [um]</definition>
    <label>Offset #2 [um]</label>
    <value>%(offset_2)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_domain_x</symbol>
    <definition>Size of mesh domain in X direction [um]</definition>
    <label>Mesh domain size (X) [um]</label>
    <value>%(domain_x)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>description</symbol>
    <definition>Brief description of the analysis</definition>
    <label>Description</label>
    <value></value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_c0</symbol>
    <definition>Resting Ca2+ concentration [uM]</definition>
    <label>Resting Ca2+ [uM]</label>
    <value>%(c0)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_n</symbol>
    <definition>Total number of scans to simulate</definition>
    <label>Number of scans</label>
    <value>%(num_scans)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>states_file</symbol>
    <definition>Path to simulation states file</definition>
    <label>States file</label>
    <value>states_0_0.out</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_offset_1</symbol>
    <definition>First scan line offset parameter (Y if scan axis is X, X if scan axis is Y or Z) [um]</definition>
    <label>Offset #1 [um]</label>
    <value>%(offset_1)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_axis</symbol>
    <definition>Axis along which to perform line scans</definition>
    <label>Scan Line Axis</label>
    <value>%(axis)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_gridfile_base</symbol>
    <definition>Simulation's grid output file base path/name</definition>
    <label>Grid file base</label>
    <value>%(grid_file)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>analysis_title</symbol>
    <definition>Title of analysis</definition>
    <label>Title</label>
    <value>DefaultLinescan</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_x_res_f</symbol>
    <definition>Length per pixel for fluorescence scans [um]</definition>
    <label>Spatial Fluorescence Resolution [um]</label>
    <value>%(pixel_size)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_padding_min</symbol>
    <definition>Lower-bound coordinate for domain padding</definition>
    <label>Lower-bound Domain Padding</label>
    <value>-4.0</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_x_min</symbol>
    <definition>Starting coordinate for scan</definition>
    <label>Scan Start Position</label>
    <value>%(x_start)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_x_max</symbol>
    <definition>Ending coordinate for scan</definition>
    <label>Scan End Position</label>
    <value>%(x_end)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_x_res_ca</symbol>
    <definition>Length per pixel for Ca2+ concentration profile [um]</definition>
    <label>Spatial Ca2+ Resolution [um]</label>
    <value>0.01</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_domain_y</symbol>
    <definition>Size of mesh domain in Y direction [um]</definition>
    <label>Mesh domain size (Y) [um]</label>
    <value>%(domain_y)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_domain_z</symbol>
    <definition>Size of mesh domain in Z direction [um]</definition>
    <label>Mesh domain size (Z) [um]</label>
    <value>%(domain_z)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_var_z</symbol>
    <definition>Variance of point spread function in Z direction</definition>
    <label>PSF Z Variance</label>
    <value>%(var_z)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>analysis_type</symbol>
    <definition>Type of analysis to perform</definition>
    <label>Analysis Type</label>
    <value>linescan</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_n_per_line</symbol>
    <definition>Number of grid files that correspond to a single scan (determines temporal resolution, e.g. if grid output interval is 1.0ms and files per scan is 2, then each scan will correspond to 2.0ms)</definition>
    <label>Grid files per scan</label>
    <value>%(scan_interval)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>scan_f0</symbol>
    <definition>Baseline fluorescence signal (resting Ca-Fluo4n concentration)</definition>
    <label>F0</label>
    <value>%(F0)s</value>
    <listable>false</listable>
  </parameter>
</analysis>

""" % {'domain_x':args.domain_x,
       'domain_y':args.domain_y,
       'domain_z':args.domain_z,
       'c0':args.c0,
       'num_scans':args.num_scans,
       'scan_interval':args.scan_interval,
       'pixel_size':args.pixel_size,
       'x_start':args.x_start,
       'x_end':args.x_end,
       'axis':args.axis,
       'offset_1':args.offset_1,
       'offset_2':args.offset_2,
       'var_x':args.var_x,
       'var_y':args.var_y,
       'var_z':args.var_z,
       'F0':F0,
       'out_dir':args.files_path,
       'grid_file':args.simulation+'/grid_0_0'}

#Write xml file
xml_path = args.files_path+"/Linescan_Parameters.xml"
f = open(xml_path,'w')
f.write(xml_content)
f.close()

#Run linescan program
base_dir = args.root_dir + '/tools/CRU3D/model/linescan';
proc = subprocess.Popen([base_dir+'/linescan','-param',xml_path],cwd=base_dir,stdout=subprocess.PIPE)

#Send output to file
output = proc.stdout.read()
f = open(args.files_path+"/log.txt",'w')
f.write(output)
f.close()
#outfile = open(args.stdout,'w')
#outfile.write(output);
#outfile.close();

#Create HTML file
f = open(args.html_file,'w')
#f.write(xml_content.replace("\n","<br>"))
#f.write(output.replace("\n","<br>"))
f.write('<a href="log.txt">Linescan Log</a><br>')
f.write('<a href="Linescan_Parameters.xml">Linescan Parameters</a><br>')
f.write('Fluorescence with optical blurring (<a href="scan_output.txt">txt</a>, no noise) (<a href="linescan.tif">tif</a>, with noise)<br>')
f.write('<img src="linescan.png" width="300px">')
#f.write('<a href="scan_output_CaF.txt">Output - Fluorescence without blurring</a><br>')
#f.write('<a href="scan_output_Ca.txt">Output - Ca2+ concentration</a><br>')
f.close()

#Read in output data
#I = numpy.loadtxt(args.files_path+'/scan_output.txt', delimiter=" ")
file = open(args.files_path+'/scan_output.txt','r')
lines = []
for line in file.xreadlines():
  lines.append(string.split(line.strip(),' '))
I = numpy.matrix(lines,float)
#I = I.transpose()
dims = I.shape
I0 = I[0,:].copy()
for i in range(0,dims[0]):
  I[i,:] = I[i,:]/I0

I = I + numpy.random.normal(loc=0,scale=float(args.noise_sigma),size=I.shape)

scipy.misc.imsave(args.files_path+'/linescan.tif', I)
scipy.misc.imsave(args.files_path+'/linescan.png', I)
