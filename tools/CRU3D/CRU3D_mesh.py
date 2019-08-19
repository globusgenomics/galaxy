#!/usr/bin/env python
"""
Generates CRU3D mesh.
"""

import sys
import os
import subprocess
import argparse

#Parse arguments
parser = argparse.ArgumentParser(description='Generate CRU3D mesh.')
parser.add_argument('--root_dir')
parser.add_argument('--html_file')
parser.add_argument('--files_path')
parser.add_argument('--size_x')
parser.add_argument('--size_y')
parser.add_argument('--size_z')
parser.add_argument('--bTT')
parser.add_argument('--tt_seg_circ')
parser.add_argument('--tt_seg_z')
parser.add_argument('--jsr_seg_circ')
parser.add_argument('--jsr_seg_z')
parser.add_argument('--jsr_r')
parser.add_argument('--subspace_r')
parser.add_argument('--trp_r')
parser.add_argument('--ryr')
parser.add_argument('--lcc')
args = parser.parse_args()

#Create output directory
#print "Creating directory: " + args.files_path
if not os.path.exists(args.files_path): 
    os.makedirs(args.files_path)

#Default parameters when TT not added
if not args.tt_seg_circ:
    args.tt_seg_circ = "23";
if not args.tt_seg_z:
    args.tt_seg_z = "47";
if not args.subspace_r:
    args.subspace_r = "0.015";

#Determine if dyad guide should be used
flag_ss_guide = "false"
if float(args.subspace_r) >= 0.040:
    flag_ss_guide = "true"

#Create XML
xml_content = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<mesh>
  <parameter>
    <symbol>meshset_title</symbol>
    <definition>Title of mesh parameter set</definition>
    <label>Title</label>
    <value>Mesh</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>flag_tt</symbol>
    <definition>Flag to include TT</definition>
    <label>Include TT</label>
    <value>%(bTT)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>dyad_r</symbol>
    <definition>Dyad height (distance from TT to jSR [um]</definition>
    <label>Dyad height (um)</label>
    <value>%(subspace_r)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>jsr_r</symbol>
    <definition>jSR thickness [um]</definition>
    <label>jSR Thickness [um]</label>
    <value>%(jsr_r)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>flag_dyad_guide</symbol>
    <definition>Flag to include guide points below jSR membrane (use when dyad height is at least 30nm)</definition>
    <label>Use dyad guide points</label>
    <value>%(flag_ss_guide)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>grid_x</symbol>
    <definition>Domain size along X axis in microns</definition>
    <label>Domain size (X) [um]</label>
    <value>%(size_x)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>n_jsr_seg_z</symbol>
    <definition>Number of 31nm sides along jSR height</definition>
    <label># of jSR segments (Z)</label>
    <value>%(jsr_seg_z)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>grid_y</symbol>
    <definition>Domain size along Y axis in microns</definition>
    <label>Domain size (Y) [um]</label>
    <value>%(size_y)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>grid_z</symbol>
    <definition>Domain size along Z axis in microns</definition>
    <label>Domain size (Z) [um]</label>
    <value>%(size_z)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>r_trpn</symbol>
    <definition>Radius from center of TT beyond which troponin and SERCa are present [um]</definition>
    <label>Troponin/SERCA radius [um]</label>
    <value>%(trp_r)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>description</symbol>
    <definition>Brief description of mesh parameter set</definition>
    <label>Description</label>
    <value></value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>n_jsr_seg</symbol>
    <definition>Number of 31nm sides along jSR arc</definition>
    <label># of jSR segments (circumference)</label>
    <value>%(jsr_seg_circ)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>cluster_lcc</symbol>
    <definition>LCC channel placement where 0's indicate empty space and 1's indicate channels on a grid with 31nm x 31nm spacing. Must be rectangular. This grid is centered in the dyadic space.</definition>
    <label>LCC (0/1-valued rectangle)</label>
    <value>%(lcc)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>n_tt_seg</symbol>
    <definition>Number of (~30nm) sides around the TT</definition>
    <label># of TT segments (circumference)</label>
    <value>%(tt_seg_circ)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>cluster_ryr</symbol>
    <definition>RyR channel placement where 0's indicate empty space and 1's indicate channels on a grid with 31nm x 31nm spacing. Must be rectangular. This grid is centered in the dyadic space.</definition>
    <label>RyR (0/1-valued rectangle)</label>
    <value>%(ryr)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>n_tt_seg_z</symbol>
    <definition>Number of 31nm sides along length of the TT</definition>
    <label># of TT segments (Z)</label>
    <value>%(tt_seg_z)s</value>
    <listable>false</listable>
  </parameter>
</mesh> """ % {'bTT':args.bTT, 
    'subspace_r':args.subspace_r, 
    'jsr_r':args.jsr_r,
    'jsr_seg_z':args.jsr_seg_z,
    'jsr_seg_circ':args.jsr_seg_circ,
    'tt_seg_z':args.tt_seg_z,
    'tt_seg_circ':args.tt_seg_circ,
    'size_x':args.size_x,
    'size_y':args.size_y,
    'size_z':args.size_z,
    'ryr':args.ryr,
    'lcc':args.lcc,
    'trp_r':args.trp_r,
    'flag_ss_guide':flag_ss_guide}

#Write xml file
xml_path = args.files_path+"/Mesh_Parameters.xml"
f = open(xml_path,'w')
f.write(xml_content)
f.close()

#Run mesh program
base_dir = args.root_dir + '/tools/CRU3D/model/mesh';
proc = subprocess.Popen([base_dir+'/mesh','-param',xml_path,'-out',args.files_path+"/mesh"],cwd=base_dir,stdout=subprocess.PIPE)

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
f.write('<a href="log.txt">Mesh Log</a><br>')
f.write('<a href="Mesh_Parameters.xml">Mesh Parameters</a><br>')
f.write('<a href="mesh.node">Node file</a><br>')
f.write('<a href="mesh.face">Face file</a><br>')
f.write('<a href="mesh.neigh">Neighbor file</a><br>')
f.write('<a href="mesh.ele">Ele file</a><br>')
f.write('<a href="mesh_ryr.txt">RyR file</a><br>')
f.write('<a href="mesh_lcc.txt">LCC file</a><br>')
f.write('<a href="mesh_face_properties.txt">Face properties file</a><br>')
f.write('<a href="mesh_tet_properties.txt">Tetrahedra properties file</a><br>')
f.close()
