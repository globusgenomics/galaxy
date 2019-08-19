#!/usr/bin/env python
"""
Runs CRU3D simulation.
"""

import sys
import os
import subprocess
import argparse

#Parse arguments
parser = argparse.ArgumentParser(description='Run CRU3D simulation.')
parser.add_argument('--root_dir')
parser.add_argument('--html_file')
parser.add_argument('--files_path')
parser.add_argument('--mesh')
parser.add_argument('--protocol')
parser.add_argument('--ryr_init')
parser.add_argument('--ryr_init_time')
parser.add_argument('--t_final')
parser.add_argument('--grid_interval')
parser.add_argument('--num_proc')
parser.add_argument('--num_sims')
parser.add_argument('--c_min')
parser.add_argument('--ryr_open_max')
parser.add_argument('--output_ryr_states')
parser.add_argument('--t_clamp')
parser.add_argument('--v_clamp')
parser.add_argument('--t_step')
parser.add_argument('--states_interval')
parser.add_argument('--A_P')
parser.add_argument('--B_TOT_ATP')
parser.add_argument('--B_TOT_CMDN')
parser.add_argument('--B_TOT_CSQN')
parser.add_argument('--B_TOT_DYE')
parser.add_argument('--B_TOT_SL')
parser.add_argument('--B_TOT_TRPN')
parser.add_argument('--CA_I')
parser.add_argument('--CA_O')
parser.add_argument('--CA_SR')
parser.add_argument('--D_ATP')
parser.add_argument('--D_CA')
parser.add_argument('--D_CAJSR')
parser.add_argument('--D_CMDN')
parser.add_argument('--D_CSQN')
parser.add_argument('--D_DYE')
parser.add_argument('--K_D_I')
parser.add_argument('--K_D_SL')
parser.add_argument('--K_D_SR')
parser.add_argument('--K_OFF_ATP')
parser.add_argument('--K_OFF_CMDN')
parser.add_argument('--K_OFF_CSQN')
parser.add_argument('--K_OFF_DYE')
parser.add_argument('--K_OFF_TRPN')
parser.add_argument('--K_ON_ATP')
parser.add_argument('--K_ON_CMDN')
parser.add_argument('--K_ON_CSQN')
parser.add_argument('--K_ON_DYE')
parser.add_argument('--K_ON_TRPN')
parser.add_argument('--RYR_ETA')
parser.add_argument('--RYR_K_MINUS')
parser.add_argument('--RYR_K_PLUS')
parser.add_argument('--RYR_PHI_KD')
parser.add_argument('--RYR_PHI_N')
parser.add_argument('--V_M')
parser.add_argument('--V_REFILL')
parser.add_argument('--V_RYR')
args = parser.parse_args()

#Create output directory
#print "Creating directory: " + args.files_path
if not os.path.exists(args.files_path): 
    os.makedirs(args.files_path)

#Determine extra parameters
output_grid = "false"
if args.protocol=="spark":
  output_grid = "true"

#Create XML parameters file the old-fashioned way
xml_content = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<simulation>
  <parameter>
    <symbol>V_REFILL</symbol>
    <definition>jSR refill rate [1/ms]</definition>
    <label/>
    <value>%(V_REFILL)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>ryr_open_max</symbol>
    <definition>Maximum number of RyRs that can be open before terminating simulation</definition>
    <label>Max # open RyRs</label>
    <value>%(ryr_open_max)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>K_D_SL</symbol>
    <definition>Sarcolemma binding site affinity [uM]</definition>
    <label/>
    <value>%(K_D_SL)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>K_ON_CSQN</symbol>
    <definition>Calsequestrin on rate [1/uM-ms]</definition>
    <label/>
    <value>%(K_ON_CSQN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>K_D_SR</symbol>
    <definition>SERCA K_D_SR constant [uM]</definition>
    <label/>
    <value>%(K_D_SR)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>B_TOT_ATP</symbol>
    <definition>ATP total concentration [uM]</definition>
    <label/>
    <value>%(B_TOT_ATP)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>K_ON_ATP</symbol>
    <definition>ATP on rate [1/uM-ms]</definition>
    <label/>
    <value>%(K_ON_ATP)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>RYR_K_MINUS</symbol>
    <definition>RyR close rate constant [1/ms]</definition>
    <label/>
    <value>%(RYR_K_MINUS)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>description</symbol>
    <definition>Brief description of this parameter set</definition>
    <label>Description</label>
    <value></value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>D_ATP</symbol>
    <definition>ATP diffusion coefficient [um^2/ms]</definition>
    <label/>
    <value>%(D_ATP)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>K_OFF_CMDN</symbol>
    <definition>Calmodulin off rate [1/ms]</definition>
    <label/>
    <value>%(K_OFF_CMDN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>RYR_A_STAR</symbol>
    <definition>RyR coupling strength</definition>
    <label/>
    <value>0</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>flag_output_ryr_states</symbol>
    <definition>Flag to output RyR states</definition>
    <label>Output RyR states</label>
    <value>%(output_ryr_states)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>t_final</symbol>
    <definition>Simulation end time (ms)</definition>
    <label>Final time (ms)</label>
    <value>%(t_final)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>K_ON_TRPN</symbol>
    <definition>Troponin on rate [1/uM-ms]</definition>
    <label/>
    <value>%(K_ON_TRPN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>seed_A</symbol>
    <definition>PRNG seed multiplied with process index</definition>
    <label>Device seed</label>
    <value>12783</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>K_D_I</symbol>
    <definition>SERCA K_D_i constant [uM]</definition>
    <label/>
    <value>%(K_D_I)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>seed_B</symbol>
    <definition>PRNG seed multiplied with simulation index</definition>
    <label>Ensemble seed</label>
    <value>83</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>V_RYR</symbol>
    <definition>RyR release rate [1/ms]</definition>
    <label/>
    <value>%(V_RYR)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>flag_nogating</symbol>
    <definition>Flag to disable channel gating</definition>
    <label>No channel gating</label>
    <value>false</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>RYR_K_PLUS</symbol>
    <definition>RyR open rate constant [1/ms]</definition>
    <label/>
    <value>%(RYR_K_PLUS)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>ryr_open_time</symbol>
    <definition>Amount of time to hold first RyR open for (ms) (optional)</definition>
    <label>Initial open RyR time (ms)</label>
    <value>%(ryr_init_time)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>flag_ffwd</symbol>
    <definition>Flag to enable simulation acceleration until first channel opening (assumes system initially at steady-state)</definition>
    <label>Skip to first opening</label>
    <value>false</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>B_TOT_CMDN</symbol>
    <definition>Calmodulin total concentration [uM]</definition>
    <label/>
    <value>%(B_TOT_CMDN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>K_OFF_CSQN</symbol>
    <definition>Calsequestrin off rate [1/ms]</definition>
    <label/>
    <value>%(K_OFF_CSQN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>K_ON_DYE</symbol>
    <definition>Fluo-4N on rate [1/uM-ms]</definition>
    <label/>
    <value>%(K_ON_DYE)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>B_TOT_SL</symbol>
    <definition>Sarcolemma binding site density [umol/um^2]</definition>
    <label/>
    <value>%(B_TOT_SL)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>K_OFF_ATP</symbol>
    <definition>ATP off rate [1/ms]</definition>
    <label/>
    <value>%(K_OFF_ATP)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>RYR_EPS_OO</symbol>
    <definition>RyR open-open coupling interaction</definition>
    <label/>
    <value>-0.85</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>K_OFF_DYE</symbol>
    <definition>Fluo-4N off rate [1/ms]</definition>
    <label/>
    <value>%(K_OFF_DYE)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>states_interval</symbol>
    <definition>Time between state outputs (ms)</definition>
    <label>State output interval (ms)</label>
    <value>%(states_interval)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>sims_per_proc</symbol>
    <definition>Number of simulations per processor</definition>
    <label>Simulations per processor</label>
    <value>%(num_sims)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>flag_output_grid</symbol>
    <definition>Flag to output mesh state</definition>
    <label>Output mesh state</label>
    <value>%(output_grid)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>K_OFF_TRPN</symbol>
    <definition>Troponin off rate [1/ms]</definition>
    <label/>
    <value>%(K_OFF_TRPN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>t_clamp</symbol>
    <definition>Clamp holding time (ms)</definition>
    <label>Clamp duration (ms)</label>
    <value>%(t_clamp)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>K_ON_CMDN</symbol>
    <definition>Calmodulin on rate [1/uM-ms]</definition>
    <label/>
    <value>%(K_ON_CMDN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>D_CMDN</symbol>
    <definition>Calmodulin diffusion coefficient [um^2/ms]</definition>
    <label/>
    <value>%(D_CMDN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>CA_SR</symbol>
    <definition>Calcium concentration (SR) [uM]</definition>
    <label/>
    <value>%(CA_SR)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>V_M</symbol>
    <definition>Membrane potential (rest) [mV]</definition>
    <label/>
    <value>%(V_M)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>t_step</symbol>
    <definition>Euler time step size (ms)</definition>
    <label>Time step (ns)</label>
    <value>%(t_step)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>D_CAJSR</symbol>
    <definition>Calcium diffusion coefficient (jSR) [um^2/ms]</definition>
    <label/>
    <value>%(D_CAJSR)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>D_DYE</symbol>
    <definition>Fluo-4N diffusion coefficient [um^2/ms]</definition>
    <label/>
    <value>%(D_DYE)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>RYR_EPS_CC</symbol>
    <definition>RyR closed-closed coupling interaction</definition>
    <label/>
    <value>-0.92</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>flag_noflux</symbol>
    <definition>Flag to enable no-flux boundary conditions</definition>
    <label>No flux boundaries</label>
    <value>true</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>D_CA</symbol>
    <definition>Calcium diffusion coefficient [um^2/ms]</definition>
    <label/>
    <value>%(D_CA)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>B_TOT_CSQN</symbol>
    <definition>Calsequestrin total concentration (jSR) [uM]</definition>
    <label/>
    <value>%(B_TOT_CSQN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>paramset_title</symbol>
    <definition>Title for parameter set</definition>
    <label>Title</label>
    <value>CRU3D</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>A_P</symbol>
    <definition>SERCA A_P constant [uM]</definition>
    <label/>
    <value>%(A_P)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>v_clamp</symbol>
    <definition>Clamp membrane potential (mV)</definition>
    <label>Clamp potential (mV)</label>
    <value>%(v_clamp)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>CA_I</symbol>
    <definition>Calcium concentration (cytosol) [uM]</definition>
    <label/>
    <value>%(CA_I)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>protocol</symbol>
    <definition>Simulation protocol type (spark/fidelity/gain)</definition>
    <label>Select protocol</label>
    <value>%(protocol)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>RYR_PHI_N</symbol>
    <definition>RyR luminal dependence power</definition>
    <label/>
    <value>%(RYR_PHI_N)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>num_procs</symbol>
    <definition>Number of processors to use</definition>
    <label># of processors</label>
    <value>%(num_proc)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>seed_base</symbol>
    <definition>Base PRNG seed</definition>
    <label>Base seed</label>
    <value>38376427</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>c_min</symbol>
    <definition>Minimum calcium concentration before terminating simulation (when all RyRs closed)</definition>
    <label>Min Ca2+ concentration (uM)</label>
    <value>%(c_min)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>CA_O</symbol>
    <definition>Calcium concentration (extracellular) [uM]</definition>
    <label/>
    <value>%(CA_O)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>ryr_open_init</symbol>
    <definition>Index of RyR to open initially (optional)</definition>
    <label>Initially-open RyR</label>
    <value>%(ryr_init)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>flag_use_clock</symbol>
    <definition>Use clock to initialize the PRNG seed</definition>
    <label>Use clock seed</label>
    <value>true</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>grid_interval</symbol>
    <definition>Time between mesh state outputs (ms)</definition>
    <label>Mesh state output interval (ms)</label>
    <value>%(grid_interval)s</value>
    <listable>false</listable>
  </parameter>
  <parameter>
    <symbol>RYR_PHI_KD</symbol>
    <definition>RyR luminal dependence constant [uM]</definition>
    <label/>
    <value>%(RYR_PHI_KD)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>B_TOT_DYE</symbol>
    <definition>Fluo-4N total concentration [uM]</definition>
    <label/>
    <value>%(B_TOT_DYE)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>B_TOT_TRPN</symbol>
    <definition>Troponin total concentration [uM]</definition>
    <label/>
    <value>%(B_TOT_TRPN)s</value>
    <listable>true</listable>
  </parameter>
  <parameter>
    <symbol>RYR_ETA</symbol>
    <definition>RyR Hill coefficient</definition>
    <label/>
    <value>%(RYR_ETA)s</value>
    <listable>true</listable>
  </parameter>
</simulation>
""" % {'ryr_init':args.ryr_init,
       'ryr_init_time':args.ryr_init_time,
       't_final':args.t_final,
       'grid_interval':args.grid_interval,
       'num_proc':args.num_proc,
       'num_sims':args.num_sims,
       'c_min':args.c_min,
       'ryr_open_max':args.ryr_open_max,
       'output_ryr_states':args.output_ryr_states,
       't_clamp':args.t_clamp,
       'v_clamp':args.v_clamp,
       't_step':args.t_step,
       'states_interval':args.states_interval,
       'A_P':args.A_P,
       'B_TOT_ATP':args.B_TOT_ATP,
       'B_TOT_CMDN':args.B_TOT_CMDN,
       'B_TOT_CSQN':args.B_TOT_CSQN,
       'B_TOT_DYE':args.B_TOT_DYE,
       'B_TOT_SL':args.B_TOT_SL,
       'B_TOT_TRPN':args.B_TOT_TRPN,
       'CA_I':args.CA_I,
       'CA_O':args.CA_O,
       'CA_SR':args.CA_SR,
       'D_ATP':args.D_ATP,
       'D_CA':args.D_CA,
       'D_CAJSR':args.D_CAJSR,
       'D_CMDN':args.D_CMDN,
       'D_DYE':args.D_DYE,
       'K_D_I':args.K_D_I,
       'K_D_SL':args.K_D_SL,
       'K_D_SR':args.K_D_SR,
       'K_OFF_ATP':args.K_OFF_ATP,
       'K_OFF_CMDN':args.K_OFF_CMDN,
       'K_OFF_CSQN':args.K_OFF_CSQN,
       'K_OFF_DYE':args.K_OFF_DYE,
       'K_OFF_TRPN':args.K_OFF_TRPN,
       'K_ON_ATP':args.K_ON_ATP,
       'K_ON_CMDN':args.K_ON_CMDN,
       'K_ON_CSQN':args.K_ON_CSQN,
       'K_ON_DYE':args.K_ON_DYE,
       'K_ON_TRPN':args.K_ON_TRPN,
       'RYR_ETA':args.RYR_ETA,
       'RYR_K_MINUS':args.RYR_K_MINUS,
       'RYR_K_PLUS':args.RYR_K_PLUS,
       'RYR_PHI_KD':args.RYR_PHI_KD,
       'RYR_PHI_N':args.RYR_PHI_N,
       'V_M':args.V_M,
       'V_REFILL':args.V_REFILL,
       'V_RYR':args.V_RYR,
       'output_grid':output_grid,
       'protocol':args.protocol}

#Write xml file. This will be read in by the cru3d program
xml_path = args.files_path+"/Simulation_Parameters.xml"
f = open(xml_path,'w')
f.write(xml_content)
f.close()

#Run mesh program with a call to mpirun
base_dir = args.root_dir + '/tools/CRU3D/model';
proc = subprocess.Popen(['mpirun','-np',args.num_proc,base_dir+'/cru3d','-param',xml_path,'-out',args.files_path,'-mesh',args.mesh+"/mesh"],cwd=base_dir,stdout=subprocess.PIPE)

#Send output to file
output = proc.stdout.read()
f = open(args.files_path+"/log.txt",'w')
f.write(output)
f.close()

#Create HTML file
f = open(args.html_file,'w')
#f.write(xml_content.replace("\n","<br>"))
f.write('<a href="log.txt">Simulation Log</a><br>')
f.write('<a href="Simulation_Parameters.xml">Simulation Parameters</a><br>')
f.close()
