#!/usr/bin/env python
"""
Runs CellML simulation.
"""

import sys
import subprocess
import xml.dom.minidom as MD

output_file = sys.argv[1]
output_cellsim = sys.argv[2]
input_cellml = sys.argv[3]
t0 = sys.argv[4]
tf = sys.argv[5]
maxStep = sys.argv[6]
tabulationStep = sys.argv[7]

#Add simulation tag for CellMLSimulator
doc = MD.parse(input_cellml)
sim_ele = doc.createElementNS('http://cellml.sourceforge.net/csim/simulation/0.1#','simulation')
sim_ele.setAttribute('xmlns','http://cellml.sourceforge.net/csim/simulation/0.1#')
sim_ele.setAttribute('id','simulation1')
doc.documentElement.appendChild(sim_ele)

bound_var = doc.createElement('boundVariable')
bound_var.setAttribute('start',t0)
bound_var.setAttribute('end',tf)
bound_var.setAttribute('maxStep',maxStep)
bound_var.setAttribute('tabulationStep',tabulationStep)
sim_ele.appendChild(bound_var)
	

#Modify parameters:
col = 1
for i in range(8,len(sys.argv),2):
	output_var = doc.createElement('outputVariable')
	output_var.setAttribute('component',sys.argv[i])
	output_var.setAttribute('variable',sys.argv[i+1])
	output_var.setAttribute('column',str(col))
	sim_ele.appendChild(output_var)
	col = col + 1

#Write modified model to XML for CellMLSimulator
xml_out = output_cellsim
f = open(xml_out,'w')
f.write(doc.toprettyxml())
f.close()

#Run model
proc = subprocess.Popen(['csim',xml_out], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

#Send output to file
output = proc.stdout.read()
#f = open(output_dir+"/raw_output.txt",'w')
#f.write(output)
#f.close()

#Check for errors
if output.find("ERROR")!=-1:
	sys.stderr.write("An error occurred. See raw output for details:\n")
	sys.stderr.write(output)

#Parse output variables
str_start = 'Running the simulation: simulation1';
i_start = output.find(str_start)
i_end = output.find('Wall clock time')
if i_start==-1 or i_end==-1:
	sys.stderr.write("Error: Could not parse output variables. See raw output for details:\n")
	sys.stderr.write(output)
out_vars = output[(i_start+len(str_start)):i_end]

out_final = ""
for line in out_vars.split('\n'):
	out_final = out_final + line.strip() + "\n"

#Send output to file
f = open(output_file,'w')
f.write(out_final)
f.close()

#Write HTML file
#f = open(html_file,'w')
#f.write('<a href="model.cellml.txt">CellML Model</a><br>')
#f.write('<a href="raw_output.txt">Raw Output</a><br>')
#f.write('<a href="vars_output.txt">Variable Output</a><br>')
#f.close()
