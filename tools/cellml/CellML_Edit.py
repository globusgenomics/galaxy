#!/usr/bin/env python
"""
Runs CellML simulation.
"""

import sys
import cgrspy.bootstrap

output_file = sys.argv[1]
input_file = sys.argv[2]

cgrspy.bootstrap.loadGenericModule('cgrs_cellml')
cellmlBootstrap = cgrspy.bootstrap.fetch('CreateCellMLBootstrap')
model = cellmlBootstrap.modelLoader.loadFromURL(input_file)

#Modify parameters:
for i in range(3,len(sys.argv),3):
	c = model.allComponents.getComponent(sys.argv[i])
	if c:
		v = c.variables.getVariable(sys.argv[i+1])
		if v:
			if v.initialValue:
				if v.publicInterface != "in":
					v.initialValue = sys.argv[i+2]
				else:
					sys.stderr.write("Error: Variable %s in component %s has public interface of type \"in\".\n" % (sys.argv[i+1],sys.argv[i]))
			else:
				sys.stderr.write("Error: Variable %s in component %s does not have an initial value attribute.\n" % (sys.argv[i+1],sys.argv[i]))
		else:
			sys.stderr.write("Error: Could not find model variable %s in component %s\n" % (sys.argv[i+1],sys.argv[i]))
	else:
		sys.stderr.write("Error: Could not find model component %s\n" % sys.argv[i])

#Write modified model to XML
f = open(output_file,'w')
f.write(model.serialisedText)
f.close()
