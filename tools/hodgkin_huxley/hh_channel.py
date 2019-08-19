import sys

file_out = sys.argv[1]
V_rest = float(sys.argv[2])
G = float(sys.argv[3])

exponents = []
equations = []
for i in range(4,len(sys.argv),2):
	equations.append(sys.argv[i])
	exponents.append(float(sys.argv[i+1]))

with open(file_out,'w') as f:
	f.write(str(V_rest) + "\n")
	f.write(str(G) + "\n")
	for i in range(0,len(exponents)):
		f.write(str(exponents[i])+";"+equations[i]+"\n")
f.close()