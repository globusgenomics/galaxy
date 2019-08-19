from scipy import *
from scipy.integrate import odeint
import numpy
import parser
import sys

file_out = sys.argv[1]
t_final = float(sys.argv[2])
t_step = float(sys.argv[3])
t_delay = float(sys.argv[4])
t_stim = float(sys.argv[5])
C_m = float(sys.argv[6])
I = float(sys.argv[7])

N_Channels = len(sys.argv)-8

V_rest = []
equations = []
exponents = []
conductances = []
N_Equations = 0
for i in range(0,N_Channels):
	f = open(sys.argv[8+i],'r')
	lines = f.readlines()
	V_rest.append(double(lines[0]))
	conductances.append(double(lines[1]))
	exp_i = []
	eq_i = []
	for line in lines[2:]:
		tok = line.split(';')
		exp_i.append(float(tok[0]))
		eq_i.append(parser.expr(tok[1].strip()).compile())
		N_Equations = N_Equations + 1
	equations.append(eq_i)
	exponents.append(exp_i)
	f.close()


def hodgkinHuxley(yy,t, p):

	#Name the variables
	V = yy[0]
	
	#Name the parameters
	C_m = p[0]
	I = p[1]
	equations = p[2]
	exponents = p[3]
	V_rest = p[4]
	conductances = p[5]
	t_delay = p[6]
	
	V_dot = 0
	if t >= t_delay and t < t_delay+t_stim:
		V_dot = -I
	
	dx_idx = 1
	for i in range(0,len(equations)):
		x = conductances[i]*(V-V_rest[i])
		for j in range(0,len(equations[i])):
			x = x * pow(yy[dx_idx],exponents[i][j])
			dx_idx = dx_idx + 1
		V_dot = V_dot + x
	V_dot = (-1/C_m)*V_dot
	
	y_dot = [V_dot]
	
	dx_idx = 1
	for i in range(0,len(equations)):
		for j in range(0,len(equations[i])):
			x = yy[dx_idx]
			y_dot.append(eval(equations[i][j]))
			dx_idx = dx_idx + 1
	
	return y_dot


#Define initial conditions for the state variables
y0 = [0]
for i in range(0,N_Equations):
	y0.append(0)

#Define time interval and spacing for the solutions
t = arange(0,t_final,t_step)

#Pack the parameters in a single vector
p = [C_m,I,equations,exponents,V_rest,conductances,t_delay]

#Call the integrator 
y = odeint(hodgkinHuxley, y0, t, args=(p,))

V = y[:,0]
data = hstack([mat(t).H, mat(V).H])

idx = 1
for i in range(0,len(equations)):
  x = conductances[i]*(V-V_rest[i])
  for j in range(0,len(equations[i])):
    x = multiply(x, power(y[:,idx],exponents[i][j]))
    idx = idx + 1
  data = hstack([data, mat(x).H])

numpy.savetxt(file_out,data,delimiter="\t")
