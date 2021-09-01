import string
from math import *
import fkepler

Kepler_c = 0.
Kepler_s = 0.

def Kepler_setup(tau, e, T0, t):
	global Kepler_c, Kepler_s
	fkepler.setup(tau, e, T0)
	Kepler_c, Kepler_s = fkepler.t2TA(t)

def Kepler0(tau, e, T0, t):
	return 1.

def Kepler1(tau, e, T0, t):
	return e + Kepler_c

def Kepler2(tau, e, T0, t):
	return Kepler_s

def v_rad(t, K, w, tau, e, T0):
	Kepler_setup(tau, e, T0, t)
	return -20.*Kepler0(tau, e, T0, t) + K*cos(w)*Kepler1(tau, e, T0, t) -\
	 K*sin(w)*Kepler2(tau, e, T0, t)

def v_loop():
	K = 30.
	w = 0.4*pi
	tau = 17.1
	e = 0.6
	T0 = 37.5

	ofname = 'kepler.dat'
	ofile = open(ofname, 'w')

	lo = 30.
	hi = 60.
	n = 6001
	dt = (hi-lo)/(n-1)
	for i in range(n):
		t = lo + i*dt
		line = "%-6.4f %-6.4f\n" % (t, v_rad(t, K, w, tau, e, T0))
		ofile.write(line)

	ofile.close()

#import profile
#profile.run("v_loop()")
#profile.run("v_loop()")
