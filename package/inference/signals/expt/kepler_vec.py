import string
from math import *
from scipy import zeros, ones, Float, sum
from scipy import random as r
import fkepler

fkepler.settol(1.e-3)

def Kepler_setup(tau, e, T0, times):
	global cs_vals
	fkepler.setup(tau, e, T0)
	return fkepler.vt2TA(times)

def Kepler0(tau, e, T0, times, cs_vals):
	return ones(len(times), Float)

def Kepler1(tau, e, T0, times, cs_vals):
	return e + cs_vals[:,0]

def Kepler2(tau, e, T0, times, cs_vals):
	return cs_vals[:,1]

def phaseCurve(ofile, tau, e, T0, amps, data, n, sigma=None):
	"""Write data for plotting a velocity curve against observations.
	The data are written as phases in [0,1), and the curve is written
	as n phases in [-.5, 1.5] and the corresponding velocities.
	Return chi**2 for the curve."""

	# First, make the curve.
	dt = 2./(n-1)
	phases = zeros(n, Float)
	for i in range(n):
		phases[i] = -.5 + i*dt
	times = tau*phases
	cs = Kepler_setup(tau, e, T0, times)
	K0 = Kepler0(tau, e, T0, times, cs)
	K1 = Kepler1(tau, e, T0, times, cs)
	K2 = Kepler2(tau, e, T0, times, cs)
	v_rad = amps[0]*K0 + amps[1]*K1 + amps[2]*K2
	for i in range(n):
		line = 'vcurve  %-6.4f  %-6.4f\n' % (phases[i], v_rad[i])
		ofile.write(line)
	ofile.write('\n')

	# Now convert the data, gathering times for chi**2 along the way.
	times = zeros(len(data), Float)
	for i in range(len(data)):
		times[i] = data[i,0]
		np, phase = divmod(data[i,0], tau)
		phase = phase / tau
		# print i, tau, data[i,0], np, phase
		v = data[i,1]
		if sigma:
			sig = sigma
		else:
			sig = data[i,2]
		line = 'datum  %-6.4f  %-6.4f  %-6.4f\n' % (phase, v, sig)
		ofile.write(line)
	ofile.write('\n')

	# Calculate chi**2.
	cs = Kepler_setup(tau, e, T0, times)
	K0 = Kepler0(tau, e, T0, times, cs)
	K1 = Kepler1(tau, e, T0, times, cs)
	K2 = Kepler2(tau, e, T0, times, cs)
	v_rad = amps[0]*K0 + amps[1]*K1 + amps[2]*K2
	resid = data[:,1] - v_rad
	if sigma:
		resid = resid / sigma
	else:
		resid = resid / data[:,2]
	chisq = sum(resid*resid)
	line = 'chisq  %-6.4f\n' % chisq
	ofile.write(line)
	return chisq



def v_loop():
	K = 30.
	w = 0.4*pi
	tau = 17.1
	e = 0.6
	T0 = 37.5

	ofname = 'kepler_vec.dat'
	ofile = open(ofname, 'w')

	lo = 30.
	hi = 60.
	n = 6001
	dt = (hi-lo)/(n-1)
	times = zeros(n, Float)
	for i in range(n):
		times[i] = lo + i*dt
	cs = Kepler_setup(tau, e, T0, times)
	K0 = Kepler0(tau, e, T0, times, cs)
	K1 = Kepler1(tau, e, T0, times, cs)
	K2 = Kepler2(tau, e, T0, times, cs)
	v_rad = -20.*K0 + K*cos(w)*K1 - K*sin(w)*K2
	for i in range(n):
		line = "%-6.4f %-6.4f\n" % (times[i], v_rad[i])
		ofile.write(line)

	ofile.close()
	print "Done!"

#import profile
#profile.run("v_loop()")
#profile.run("v_loop()")

def v_rad(K, w, tau, e, T0, times):
	cs = Kepler_setup(tau, e, T0, times)
	K1 = Kepler1(tau, e, T0, times, cs)
	K2 = Kepler2(tau, e, T0, times, cs)
	return K*cos(w)*K1 - K*sin(w)*K2


def keplerSim(tau, e, T0, K, w, sig, tlo, thi, n):
	dt = (thi-tlo)/(n-1)
	data = zeros((n,2), Float)
	data[:,0] = r.uniform(tlo, thi, (n))
#	for i in range(n):
#            data[i,0] = tlo + i*dt
        data[:,1] = v_rad(K, w, tau, e, T0, data[:,0])+r.normal(0.,sig,(n))
        print "Created data."
	return data

if 1:
	import ioutils
	from scipy import array
	v0 = 15.
	K = 10.
	w = 0.4*pi
	tau = 17.1
	e = 0.6
	T0 = 37.5
	sig = 5.
	if 1:
		data = keplerSim(tau, e, T0, K, w, sig, 30., 60., 21)
		data[:,1] = data[:,1] + v0
		print repr(data)
	else:
		import kdata
		data = kdata.data10r2
	ofile = ioutils.openText('kepler_vec.dat')
	amps = (v0, K*cos(w), -K*sin(w))
	print 'amps: ', amps
	print 'chisq = ', phaseCurve(ofile, tau, e, T0, amps, data, 200, sigma=5.)
	# amps = [  3.13680476,  68.23499664, -60.6215832 ]
	# print phaseCurve(ofile, 16.53, .9, 4.314, amps, data, 200, sigma=5.)
	ofile.close()
