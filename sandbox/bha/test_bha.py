from scipy import *
from USTimeSeries import USTimeSeries
import bhacalc

# Specify the parameters for the underlying signal & sampling.
dt = 1/48000.
f = 984.3750	# The value near 1k for 1024 samples
f = 1007.8125	# The value near 1k for 2048 samples (1/2 way between for 1024)
#f = 1031.25
#f = 1017.	# A value between grid pts
#f = 1000.
A1 = 1000.
phi1 = 0.
A2 = 0.01*A1
phi2 = pi/2.
A3 = 0.0001*A1
phi3 = pi/4.
sig = 1.
N_s = 1024

# Simulate some real-valued data from a sinusoid with white noise.
times = arange(N_s)*dt
data = A1*cos(2*pi*f*times - phi1) + A2*cos(4*pi*f*times - phi2) + \
	A3*cos(6*pi*f*times - phi3) + random.normal(0.,sig,(N_s))

usts = USTimeSeries(data, dt)
usts.transform(8)
print 'f range: ', usts.f_range()

m = 3
prjxns = zeros((2*m))
fund, rl, im, p_fund = usts.f_info(f)
for i in range(1,m+1):
	fh = i*f
	fh, rl, im, p = usts.f_info(fh)
	print fh, i, 2*i-2, 2*i-1
	prjxns[2*i-2] = rl
	prjxns[2*i-1] = -im
print 'Prjxns: ', prjxns

gamma = 2*pi*fund*dt
gammas = arange(1,m+1)*gamma
print 'gammas: ', gammas
metric, L, J, amps, S = bhacalc.harmonicAnalysis(m, N_s, gamma, prjxns)
#print "metric: \n", metric
#print "Lower triangle:\n", L
print "Det: ", J
print 'Amps: ', amps
hamps = zeros((m))
for i in range(m):
	hamps[i] = sqrt( amps[2*i]**2 + amps[2*i+1]**2 )
print 'HAmps: ', hamps
A_fund = sqrt( amps[0]**2 + amps[1]**2 )
print 'Fundamental: ', A_fund
print '.............'

if 1:
	from bha import Harmonics
	harm = Harmonics(m, usts, sig)
	harm.analyze(f)
	print 'prjxns: ', harm.prjxns
	#print "metric: \n", harm.metric
	#print "Lower triangle:\n", harm.L
	print "Det: ", harm.jac
	print 'Suff. stat: ', harm.suf
	amps = harm.amplitudes()
	print 'Amps:  ', amps
	print 'Sigs:  ', harm.sigmas()
	print 'Ratio: ', amps/harm.sigmas()
else:
	usts.transform()
	for i in range(1,m+1):
		f, rl, im, p = usts.f_info(i*fund)
		print i, i*fund, sqrt( amps[2*i-2]**2 + amps[2*i-1]**2 )/A_fund, sqrt(p/p_fund)
	print 'Suff. stat: ', S


if 0:
	print '.............'
	metric, L, J, amps, S = bhacalc.multipleAnalysis(N_s, gammas, prjxns)
	#print "metric: \n", metric
	#print "Lower triangle:\n", L
	print "Det: ", J
	print 'Amps:  ', amps
	hamps = zeros((m))
	for i in range(m):
		hamps[i] = sqrt( amps[2*i]**2 + amps[2*i+1]**2 )
	print 'HAmps: ', hamps
	A_fund = sqrt( amps[0]**2 + amps[1]**2 )
	print 'Fundamental: ', A_fund
	usts.transform()
	for i in range(1,m+1):
		f, rl, im, p = usts.f_info(i*fund)
		print i, i*fund, sqrt( amps[2*i-2]**2 + amps[2*i-1]**2 )/A_fund, sqrt(p/p_fund)
	print 'Suff. stat: ', S
print "============="
