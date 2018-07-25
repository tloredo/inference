from numpy import *
import bhacalc

__all__ = ['Harmonics', 'Multiplet']

def qinterp(x, xvals, yvals):
	"""Quadratically interpolate f(x) from 3 values yvals[i] = f(xvals[i])."""
	if len(xvals) != 3 or len(yvals) !=3:
		raise ValueError, 'Need exactly 3 values!'
	if x < xvals[0] or x > xvals[2]:
		raise ValueError, 'x out of range---dangerous to extrapolate!'
	d = x - xvals
	x01 = xvals[0] - xvals[1]
	x02 = xvals[0] - xvals[2]
	x12 = xvals[1] - xvals[2]
	return yvals[0]*d[1]*d[2]/(x01*x02) - yvals[1]*d[0]*d[2]/(x01*x12) \
		+ yvals[2]*d[0]*d[1]/(x02*x12)

def test_qinterp():
	def f(x): return 3. + 2*x - 7*x**2
	xvals = array([.4, 3., 8.])
	yvals = []
	for x in xvals:
		yvals.append(f(x))
	yvals = array(yvals)
	xi = .4 + .5*arange(15)
	for x in xi:
		print x, f(x), qinterp(x, xvals, yvals)
	print '===== test_qinterp ====='

#test_qinterp()


################################################################################
class Harmonics:

	"""Analyze using a model consisting of a fundamental & integer harmonics."""

#...............................................................................
	def __init__(self, nharm, usts, sigma):
		"""Initialize with number of harmonics (incl. fundamental) & a USTimeSeries 
		instance.  The time series should already be transformed (preferably 
		significantly oversampled)."""
		if not usts.did_fft: raise ValueError, 'Must transform time series first!'
		self.nharm = nharm
		self.usts = usts
		self.var = sigma**2
		self.f = None

#...............................................................................
	def analyze(self, f):
		"""Perform an analysis using as a fundamental the FFT frequency nearest
		to f."""
		self.f_in = f
		self.f, rl, im, p = self.usts.f_info(f)
		self.prjxns = zeros((2*self.nharm))
		self.prjxns[0], self.prjxns[1] = rl, -im
		for i in range(2,self.nharm+1):
			fh = i*self.f
			fh, rl, im, p = self.usts.f_info(fh)
			print 'harm: ', i, i*self.f, fh
			self.prjxns[2*i-2], self.prjxns[2*i-1] = rl, -im
		gamma = 2*pi*self.f*self.usts.dt
		self.metric, self.L, self.jac, self.amps, self.suf = \
			bhacalc.harmonicAnalysis(self.nharm, self.usts.n, gamma, self.prjxns)
		#self.Q = self.dsqr - self.suf
		self.jac = 1. / self.jac
		#return (self.suf, self.Jac)
		self.hamps = zeros((self.nharm))
		for i in range(self.nharm):
			self.hamps[i] = sqrt(self.amps[2*i]**2 + self.amps[2*i+1]**2)

#...............................................................................
	def amplitudes(self):
		"""Return the amplitudes of the harmonics."""
		
		if self.f == None: raise RuntimeError, 'Need to call analyze first!'
		return self.hamps

#...............................................................................
	def covar(self):
		"""Return the covariance matrix for the best-fit cos & sin amplitudes."""

#...	If we haven't already done the metric calculations, complain!
		if self.f == None:
			raise RuntimeError, 'Must call analyze first!'

		return bhacalc.covar(self.L) * self.var

#...............................................................................
	def sigmas(self):
		"""Return the standard deviations for the best-fit harmonic amplitudes."""

#...	If we haven't already done the metric calculations, complain!
		if self.f == None:
			raise RuntimeError, 'Must call analyze first!'

		covars = bhacalc.covar(self.L) * self.var
		sigs = zeros((self.nharm))
		for i in range(self.nharm):
			i1, i2 = 2*i, 2*i+1
			sigs[i] = (self.amps[i1]/self.hamps[i])**2 * covars[i1][i1] + \
					  (self.amps[i2]/self.hamps[i])**2 * covars[i2][i2]
		return sqrt(sigs)

################################################################################
class Multiplet:

	"""Analyze using a model consisting of a multiplet of sinusoids."""

#...............................................................................
	def __init__(self, usts, sigma=None, interp=None, exact=None):
		"""Initialize with number of harmonics (incl. fundamental) & a USTimeSeries 
		instance.  The time series should already be transformed (preferably 
		significantly oversampled)."""
		if not usts.did_fft: raise ValueError, 'Must transform time series first!'
		self.usts = usts
		self.dsqr = sum(usts.samples**2)
		self.n = usts.n
		if sigma:
			self.var = sigma**2
		else:
			self.var = None
		self.f_in = None
		self.interp = interp
		self.exact = exact

#...............................................................................
	def analyze(self, flist):
		"""Perform an analysis using as a fundamental the FFT frequency nearest
		to f."""
		self.nf = len(flist)
		self.f_in = flist
		self.dof = self.n - self.nf

		# Get projections on the grid or by interpolating; adjust flist if
		# fixed to the grid.
		gammas = array(flist)
		self.prjxns = zeros((2*self.nf))
		if self.exact:
			for i in range(self.nf):
				self.prjxns[2*i] = self.usts.realDTFT(gammas[i])
				self.prjxns[2*i+1] = -self.usts.imagDTFT(gammas[i])
		elif self.interp:
			for i in range(self.nf):
				fvals, rvals, ivals = self.usts.bracket3(gammas[i])
				ivals = - ivals
				self.prjxns[2*i] = qinterp(gammas[i], fvals, rvals)
				self.prjxns[2*i+1] = qinterp(gammas[i], fvals, ivals)
				print i, gammas[i], fvals
				print '  --> ', rvals, self.prjxns[2*i]
				print '  --> ', ivals, self.prjxns[2*i+1]
		else:
			for i in range(self.nf):
				f, rl, im, p = self.usts.f_info(gammas[i])
				gammas[i] = f	# Reset gamma to the grid value
				self.prjxns[2*i], self.prjxns[2*i+1] = rl, -im
		self.flist = gammas
		gammas = 2*pi*gammas*self.usts.dt

		self.metric, self.L, self.jac, self.amps, self.suf = \
			bhacalc.multipletAnalysis(self.usts.n, gammas, self.prjxns)
		self.jac = 1. / self.jac
		# The variance estimate and its error, for cases when sigma is unknown.
		self.varest = (self.dsqr - self.suf)/(self.dof-2.)
		self.varerr = self.varest * sqrt(2./(self.dof-4.))
		#return (self.suf, self.Jac)
		self.hamps = zeros((self.nf))
		for i in range(self.nf):
			try:
				self.hamps[i] = sqrt(self.amps[2*i]**2 + self.amps[2*i+1]**2)
			except:
				print flist, self.amps[2*i], self.amps[2*i+1]

#...............................................................................
	def amplitudes(self):
		"""Return the amplitudes of the harmonics."""
		
		if self.f_in == None: raise RuntimeError, 'Need to call analyze first!'
		return self.hamps

#...............................................................................
	def noiseEst(self):
		"""Return the estimated noise rms level and its error."""
		
		if self.f_in == None: raise RuntimeError, 'Need to call analyze first!'
		return sqrt(self.varest), 0.5*sqrt(self.varest)*sqrt(2./(self.dof-4.))

#...............................................................................
	def covar(self):
		"""Return the covariance matrix for the best-fit cos & sin amplitudes."""

#...	If we haven't already done the metric calculations, complain!
		if self.f_in == None:
			raise RuntimeError, 'Must call analyze first!'

		if self.var:
			return bhacalc.covar(self.L) * self.var
		else:
			return bhacalc.covar(self.L) * self.varest

#...............................................................................
	def sigmas(self):
		"""Return the standard deviations for the best-fit harmonic amplitudes."""

#...	If we haven't already done the metric calculations, complain!
		if self.f_in == None:
			raise RuntimeError, 'Must call analyze first!'

		if self.var:
			covars = bhacalc.covar(self.L) * self.var
		else:
			covars = bhacalc.covar(self.L) * self.varest
		sigs = zeros((self.nf))
		for i in range(self.nf):
			i1, i2 = 2*i, 2*i+1
			sigs[i] = (self.amps[i1]/self.hamps[i])**2 * covars[i1][i1] + \
					  (self.amps[i2]/self.hamps[i])**2 * covars[i2][i2]
		return sqrt(sigs)

################################################################################
def test_harmonics(m):
	from USTimeSeries import USTimeSeries
	from numpy.random import normal

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
		A3*cos(6*pi*f*times - phi3) + normal(0.,sig,(N_s))
	
	usts = USTimeSeries(data, dt)
	usts.transform(8)
	print 'f range: ', usts.f_range()

	harm = Harmonics(m, usts, sig)
	harm.analyze(f)
	print 'prjxns: ', harm.prjxns
	#print "metric: \n", harm.metric
	#print "Lower triangle:\n", harm.L
	print "Jac: ", harm.jac
	print 'Suff. stat: ', harm.suf
	amps = harm.amplitudes()
	print 'Amps:  ', amps
	print 'Sigs:  ', harm.sigmas()
	print 'Ratio: ', amps/harm.sigmas()

	# Now try it with multiplet.
	print '... Multiplet version:'
	mult = Multiplet(usts, sig)
	flist = []
	for i in range(m):
		flist.append((i+1)*f)
	print 'flist: ', flist
	mult.analyze(flist)
	print "Jac: ", mult.jac
	print 'Suff. stat: ', mult.suf
	amps = mult.amplitudes()
	print 'Amps:  ', amps
	print 'Sigs:  ', mult.sigmas()
	print 'Ratio: ', amps/mult.sigmas()

	# Now estimate the amplitudes using a nonpadded FFT.
	usts.transform()
	fund, rl, im, p_fund = usts.f_info(f)
	A_fund = sqrt( amps[0]**2 + amps[1]**2 )
	for i in range(1,m+1):
		f, rl, im, p = usts.f_info(i*fund)
		print i, i*harm.f, amps[i-1]/amps[0], f, sqrt(p/p_fund)
	print "Done!\n==============="

#test_harmonics(3)

def test_multiplet():
	from USTimeSeries import USTimeSeries
	from numpy.random import normal

	# Specify the parameters for the underlying signal & sampling.
	dt = 1/48000.
	f1 = 984.3750	# The value near 1k for 1024 samples
	f1 = 1007.8125	# The value near 1k for 2048 samples (1/2 way between for 1024)
	#fs = 1031.25
	#fs = 1017.	# A value between grid pts
	#fs = 1000.
	A1 = 100.
	phi1 = 0.
	f2 = f1 - 30.
	A2 = 0.01*A1
	phi2 = pi/2.
	sig = 1.
	N_s = 1024
	nf = 2
	
	# Simulate some real-valued data from a sinusoid with white noise.
	times = arange(N_s)*dt
	data = A1*cos(2*pi*f1*times - phi1) + A2*cos(2*pi*f2*times - phi2) \
		+ 0*normal(0.,sig,(N_s))
	
	usts = USTimeSeries(data, dt)
	usts.transform()
	print 'unpadded f range: ', usts.f_range()
	usts.transform(16)
	print 'padded f range: ', usts.f_range()

	# Try it with multiplet.
	print '... Multiplet version:'
	mult = Multiplet(usts, sig, exact=1)
	flist = (f1, f2)
	print 'flist: ', flist
	mult.analyze(flist)
	print 'used:  ', mult.flist
	print "Jac: ", mult.jac
	print 'Suff. stat: ', mult.suf
	amps = mult.amplitudes()
	print 'Amps:  ', amps
	print 'Sigs:  ', mult.sigmas()
	print 'Ratio: ', amps/mult.sigmas()

	# Now estimate the amplitudes using a nonpadded FFT.
	usts.transform()
	fund, rl, im, p_fund = usts.f_info(flist[0])
	for i in range(nf):
		f, rl, im, p = usts.f_info(flist[i])
		print i, flist[i], mult.flist[i], amps[i]/amps[0], f, sqrt(p/p_fund)
	print "Done!\n==============="

#test_multiplet()


