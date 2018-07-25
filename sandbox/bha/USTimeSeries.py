from scipy import *

__all__ = ['FreqList', 'USTimeSeries']

class FreqList:

	"""A list of frequencies associated with a real FFT."""

	def __init__(self, n_fft, dt):
		self.n_fft = n_fft
		self.dt = dt
		self.ndt = float(n_fft*dt)

	def __len__(self):
		return self.n_fft

	def __getitem__(self, i):
		return i/self.ndt

class USTimeSeries:
	"""
	Uniformly Sampled Time Series.
	"""

	def __init__(self, samples, dt, t0=None):
		"""
		Store the sampled values and the sampling info.
		"""
		if len(shape(samples)) != 1:
			raise NotImplementedError, 'Must be 1-d time series!'
		self.samples = samples
		self.dt = dt
		self.t0 = t0
		if t0: raise ValueError, 'Nonzero t0 not yet implemented!'
		self.n = len(samples)
		self.did_fft = False

	def transform(self, oversamp=1):
		"""
		Perform an FFT of the samples, oversampling by the factor oversamp
		(by zero-padding the data before the FFT).
		"""
		if oversamp < 1: raise ValueError, "Oversample factor must be > 1!"
		self.n_fft = int(self.n*oversamp)
		self.n_freq = self.n_fft/2 + 1	# This count includes f=0
		padded = zeros((self.n_fft))
		padded[0:self.n] = self.samples
		self.fft = fftpack.real_fft(padded)
		self.real = self.fft.real
		self.imag = self.fft.imag
		self.f_vals = FreqList(self.n_fft, self.dt)
		self.did_fft = True

	def f_range(self):
		"""
		Return info about the frequency range spanned by the FFT:
			(n_freq, df, f_max).
		"""
		if not self.did_fft: raise RuntimeError, 'Must call transform first!'
		return self.n_freq, 1./(self.n_fft*self.dt), 0.5/self.dt

	def f_index(self, f):
		"""
		Return the index for the frequency nearest to f in the transform.
		"""
		if not self.did_fft: raise RuntimeError, 'Must call transform first!'
		return int(round(self.n_fft*f*self.dt))

	def f_info(self, f):
		"""
		Return info about the frequency nearest to f in the transform:
			(f_fft, real, imag, PSD).
		"""
		if not self.did_fft: raise RuntimeError, 'Must call transform first!'
		i = int(round(self.n_fft*f*self.dt))
		f_fft = i/(self.dt*self.n_fft)
		return f_fft, self.real[i], self.imag[i], \
			2.*(self.real[i]**2 + self.imag[i]**2)

	def bracket3(self, f):
		"""
		Return 3 frequencies that bracket f and their transform values.
		"""
		i0 = int(round(self.n_fft*f*self.dt))
		if i0 == 0:
			i0 = 0
		elif i0 == self.n_freq-1:
			i0 = self.n_freq-3
		else:
			i0 = i0 - 1
		fvals = arange(i0, i0+3)/(self.dt*self.n_fft)
		rvals = self.real[i0:i0+3]
		ivals = self.imag[i0:i0+3]
		return fvals, rvals, ivals

	def realDTFT(self, f):
		"""
		The real part of the DTFT at f, calculated directly (slowly!).
		"""
		vec = cos(2*pi*f*arange(self.n)*self.dt)
		return sum(vec*self.samples)
		
	def imagDTFT(self, f):
		"""
		The imaginary part of the DTFT at f, calculated directly (slowly!).
		This uses the *negative* exponent FT convention of FFTPACK/FFTW, so
		this is the *negative* of the projection of sine on the samples.
		"""
		vec = sin(2*pi*f*arange(self.n)*self.dt)
		return -sum(vec*self.samples)

	def write_fft(self, file):
		"""
		Write the FFT to file as (f, real, imag, dB(PSD)).
		"""
		if not self.did_fft: raise RuntimeError, 'Must call transform first!'
		for i in range(self.n_freq):
			f_fft = i/(self.dt*self.n_fft)
			rl, im = self.real[i], self.imag[i]
			p = 10. * log10( 2.*(self.real[i]**2 + self.imag[i]**2) )
			s = "%4.4f  %3.3e  %3.3e  %3.3e\n" % (f_fft, rl, im, p)
			file.write(s)

	def __len__(self):
		return len(self.samples)

def test():

	# Specify the parameters for the underlying signal & sampling.
	dt = 1/48000.
	f = 984.3750	# The value near 1k for 1024 samples
	f = 1007.8125	# The value near 1k for 2048 samples (1/2 way between for 1024)
	f = 1031.25
	#f = 1017.	# A value between grid pts
	f = 1000.
	phi = 0.
	A = 1.
	N_s = 1024*16

	# Simulate some real-valued data from a sinusoid with white noise.
	times = arange(N_s)*dt
	data = A*cos(2*pi*f*times - phi)

	usts = USTimeSeries(data, dt)
	usts.transform(8)
	print 'f range: ', usts.f_range()
	f, c, s, p = usts.f_info(1.6e3)
	print 'Near 1.k: ', f
	print 'Fast:  ', c, s, p
	rd, id = usts.realDTFT(f), usts.imagDTFT(f)
	pd = 2*(c**2 + s**2)
	print 'Slow:  ', rd, id, pd
	print 'First f: ', usts.f_vals[0], usts.fft[0]
	print usts.f_info(usts.f_vals[0])
	n = usts.f_range()[0]
	print 'Last f:  ', usts.f_vals[n-1], usts.fft[n-1]
	print usts.f_info(usts.f_vals[n-1])
	print '==============='

if __name__ == '__main__':
	test()

