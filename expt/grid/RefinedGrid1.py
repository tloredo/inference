# Created:  21 Apr 2000 by Tom Loredo
# Last mod: 23 Apr 2000

from math import sqrt, log, floor, ceil, pi
from string import lower
from numpy import fromfunction, array, zeros, sum, exp, clip, shape
from types import TupleType
import ioutils

Lin = "lin"
Log = "log"
rt2pi = sqrt(2*pi)


#===============================================================================
def IsPeak(y):
	"""Return true if the triplet of y values represents a local maximum."""

	return (y[0] < y[1] > y[2]) or (y[0] < y[1] == y[2])

#===============================================================================
def PeakParams(x, y):
	"""Given log f(x) at three points, estimate the parameters of a 
	Gaussian function through those points, (A, mean, sig).
	
	This generalizes the regular grid algorithm of H. J. Sanchez,
	Comp. in Phys., Jul/Aug 1991, 407."""

	w = x[2]-x[0]
	alpha = (x[1]-x[0]) / w
	alpha2 = alpha**2
	fac = 1. - 2.*alpha2
	a12 = y[0] - y[1]
	C = alpha2/fac + a12/(a12-y[1]+y[2])
	d = C*w*fac**2 / (2.*alpha*(1.-alpha) + C*(2.-4.*alpha)*fac)
	mean = x[0] + d
	alpha2 = w*alpha
	sig = alpha2*(0.5*alpha2 - d) / a12
	logA = y[1] + 0.5*(x[1]-mean)**2/sig
	sig = sqrt(sig)
	return (logA, mean, sig)

#x = [0, 1., 3.]
#y = [-8., -2., -2.]
#print PeakParams(x, y)

################################################################################
class LabeledMaxList:
	"""Maintain a sorted list of a fixed number of maximum values of
	a quantity, along with labels of the maxima."""
 

#-------------------------------------------------------------------------------
	def __init__(self, nmax=None, dynrange=None):
		"""Sets the length of the list of labeled maxima."""
		self.nmax = nmax
		if dynrange is None:
			self.dynrange = None
		else:
			self.dynrange = log(dynrange)
		self.list = []

#-------------------------------------------------------------------------------
	def __len__(self):
		return len(self.list)

#-------------------------------------------------------------------------------
	def __getitem__(self, n):
		return self.list[n]

#-------------------------------------------------------------------------------
	def __setitem__(self, n, val):
		self.list[n] = val

#-------------------------------------------------------------------------------
	def __delitem__(self, n, val):
		del self.list[n]

#-------------------------------------------------------------------------------
	def __repr__(self):
		return repr(self.list)

#-------------------------------------------------------------------------------
	def __str__(self):
		return str(self.list)

#-------------------------------------------------------------------------------
	def update(self, q, label):
		"""Update the list of maxima by considering a new value of the
		sorted quantity, q, whose label tuple is given by label."""

#---	Don't waste time if q is lower than the last element of a full list.
		hi = len(self.list)
		if (hi == self.nmax):
			if (q <= self.list[-1][0]):
				return

#---	We'll do this by bisection, as in the bisect.py library module. 
		lo = 0
		while lo < hi:
			mid = (lo+hi)/2
			if (q > self.list[mid][0]):
				hi = mid
			else:
				lo = mid+1
		self.list.insert(lo, (q, label))
		if self.nmax != None:
			if (len(self.list) > self.nmax):
				del self.list[self.nmax]

#---	Check the dynamic range spanned by list members if necessary.
		if self.dynrange != None:
			top = self.list[0][0]
			i = lo		# We must start at lo in case we inserted at the end.
			while i < len(self.list):
				if top - self.list[i][0] > self.dynrange:
					del self.list[i]
				else:
					i = i+1

#-------------------------------------------------------------------------------
	def limitRange(self, dynrange):
		"""Ensure that elements of the current list span a limited dynamic
		range.  This presumes the elements are the *logarithms* of the
		quantity of interest."""

		logrange = log(dynrange)
		top = self.list[0][0]
		i = 1
		while i < len(self.list):
			if top - self.list[i][0] > logrange:
				del self.list[i]
			else:
				i = i+1

if 0:
	print "Testing LML:"
	lml = LabeledMaxList(3)
	lml.update(3.,"three")
	print lml
	lml.update(3.1,"three.one")
	print lml
	lml.update(1.,"one")
	print lml
	lml.update(10.,"ten")
	print lml
	lml.update(20.,"twenty")
	print lml
	lml.update(4.,"four")
	print lml
	from copy import copy
	lml2 = copy(lml.list)
	lml.limitRange(3.)
	print lml
	print lml2

################################################################################
class BasicInterval:
	"""A basic, positively-directed (closed) interval of equispaced samples."""

#-------------------------------------------------------------------------------
	def __init__(self, lo, hi, n):
		"""Define the interval and sample density.
		This assumes the direction has been verified in the calling code."""

		self.lo = lo
		self.hi = hi
		if (n <= 1): raise ValueError, "Need n > 1!"
		self.n = n
		self.step = (hi-lo) / (n-1.)

#-------------------------------------------------------------------------------
	def __str__(self):
		return str( (self.lo, self.hi, self.n, self.step) )

#-------------------------------------------------------------------------------
	def locate(self, x):
		"""Locate x in the interval, returning the subinterval
		number and endpoints (i, l, h).  This assumes x is
		in fact in the interval.  If x=hi, then i=n, and l=h=hi."""

		i = int( floor((x-self.lo)/self.step) )
		l = self.lo+i*self.step
		if (i == self.n-1):
			h = l
		else:
			h = l + self.step
		return (i, l, h)

#-------------------------------------------------------------------------------
	def resample(self, n):
		"""Change the # of samples in the interval."""

		self.n = n
		self.step = (self.hi-self.lo) / (n-1.)

if 0:
	print "Testing BasicInterval:"
	bi = BasicInterval(5., 10., 6)
	print bi.locate(5.)
	print bi

################################################################################
class RefinedIntervalList1:
	"""Maintain a list of continguous 1-d sampled intervals, bounded by the
	initial interval.  Ensure that all subintervals are sampled more
	densely than the base interval."""

#-------------------------------------------------------------------------------
	def __init__(self, lo, hi, n, type="linear"):
		"""Define the overall interval and base sample density."""

#---  Save some defining info.  Note that the internal representation
#---  is always a positively directed interval, and if the list is of
#---  log type, the log boundaries are stored.
		if (n <= 1):
			raise ValueError, "Number of samples must be >1!"
		if (lo == hi):
			raise ValueError, "Null base interval!"
		elif (lo < hi):
			self.posQ = 1
			self.lo = lo
			self.hi = hi
		else:
			self.posQ = 0
			self.lo = hi
			self.hi = lo
		type = lower(type)
		if (type == "lin" or type == "linear"):
			self.linQ = 1
		elif (type == "log" or type == "logarithmic"):
			self.linQ = 0
			if (lo == 0. or hi == 0.):
				raise ValueError, "Zero boundary for log interval!"
			self.lo = log(lo)
			self.hi = log(hi)
			raise ValueError, "***Log intervals are untested!***"
		else:
			raise ValueError, "Interval type must be Lin or Log!"
#		self.base = ( self.lo, self.hi, n, (self.hi-self.lo) / (n-1.) )
		self.base = BasicInterval(self.lo, self.hi, n)

#---  Initialize the list with the base interval.
		self.ilist = [ self.base ]
		self.nint = 1
		self.npts = n		# This is correct only at the end of updates/merges.
		self.cumpts = [ n ]	# Ditto for this.

#-------------------------------------------------------------------------------
	def __str__(self):

		s = "[\n"
		for i in self.ilist:
			s = s + str(i) + "\n"
		s = s + "]\n"
		return s

#-------------------------------------------------------------------------------
	def __getitem__(self, n):

		return self.ilist[n]

#-------------------------------------------------------------------------------
	def _del(self, i):
		"""Insert an interval into the interval list.
		This is private because it corrupts the adjacency and point count."""

		del self.ilist[i]
		self.nint = self.nint - 1

#-------------------------------------------------------------------------------
	def _insert(self, i, intvl):
		"""Insert an interval into the interval list.
		This is private because it can corrupt the adjacency and point count."""

		self.ilist.insert(i, intvl)
		self.nint = self.nint + 1

#-------------------------------------------------------------------------------
	def update(self, lo, hi, n):
		"""Update the interval list to include a new interval."""

#---  First, make sure we are in the base interval.  Truncate if bigger;
#---  complain if it doesn't overlap at all.
		if (not self.linQ):
			lo = log(lo)
			hi = log(hi)
		if (self.posQ):
			lo = max(lo, self.lo)
			hi = min(hi, self.hi)
		else:
			lo = max(hi, self.lo)
			hi = min(lo, self.hi)
		if (hi <= lo):
			raise ValueError, "Interval does not overlap base interval!"

#---  Identify which intervals include the new boundaries.
		for i in range(self.nint):
			if (self.ilist[i].lo <= lo and self.ilist[i].hi > lo):
				nlo = i
				int_lo = self.ilist[i]
				break
		for i in range(self.nint):
			if (self.ilist[i].lo < hi and self.ilist[i].hi >= hi):
				nhi = i
				int_hi = self.ilist[i]
				break

#---  Now expand the new interval to abutt grid pts in existing intervals,
#---  updating the # of points to maintain the density.
		step = (hi-lo) / (n-1.)
		lo = int_lo.locate(lo)[1]
		hi_int = int_hi.locate(hi)
		if (hi == hi_int[1]):
			hi = hi_int[1]
		else:
			hi = hi_int[2]
		n = int( ceil( (hi-lo)/step ) ) + 1
		step = (hi-lo) / (n-1.)

#---  If the new intvl is entirely in an existing one, we have to insert it,
#---  unless it undersamples the interval.  Be careful of it exactly matching
#---  a boundary, in which case fewer inserts are needed.
		if (nlo == nhi):
			if (step >= int_lo.step): return
			self._del(nlo)
			if (lo > int_lo.lo):
				m = int( ceil((lo-int_lo.lo)/int_lo.step) ) + 1
				self._insert(nlo, BasicInterval(int_lo.lo, lo, m))
				nlo = nlo + 1
			self._insert(nlo, BasicInterval(lo, hi, n))
			if (int_lo.hi > hi):
				m = int( ceil((int_lo.hi-hi)/int_lo.step) ) + 1
				self._insert(nlo+1, BasicInterval(hi, int_lo.hi, m))
			self.merge()
			return

#---  Otherwise, we have to break up the intervals containing the boundaries,
#---  and possibly resample intervals completely covered by the new one.
		if (step < int_lo.step):
			self._del(nlo)
			if (lo > int_lo.lo):
				m = int( ceil((lo-int_lo.lo)/int_lo.step) ) + 1
				self._insert(nlo, BasicInterval(int_lo.lo, lo, m))
				nlo = nlo + 1
				nhi = nhi + 1		# Shift nhi since we're adding an intvl.
			m = int( ceil((int_lo.hi-lo)/step) ) + 1
			self._insert(nlo, BasicInterval(lo, int_lo.hi, m))
			nlo = nlo + 1
		if (step < int_hi.step):
			self._del(nhi)
			m = int( ceil((hi-int_hi.lo)/step) ) + 1
			self._insert(nhi, BasicInterval(int_hi.lo, hi, m))
			if (int_hi.hi > hi):
				m = int( ceil((int_hi.hi-hi)/int_hi.step) ) + 1
				self._insert(nhi+1, BasicInterval(hi, int_hi.hi, m))
				# We needn't update nhi here since range won't include it.
		for i in range(nlo, nhi):
			if (self.ilist[i].step > step):
				n = int( ceil((self.ilist[i].hi-self.ilist[i].lo)/
					      step) ) + 1
				self.ilist[i].resample(n)
		self.merge()

#-------------------------------------------------------------------------------
	def merge(self):
		"""Merge adjacent intervals with a common step size to make
		the interval list as small as possible.
		This also makes a correct point count (not counting overlapping
		points at boundaries)."""

		i = 0
		self.npts = self.ilist[0].n - 1
		self.cumpts = [ self.npts ]
		while 1:
			if (self.nint < i+2):		# Return if there is no next intvl.
				self.npts = self.npts + 1
				self.cumpts[-1] = self.npts
				return
			self.npts = self.npts + self.ilist[i+1].n - 1
			self.cumpts.append(self.npts)
			if (self.ilist[i].step == self.ilist[i+1].step):	# Merge.
				lo = self.ilist[i].lo
				hi = self.ilist[i+1].hi
				n = self.ilist[i].n + self.ilist[i+1].n - 1
				self._del(i)
				self._del(i)
				self._insert(i, BasicInterval(lo, hi, n))
				self.cumpts[i] = self.cumpts[i+1]
				del self.cumpts[-1]
			else:
				i = i + 1

#-------------------------------------------------------------------------------
	def nonbase(self):
		"""Return a list of all non-base (i.e., refined) intervals."""

		l = []
		for intvl in self.ilist:
			if (intvl.step != self.base.step):
				l.append(intvl)
		return l

#-------------------------------------------------------------------------------
	def locate(self, x):
		"""Return the interval containing x."""

#---  Make sure x is in range.
		if (not self.linQ): x = log(x)
		if (x < self.lo or x > self.hi):
			raise ValueError, "x out of range!"

#---  We'll treat the intervals as [), except the last which is [].
		if (x == self.hi):
			return self.ilist[-1]

#---  Use bisection.
		lo = 0
		hi = self.nint
		mid = (lo+hi)/2
		while 1:
			if (self.ilist[mid].lo <= x < self.ilist[mid].hi):
				break
			elif (x < self.ilist[mid].lo):
				hi = mid
			else:
				lo = mid
			mid = (lo+hi)/2
		return self.ilist[mid]

#-------------------------------------------------------------------------------
	def nlocate(self, n):
		"""Return the interval containing pt number n (0-based);
		include the intvl number and the location of n in the intvl."""

#---  Make sure n is in range.
		if (n < 0 or n > self.npts-1):
			raise ValueError, "n out of range!"

#---  We'll treat the intervals as [), except the last which is [].
		if (n == self.npts-1):
			return (self.nint-1, self.ilist[-1].n-1, self.ilist[-1])

#---  Use bisection.
		lo = 0
		hi = self.nint
		while lo < hi:
			mid = (lo+hi)/2
			if (n+1 > self.cumpts[mid]):
				lo = mid + 1
			else:
				hi = mid
		if (lo == 0):
			indx = n
		else:
			indx = n - self.cumpts[lo-1]
		return (lo, indx, self.ilist[lo])

#-------------------------------------------------------------------------------
if 0:
	print "Testing RefinedIntervalList1:"
	cil = RefinedIntervalList1(5., 10., 6)
	cil.merge()
	cil.update(7., 8., 5)
	print cil
	print "# of pts: ", cil.npts
	print "cum: ", cil.cumpts
	cil.update(6., 9., 16)
	cil.update(7., 8., 11)
	print cil
	print "# of pts: ", cil.npts
	print "cum: ", cil.cumpts
#	cil.merge()
#	print cil
	print map(str,cil.nonbase())
	print cil.locate(7.5)
	print cil.nlocate(6)

################################################################################
class BasicGrid1:
	"""Create a 1-d linear or logarithmic grid with automated evaluation of a 
	function and optional bookkeeping of local maxima."""


#-------------------------------------------------------------------------------
	def __init__(self, lo, hi, n, type="linear", lml=None):
		"""Define the grid abscissa."""

#---  Copy basic defining data and check for possible illegal values.
		self.lo = lo
		self.hi = hi
		if (n <= 1):
			raise ValueError, "Number of points must be >1!"
		self.n = n
		type = lower(type)
		if (type == "lin" or type == "linear"):
			self.linQ = 1
		elif (type == "log" or type == "logarithmic"):
			self.linQ = 0
			if (lo == 0. or hi == 0.):
				raise ValueError, "Zero boundary for log grid!"
		else:
			raise ValueError, "Grid type must be Lin or Log!"
		self.lml = lml
		self.ntupvals = None

#---  Find the appropriate step size.
		if (self.linQ):
			self.step = (hi-lo)/(n-1.)
		else:
			self.step = log(hi/lo)/(n-1.)

#---  Initialize abscissas and ordinates to None.
		self.absc = None
		self.vals = None

#-------------------------------------------------------------------------------
	def __str__(self):

		s = '[\n'
		for i in xrange(self.n):
			s = s + str(self.point(i)) + '\n'
		s = s + ']\n'
		return s

#-------------------------------------------------------------------------------
	def stepfunc(self, func, select=None, notify=None):
		"""Gather function values by applying func step by step through the grid.
		Note that we don't store the abscissas here to save space.
		If the function returns tuples, select picks the element to monitor
		in the list of maxima.
		The notify parameter allows monitoring of progress."""

		x1, x2 = None, None
		for i in xrange(self.n):
			if (self.linQ):
				x = self.lo + i*self.step
			else:
				x = self.lo * exp(i*self.step)

#---  First time thru, check whether func returns tuples and create storage
#---  appropriately.
			if i == 0:
				f = func(x)
				if ( type(f) is TupleType ):
					if select is None: raise ValueError, "must set select!"
					self.ntupvals = len(f)
					self.vals = zeros((self.n, self.ntupvals), float)
				else:
					self.vals = zeros((self.n), float)
				self.vals[i] = f
			else:
				self.vals[i] = func(x)

			if notify != None:
				if i%notify == 0: print i, x, self.vals[i]

#---  Here is the check for a local maximum, with special treatment for
#---  possible "half peaks" at the ends.  These have fudged params
#---  from a linear fit to 2 pts to find dx causing the log to drop 1/2.
#*** check the endpt trtmt.
			if (self.lml != None and i > 1):
				if self.ntupvals:
					y3 = self.vals[i-2:i+1, select]
				else:
					y3 = self.vals[i-2:i+1]
				if (IsPeak(y3)):
					logA, mean, sig = PeakParams((x1, x2, x), y3)
					area = logA + log(sig*rt2pi)
					self.lml.update(area, (i-1, logA, mean, sig))
				if (i == 2 and (y3[0] > y3[1])):
					sig = -0.5*(x2-x1)/(y3[1]-y3[0])
					logA, mean = y3[0], x1
					area = logA + log(0.5*sig*rt2pi)
					self.lml.update(area, (i-2, logA, mean, sig))
				if (i == (self.n-1) and (y3[1] < y3[2])):
					sig = 0.5*(x-x2)/(y3[2]-y3[1])
					logA, mean = y3[2], x
					area = logA + log(0.5*sig*rt2pi)
					self.lml.update(area, (i, logA, mean, sig))
			x1, x2 = x2, x
		
#-------------------------------------------------------------------------------
#*** test treatment of endpoints, which was just copied from above
	def vecfunc(self, func, select=None):
		"""Gather function values by applying func to a vector of grid abscissas."""

#---  First evaluate func on a vector of abscissa values.
		if (self.linQ):
			self.absc = self.lo + self.step*fromfunction(lambda i: i, [self.n])
		else:
			self.absc = self.lo * exp(self.step*fromfunction(lambda i: i, [self.n]))
		self.vals = func(self.absc)

#---  Note whether or not the values are "tuples".
		if len(shape(self.vals)) == 2:
			self.ntupvals = shape(self.vals)[1]
		else:
			self.ntupvals = None

#---  Now look for local maxima, with special treatment for
#---  possible "half peaks" at the ends.
		if (self.lml != None):
			for i in xrange(2,self.n):
				if self.ntupvals:
					y3 = self.vals[i-2:i+1, select]
				else:
					y3 = self.vals[i-2:i+1]
				x3 = self.absc[i-2:i+1]
				if (IsPeak(y3)):
					logA, mean, sig = PeakParams(x3, y3)
					area = logA + log(sig*rt2pi)
					self.lml.update(area, (i-1, logA, mean, sig))
				if (i == 2 and (y3[0] > y3[1])):
					#A, mean, sig = PeakParams(self.absc[i-2:i+1], y3)
					#area = 0.5*A*sig*rt2pi
					#self.lml.update(area, (i-2, A, mean, sig))
					sig = -0.5*(x3[1]-x3[0])/(y3[1]-y3[0])
					logA, mean = y3[0], x3[0]
					area = logA + log(0.5*sig*rt2pi)
					self.lml.update(area, (i-2, logA, mean, sig))
				if (i == (self.n-1) and (y3[1] < y3[2])):
					#A, mean, sig = PeakParams(self.absc[i-2:i+1], y3)
					#area = 0.5*A*sig*rt2pi
					#self.lml.update(area, (i, A, mean, sig))
					sig = 0.5*(x3[2]-x3[1])/(y3[2]-y3[1])
					logA, mean = y3[2], x3[2]
					area = logA + log(0.5*sig*rt2pi)
					self.lml.update(area, (i, logA, mean, sig))

#-------------------------------------------------------------------------------
	def point(self, i):
		"""Return the abscissa, width, and value associated with point i.
		The width is the size of the step to the next point, and is
		suitable for trapezoid rule quadrature; it vanishes for the
		right endpoint."""

#---  Make sure i is legal.
		if (i < 0 or i >= self.n):
			raise ValueError, "Illegal index!"

#---  Get abscissa value, using a stored value if available.  If the ordinates
#---  have been calculated, get the ordinate; otherwise use None.
		if (self.absc == None):
			if (self.linQ):
				x = self.lo + i*self.step
			else:
				x = self.lo * exp(i*self.step)
		else:
			x = self.absc[i]
		if (self.vals == None):
			v = None
		else:
			v = self.vals[i]

#---  Get the size of the step to the next point (0 for last point!).
		if (i == self.n-1):
			w = 0.
		else:
			if (self.linQ):
				w = self.step
			else:
				w = x * (exp(self.step)-1.)

		return (x, w, v)

#-------------------------------------------------------------------------------
	def trapzd(self, select=None):
		"""Return the trapezoid rule quadrature of the values on the grid."""

#---  An empty grid gives an error!
		if (self.vals == None): raise RuntimeError, "grid values not yet set!"

#---  Pick out the values to use if we have a grid of tuples.
		if self.ntupvals:
			if select is None: raise ValueError, "must set select!"
			vals = self.vals[:, select]
		else:
			vals = self.vals

#---  The linear case can use array math.
		if (self.linQ):
			return self.step*(sum(vals) - 0.5*(vals[0]+vals[-1]))

#---  The log case can either step thru the grid to save space, or
#---  make an array of widths to save time.  We do the latter here.
		w = self.absc
		if (w == None):
			w = self.lo * exp(self.step*fromfunction(lambda i: i, [self.n]))
		w = w[1:] - w[:-1]
		return sum(0.5 * w * (vals[1:] + vals[:-1]))


#-------------------------------------------------------------------------------
	def trapzdLog(self, select=None):
		"""Return the log of the trapezoid rule quadrature of the values on the grid, 
		treating them as the log of the integrand."""

#---  An empty grid gives an error!
		if (self.vals == None): raise RuntimeError, "grid values not yet set!"

#---  Make an array of step sizes.
		if (self.linQ):
			w = zeros(self.n-1, float) + self.step
		else:
			w = self.absc
			if (w == None):
				w = self.lo * exp(self.step*fromfunction(lambda i: i, [self.n]))
			w = w[1:] - w[:-1]

#---  Pick out the values to use if we have a grid of tuples.
		if self.ntupvals:
			if select is None: raise ValueError, "must set select!"
			vals = self.vals[:, select]
		else:
			vals = self.vals

#---  Calculate with respect to the max to avoid overflow.
		fmax = max(vals)
		fvals = exp(clip(vals-fmax,-725.,0.))
		s = sum(0.5 * w * (fvals[1:] + fvals[:-1]))
		return fmax + log(s)

#-------------------------------------------------------------------------------
	def subgrid(self, lo, hi):
		"""Create a subgrid of this grid, with data elements pointing
		to slices of this grid's data as appropriate."""

#---  Only subgrid if (lo,hi) are grid points.
		if (lo < self.lo or hi > self.hi):
			raise ValueError, "subgrid out of range!"
		if (self.linQ):
			lnum = lo-self.lo
			rlo = lnum/self.step
			hnum = hi-self.lo
			rhi = hnum/self.step
		else:
			lnum = log(lo/self.lo)
			rlo = lnum/self.step
			hnum = log(hi/self.lo)
			rhi = hnum/self.step
		nlo = int(round(rlo,2))
		if (abs(nlo*self.step-lnum) > 1.e-3*self.step):
			raise ValueError, "lo doesn't match a grid pt!"
		nhi = int(round(rhi,2))
		if (abs(nhi*self.step-hnum) > 1.e-3*self.step):
			raise ValueError, "hi doesn't match a grid pt!"

#---  Make the grid and copy the data.
		n = nhi - nlo + 1
		if (self.linQ):
			grid = BasicGrid1(lo, hi, n, Lin, self.lml)
		else:
			grid = BasicGrid1(lo, hi, n, Log, self.lml)
		if (self.absc != None):
			grid.absc = self.absc[nlo:nhi+1]
		if (self.vals != None):
			grid.vals = self.vals[nlo:nhi+1]
		grid.ntupvals = self.ntupvals
		return grid
		
#-------------------------------------------------------------------------------
if 0:
	print "Testing BasicGrid1:"
	from numpy import log
	lml = LabeledMaxList(4)
	sg1 = BasicGrid1(10, 1000, 20, lml=lml)
#	def f(x): return -((x-300.)*(x-300.))/2.e4
#	def f(x): return log(x)
	def f(x):
		lx = zeros((len(x),2), float)
		lx[:,0] = log(x)
		lx[:,1] = x
		return lx
	sg1.vecfunc(f,0)
#	print sg1
	print sg1.trapzd(0), sg1.trapzdLog(0)
	print sg1.lml, sg1.point(2)

################################################################################
class RefinedGrid1:
	"""Create a 1-d grid with automated evaluation of a function and single-pass
	refinement in neighborhoods of local maxima."""

#-------------------------------------------------------------------------------
	def __init__(self, lo, hi, n, type="linear", nmax=10, dynrange=None,
			     nsig=4, sigpts=11):
		"""Define the base grid and parameters defining the refinement.
		There will be sigpts new points per each +-sigma range, for
		nsig sigmas."""

		self.nmax = nmax
		self.lml = LabeledMaxList(nmax=nmax, dynrange=dynrange)
		self.type = type
		self.base = BasicGrid1(lo, hi, n, type, lml=self.lml)
		self.ilist = RefinedIntervalList1(lo, hi, n, type)
		self.nsig = nsig
		self.rpts = nsig*(sigpts-1) + 1
		self.totpts = n
		self.func = None
		self.glist = [ self.base ]

#-------------------------------------------------------------------------------
	def __len__(self):
		"""Total # of pts in the refined grid."""

		return self.totpts

#-------------------------------------------------------------------------------
	def firstpass(self, f, select=None, notify=None):
		"""Apply func on the base grid."""

		self.func = f
		self.select = select
		self.base.stepfunc(f, select=select, notify=notify)

#-------------------------------------------------------------------------------
	def refine(self, dynrange=None, notify=None):
		"""Refine the grid at the maxima.
		Note:  The boundary values are recomputed for each subgrid,
		resulting in a small inefficiency."""

#---  Make sure we've done a firstpass.
		if (self.func == None):
			raise RuntimeError, "Must do firstpass before refine!"

#---  Limit the dynamic range to consider for the refinement, if requested.
		if dynrange != None:
			self.lml.limitRange(dynrange)

#---  First, refine the interval list.
		for (A, label) in self.lml:
			lo = label[2] - self.nsig*label[3]
			hi = label[2] + self.nsig*label[3]
			self.ilist.update(lo, hi, self.rpts)
		self.totpts = self.ilist.npts

#---  Make a list of grids associated with each interval.  Intervals with step
#---  size equal to the base step size are assigned subgrids of the base grid.
		self.glist = []
		self.rlml = LabeledMaxList(self.nmax)
		for intvl in self.ilist:
			if (intvl.step == self.ilist.base.step):
				self.glist.append(self.base.subgrid(intvl.lo,intvl.hi))
			else:
				self.glist.append(BasicGrid1(intvl.lo, intvl.hi, intvl.n,
					self.type, lml=self.rlml))

#---  Go thru the grid list, evaluating the function and compiling a refined
#---  list of labeled maxima.  But don't duplicate the work done on the base!
		if notify != None:
			print "Refining among ",len(self.glist), " subgrids:"
		for grid in self.glist:
			if (grid.step != self.base.step):
				if notify != None:
					print ">>>Refining over [",grid.lo,", ",grid.hi,"]..."
				grid.stepfunc(self.func, select=self.select, notify=notify)

#-------------------------------------------------------------------------------
	def point(self, n):
		"""Return (absc, width, value) associated with point # n (0-based)."""

#---  Find the interval containing pt n; return the data.
		(i, indx, intvl) = self.ilist.nlocate(n)
		return self.glist[i].point(indx)

#-------------------------------------------------------------------------------
	def trapzd(self, select=None):
		"""Return the trapezoid rule quadrature based on the current grid."""

#---  If there is no refinement, just use the base grid.
		if (self.glist == None):
			return self.base.trapzd(select=select)

#---  Otherwise, sum over the refined grids.
		sum = 0.
		for grid in self.glist:
			sum = sum + grid.trapzd(select=select)
		return sum

#-------------------------------------------------------------------------------
	def trapzdLog(self, select=None):
		"""Return the log trapezoid rule quadrature based on the current grid,
		treating the function values as logs of the integrand."""

#---  If there is no refinement, just use the base grid.
		if (self.glist == None):
			return self.base.trapzdLog(select=select)

#---  Otherwise, sum over the refined grids taking care to avoid overflows.
		lsum = self.glist[0].trapzdLog(select=select)
		for grid in self.glist[1:]:
			l = grid.trapzdLog(select=select)
			if (l > lsum): l, lsum = lsum, l
			lsum = lsum + log(1. + exp(l-lsum))
		return lsum

#-------------------------------------------------------------------------------

if 0:
	print "Testing RefinedGrid1:"
	rg = RefinedGrid1(0., 100., 100, sigpts=11, nsig=5)
	def bexp(x):
		if (abs(x) < 300.):
			return exp(x)
		else:
			return 1.e-300
	def f(x): return log( bexp(-(x-50.234)**2/2.) + .2*bexp(-(x-72.432)**2/2.) )
#	def f(x): return -(x-50.234)**2/2., bexp(-(x-50.234)**2/2.)
#	def f(x): return log( 2*x )
	rg.firstpass(f,0)
	print rg.base.lml
# Use the next with log(2*x) & print the contributions to trapzdLog
#	rg.lml = rg.base.lml = [(2., (1, 0., 10., 2.)), (1., (2, 0., 70., 2.))]
	print "Base trapzd =    ", rg.trapzdLog(0)
	rg.refine()
	print "Refined trapzd = ", rg.trapzd(1)
	print rg.ilist
	print rg.point(98)
	ofile = ioutils.OpenText('rg.dat')
	print rg.totpts
	for i in range(rg.totpts):
		point = rg.point(i)
		line = "%5i  %#10.4g  %#10.4g\n" % (i, point[0], point[2])
		ofile.write(line)
	ofile.close()
