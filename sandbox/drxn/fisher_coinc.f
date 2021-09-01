c-----------------------------------------------------------------------
	subroutine sig_to_kappa (sig, kappa)

c---	"SIGMA TO KAPPA for fisher dist'n"
c
c	Sets the kappa parameter defining the Fisher dist'n for
c	the direction uncertainty for a burst.  sig is the 68.3%
c	confidence region angular radius in radians.

c+++  Arguments:
	real*8 sig, kappa

	real*8 zero, twopi
	parameter (zero = 0., twopi = 6.283185307179586d0)

cf2py intent(in) sig
cf2py intent(out) kappa

c+++  Locals:
	real*8 one, minus, C, C1, eps
	integer nrmax
	parameter (one = 1., minus = -1., C = 0.683, C1 = 1.-C)
	parameter (eps = 1.e-5, nrmax = 25)
	integer iter
	real*8 kold, f, csig, arg, ekc, ekk, num, den

c---  Find kappa burst using Newton-Raphson.  Use
c---  the Gaussian (small angle) approximation as a starting guess.
	kold = 2.297707 / sig**2
	csig = cos(sig)
	iter = 0
10	continue
	if (iter .ge. nrmax) pause 'Too many iters for kappa!'
	arg = csig - 1.
	ekc = exp(kold*arg)
	if (kold .lt. 250.) then
	    ekk = exp(-2.*kold)
	else
	    ekk = 0.
	endif
	den = 1. - ekk
	num = 1. - ekc
	f = num / den - C
	kappa = kold -
     *          f * den / (2.*num*ekk/den - arg*ekc)
	iter = iter + 1
	if (abs(kappa-kold)/kappa .gt. eps) then
	    kold = kappa
	    goto 10
	endif

	return
	end

c-----------------------------------------------------------------------
	function doublet_lbf (drxn1, kappa1, lsk1, drxn2, kappa2, lsk2)

c---	"DOUBLET Log Likelihood"
c
c	Calculate the log likelihood factor for a doublet.

c+++  Arguments:
	real*8 doublet_lbf
	real*8 drxn1(3), kappa1, lsk1
	real*8 drxn2(3), kappa2, lsk2

cf2py intent(in) drxn1, kappa1, lsk1, drxn2, kappa2, lsk2

c+++  Locals (const = 2*log(4*pi)):
	real*8 const
	parameter (const = 5.0620484939385815d0)
	real*8 tot(3), R

c+++  Functions:
	real*8 lsinhc

	tot(1) = kappa1*drxn1(1) + kappa2*drxn2(1)
	tot(2) = kappa1*drxn1(2) + kappa2*drxn2(2)
	tot(3) = kappa1*drxn1(3) + kappa2*drxn2(3)
	R = sqrt(tot(1)**2 + tot(2)**2 + tot(3)**2)

	doublet_lbf = lsinhc(R) - lsk1 - lsk2

	return
	end

c-----------------------------------------------------------------------
	function triplet_lbf (drxn1, kappa1, lsk1, drxn2, kappa2, lsk2,
     *                      drxn3, kappa3, lsk3)

c---	"TRIPLET Log Likelihood"
c
c	Calculate the log likelihood factor for a triplet.

c+++  Arguments:
	real*8 triplet_lbf
	real*8 drxn1(3), kappa1, lsk1, drxn2(3), kappa2, lsk2
	real*8 drxn3(3), kappa3, lsk3

cf2py intent(in) drxn1, kappa1, lsk1, drxn2, kappa2, lsk2
cf2py intent(in) drxn3, kappa3, lsk3

c+++  Locals (const = 3*log(4*pi)):
	real*8 const
	parameter (const = 7.5930727409078722d0)
	real*8 tot(3), R
	integer n, nd
	parameter (nd = 3)
	real*8 drxn(3,nd), kappa(nd), lsk(nd)

c+++  Functions:
	real*8 lsinhc

c---  Writing it all out is tedious beyond doublets, so put the
c---  direction info in arrays.  For speed, do it here rather than
c---  in Python.
	do 10 n=1, 3
	   drxn(n,1) = drxn1(n)
	   drxn(n,2) = drxn2(n)
	   drxn(n,3) = drxn3(n)
 10	continue
	kappa(1) = kappa1
	kappa(2) = kappa2
	kappa(3) = kappa3
	lsk(1) = lsk1
	lsk(2) = lsk2
	lsk(3) = lsk3

c---  Here's the actual calculation.
	tot(1) = 0.
	tot(2) = 0.
	tot(3) = 0.
	do 20 n=1, 3
	    tot(1) = tot(1) + kappa(n)*drxn(1,n)
	    tot(2) = tot(2) + kappa(n)*drxn(2,n)
	    tot(3) = tot(3) + kappa(n)*drxn(3,n)
20	continue
	R = sqrt(tot(1)**2 + tot(2)**2 + tot(3)**2)

	triplet_lbf = lsinhc(R)
	do 40 n=1, 3
	    triplet_lbf = triplet_lbf - lsk(n)
40	continue

	return
	end

c-----------------------------------------------------------------------
	function multiplet_lbf (nd, drxn, kappa, lsk)

c---	"MULTIPLET Log Likelihood"
c
c	Calculate the log likelihood factor for a quadruplet.

c+++  Parameters:
	real*8 multiplet_lbf
	integer nd
	real*8 drxn(3,nd), kappa(nd), lsk(nd)

cf2py intent(in) nd, drxn, kappa, lsk

c+++  Locals:
	real*8 log4pi
	parameter (log4pi = 2.5310242469692907d0)
	real*8 tot(3), R
	integer n

c+++  Functions:
	real*8 lsinhc

	tot(1) = 0.
	tot(2) = 0.
	tot(3) = 0.
	do 20 n=1, nd
	    tot(1) = tot(1) + kappa(n)*drxn(1,n)
	    tot(2) = tot(2) + kappa(n)*drxn(2,n)
	    tot(3) = tot(3) + kappa(n)*drxn(3,n)
20	continue
	R = sqrt(tot(1)**2 + tot(2)**2 + tot(3)**2)

	multiplet_lbf = lsinhc(R)
	do 40 n=1, nd
	    multiplet_lbf = multiplet_lbf - lsk(n)
40	continue

	return
	end

c-----------------------------------------------------------------------
	function lsinhc (x)

c	Compute log(sinh(x)/x).

c+++  Arguments:
	real*8 lsinhc, x

c+++  Locals:
	real*8 log2
	parameter (log2 = 0.6931471805599453)

	lsinhc = x - log2 - log(x)
	if (x .lt. 20.) lsinhc = lsinhc - exp(-2.*x)

	return
	end

