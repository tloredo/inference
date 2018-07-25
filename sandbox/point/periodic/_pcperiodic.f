c	# OF BINS CHECK!!
c-----------------------------------------------------------------------
c	pcd.f:  "Piecewise Constant w/ Derivatives"
c
c	Subprograms for Bayesian analysis of a piecewise-constant
c	model for Poisson event locations w/ frequency derivatives
c	known.
c
c	Many thanks to Phil Gregory for inspiration and help in 
c	thinking about this problem!
c
c	Adapted from pc.f, 27 Nov 92   TJL
c   Adapted for Python  Jan 2006   TJL
c   Modified:
c       30 may 2006 - Separate phases* to reduce computation
c                     when no derivs are used
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	subroutine phases0(w, data, phases, ndat)

c---	PHASES with no freqeuncy derivatives
c
c---	Calculate the phases of the arrival times in data, given
c	w.  The phases are in [0,1).

c---  Arguments:
	integer ndat
	real*8 w
	real*8 data(ndat), phases(ndat)

c	phases is an f2py "inplace" argument, so we don't waste time 
c	reallocating memory on each call.
c
cf2py intent(in) w, ndat, data
cf2py intent(inplace) phases

c---  Locals:
	integer n
	real*8 twopi
	parameter (twopi = 6.2831853071795862d0)
	real*8 w2

        w2 = w / twopi
	do n=1, ndat
c	    phases(n) = ((wddot*d/6. + 0.5*wdot)*d + w) * d / twopi
	    phases(n) = w2 * data(n)
	enddo

	return
	end

c-----------------------------------------------------------------------
	subroutine phasesd(w, wdot, data, phases, ndat)

c---	PHASES with a single frequency Derivative
c
c---	Calculate the phases of the arrival times in data, given
c	w and wdot.  The phases are in [0,1).
c
c	phases is an f2py "inplace" argument, so we don't waste time 
c	reallocating memory on each call.

c---  Arguments:
	integer ndat
	real*8 w, wdot
	real*8 data(ndat), phases(ndat)

cf2py intent(in) w, wdot, ndat, data
cf2py intent(inplace) phases

c---  Locals:
	integer n
	real*8 twopi
	parameter (twopi = 6.2831853071795862d0)
	real*8 d, wd2, w2

        w2 = w / twopi
        wd2 = 0.5*wdot / twopi
	do n=1, ndat
	    d = data(n)
c	    phases(n) = ((wddot*d/6. + 0.5*wdot)*d + w) * d / twopi
	    phases(n) = (wd2*d + w2) * d
	enddo

	return
	end

c-----------------------------------------------------------------------
	subroutine phasesdd(w, wdot, wddot, data, phases, ndat)

c---	PHASES with frequency Double Derivative
c
c---	Calculate the phases of the arrival times in data, given
c	w, wdot, and wddot.  The phases are in [0,1).
c
c	phases is an f2py "inplace" argument, so we don't waste time 
c	reallocating memory on each call.

c---  Arguments:
	integer ndat
	real*8 w, wdot, wddot
	real*8 data(ndat), phases(ndat)

cf2py intent(in) w, wdot, wddot, ndat, data
cf2py intent(inplace) phases

c---  Locals:
	integer n
	real*8 twopi
	parameter (twopi = 6.2831853071795862d0)
	real*8 d, wdd2, wd2, w2

        w2 = w / twopi
        wd2 = 0.5*wdot / twopi
        wdd2 = wddot/6. / twopi
	do n=1, ndat
	    d = data(n)
c	    phases(n) = ((wddot*d/6. + 0.5*wdot)*d + w) * d / twopi
	    phases(n) = ((wdd2*d + wd2)*d + w2) * d
	enddo

	return
	end

c-----------------------------------------------------------------------
	subroutine binbounds(nbins, phibins)

c---	BIN BOUNDarieS
c
c	Calculate the phases of the left boundaries of the bins.
c
c	phibins is an f2py "inplace" argument, so we don't waste time 
c	reallocating memory on each call.

c---  Arguments:
	integer nbins
	real*8 phibins(nbins)

cf2py intent(in) nbins
cf2py intent(inplace) phibins

c---  Locals:
	integer j
	real*8 one
	parameter (one = 1.d0)
	real*8 phib

	phib = one / nbins
	do j=1, nbins
	    phibins(j) = (j-one) * phib
	enddo

	return
	end

c-----------------------------------------------------------------------
	subroutine fbindata (ndat, phases, phi, nbins, phibins, 
     *                       bins, dphi)

c---	Fast BIN DATA
c
c---	Bin the event phase data, shifted by phi1; don't calculate the 
c	binning factors, to speed things up when they are not needed.
c
c	The phase boundaries for the bins should be pre-calculated and
c	passed in phibins (see binbounds, above).
c
c	phi and dphi are in [0,2*pi).  phases and phibins are in [0,1)
c	to speed the calculation.
c
c	This also returns dphi, the smallest change in phase that will
c	shift an event to another bin.  Bins are semiopen: [).  
c
c	bins is an f2py "inplace" argument, so we don't waste time 
c	reallocating memory on each call.

c---  Arguments:
	integer ndat, nbins, bins(nbins)
	real*8 phases(ndat), phi, phibins(nbins), dphi

cf2py intent(in) ndat, phases, phi, nbins, phibins
cf2py intent(inplace) bins
cf2py intent(out) dphi

c---  Locals:
	integer j, n
	real*8 zero, one, tiny, twopi
	parameter (zero = 0.d0, one = 1.d0, tiny = 1.d-10)
	parameter (twopi = 6.2831853071795862d0)
	real*8 phase, phi1

c---  Zero the bins.
	do 15 j=1, nbins
	    bins(j) = 0
15	continue

c---  Bin the data, keeping track of the smallest forward distance to a
c---  bin boundary.  We use phases over [0,1) rather than [0,twopi)
c---  so we can use simple subtraction, without divisions,
c---  to perform a modulo (use of the mod function is much slower!).
	dphi = one / nbins
	phi1 = phi / twopi
	do 20 n=1, ndat
	    phase = phases(n) + phi1
	    phase = phase - int(phase)
	    j = nbins*phase
	    j = j + 1
	    bins(j) = bins(j) + 1
	    if (j .lt. nbins) then
	        phase = phibins(j+1) - phase
	    else
	        phase = one - phase
	    endif
	    if (phase .lt. dphi .and. phase .gt. zero) dphi = phase
20	continue

c---  Add nbins*tiny to dphi, to make sure that it moves the event
c---  to the next bin, even with roundoff error.
c	if (dphi .lt. 1.e-10) write(*,'(4g13.5)') w,phi,dphi,dphi/nbins
	dphi = twopi * (dphi + nbins*tiny)

	return
	end

c-----------------------------------------------------------------------
	function limult (ndat, nbins, bins)

c---	Log Inverse Multiplicity
c
c	The natural log of the inverse of the multiplicity of the
c	binned distribution.

c---  Arguments:
	real*8 limult
	integer ndat, nbins, bins(nbins)

cf2py intent(in) ndat, nbins, bins

c---  Locals:
	real*8 factln
	integer n

	limult = -factln(ndat)
	do 20 n=1, nbins
	    if (bins(n) .lt. 0) write(*,*) 'lim: ', n, bins(n)
	    limult = limult + factln(bins(n))
20	continue

	return
	end

c-----------------------------------------------------------------------
	function plimult (ndat, nbins, bins)

c---	Part of the Log Inverse Multiplicity
c
c	The natural log of the inverse of the multiplicity of the
c	binned distribution, without the ndat! contribution, to
c	save a little time (one gammln call per call, if ndat>100).

c---  Arguments:
	real*8 plimult
	integer ndat, nbins, bins(nbins)

cf2py intent(in) ndat, nbins, bins

c---  Locals:
	real*8 factln
	integer n

	plimult = 0.
	do 20 n=1, nbins
	    plimult = plimult + factln(bins(n))
20	continue

	return
	end

c-----------------------------------------------------------------------
	function ll_wpb (ndat, nbins, bins)

c---	Log quasiLikelihood for W (omega), Phi, and Bin number
c
c	The log of the joint quasilikelihood for omega, phi, and the # 
c	of bins, with the fractions marginalized away.  This is NOT
c	normalized.
c
c	The w and phi dependence enter thru the binning, done by
c	subroutine fbindata.  By "w" we here mean any parameters
c	determining event phases, e.g., w, wdot, wddot.

c---  Arguments:
	real*8 ll_wpb
	integer ndat, nbins, bins(nbins)

cf2py intent(in) ndat, nbins, bins

c---  Locals:
	real*8 factln, limult
	integer arg

	arg = nbins - 1
	ll_wpb = factln(arg)
	arg = ndat + arg
	ll_wpb = ll_wpb - factln(arg)
	ll_wpb = ll_wpb + (1.*ndat)*nbins +factln(ndat) +
     +           limult(ndat, nbins, bins)

	return
	end

c-----------------------------------------------------------------------
	function lql_wp (ndat, nbins, bins)

c---	Log quasiLikelihood for W (omega) and Phi
c
c	The log of the joint quasilikelihood for omega and phi,
c	cond. on the # of bins, with the fractions marginalized away.  
c	This is NOT normalized.
c
c	The w and phi dependence enter thru the binning, done by
c	subroutine fbindata.
c
c	This is just ll_wpb, without the nbin-dependent terms which
c	may (negligibly...) slow things down.  That is, it is
c	simply limult.  It is here for completeness; to eliminate
c	the overhead of a subroutine call, limult is often called
c	directly elsewhere in this file where lp_wpb is required.

c---  Arguments:
	real*8 lql_wp
	integer ndat, nbins, bins(nbins)

cf2py intent(in) ndat, nbins, bins

c---  Locals:
	real*8 limult

	lql_wp = limult(ndat, nbins, bins)

	return
	end

c-----------------------------------------------------------------------
	subroutine trap_absc(npts, nbins, absc)

c---	TRAPezoid rule ABSCissas
c
c	Calculate the phase points for trapezoid rule quadrature over
c	phase.

c---  Arguments:
	integer npts, nbins
	real*8 absc(npts)

cf2py intent(in) npts, nbins
cf2py intent(out) absc

c---  Locals:
	real*8 twopi
	parameter (twopi = 6.2831853071795862d0)

	del = twopi / nbins / (npts-1)
	do np=1, npts
	    absc(np) = (np-1)*del
	enddo

	return
	end

c-----------------------------------------------------------------------
	subroutine qlw_trap (ndat, lndfac, phases, nbins, phibins, bins,
     +                       npts, absc, lfac, qlike, chi2)

c---	QuasiLikelihood for W, TRAPezoid rule
c
c	Here "W" denotes the event phase parameters, e.g., 
c	w, wdot, wddot.
c
c	The likelihood for omega, conditional on the # of
c	bins, with phi and the f's marginalized away.
c
c	npts determines the number of points to use for the
c	trapezoid rule.
c
c	This is NOT normalized.
c
c	Because of symmetry, the integral over phi is done from
c	phi = 0 to phi = 2*pi/nbins, and then multiplied by nbins.
c
c	lndfac = log(ndat!), required as an argument to avoid
c	repeating calls to factln.
c
c	lfac is the logarithm of an arbitrary factor by which the
c	result is scaled, to facilitate computations that might
c	otherwise be out of the range of real*8 arithmetic.
c	The log of the actual likelihood is log(qlike) - lfac.
c
c	qlike is the quasilikelihood.
c
c	Also return the phase-averaged chi**2.
c
c	phibins and bins are used as storage and declared
c	"inplace" for Python to avoid repeated memory allocation.
c
c	The "cb" commented lines would implement Phil Gregory's bin 
c	factors, using a slower version of fbindata not included here.

c---  Arguments:
	integer ndat, nbins, bins(nbins), npts
	real*8 lndfac, phases(ndat), phibins(nbins), absc(npts)
	real*8 lfac, qlike, chi2

cf2py intent(in) ndat, lndfac, phases, nbins, npts, absc, lfac
cf2py intent(out) qlike, chi2
cf2py intent(inplace) phibins, bins

c---  Locals:
	integer npmax
	real*8 twopi
	parameter (npmax = 200, twopi = 6.2831853071795862d0)
	integer np
	real*8 del, arg, dphi
cb	real*8 binf(nbmax), logS
	real*8 plimult, chisqr

c---  Do the integral.
	qlike = 0.
	chi2 = 0.
	del = absc(2) - absc(1)
	do 40 np=1, npts
	    call fbindata(ndat, phases, absc(np), nbins, phibins,
     *                    bins, dphi)
	    arg = plimult(ndat, nbins, bins) - lndfac + lfac
c	    print *, 'bins ', absc(np), bins(1), bins(2), arg-lfac
cb	    call bindata(tlo, thi, data, ndat, w, absc(np), nbins, bins,
cb     +                   binf, logS, dphi)
cb	    arg = plimult(ndat, nbins, bins) - lndfac + logS + lfac
	    if (np .eq. 1 .or. np .eq. npts) then
	        qlike = qlike + 0.5*exp(arg)
	        chi2 = chi2 + 0.5*chisqr(ndat,nbins,bins)
	    else
	        qlike = qlike + exp(arg)
	        chi2 = chi2 + chisqr(ndat,nbins,bins)
	    endif
40	continue
	qlike = qlike * nbins * del
	chi2 = chi2 * nbins * del / twopi

	return
	end

c-----------------------------------------------------------------------
	subroutine qlw_exact (ndat, lndfac, phases, nbins, phibins, 
     +                        bins, lfac, qlike, chi2)

c---	Conditional quasilikelihood for W (omega), Exact
c
c	The quasilikelihood for omega, conditional on the # of
c	bins, with phi and the f's marginalized away.
c
c	The integral is done exactly, requiring ~ndat+1 quadrature
c	points.  The number of points used is returned as npts.
c
c	This is NOT normalized.
c
c	Because of symmetry, the integral over phi is done from
c	phi = 0 to phi = 2*pi/nbins, and then multiplied by nbins.
c
c	lfac is the logarithm of an arbitrary factor by which the
c	result is scaled, to facilitate computations that might
c	otherwise be out of the range of real*8 arithmetic.
c	The log of the actual likelihood is log(qlike) - lfac.
c
c	Also return the phase-averaged chi**2.
c
c	The "cb" commented lines would implement Phil Gregory's bin 
c	factors, using a slower version of fbindata not included here.

c---  Arguments:
	integer ndat, nbins, bins(nbins)
	real*8 lndfac, phases(ndat), phibins(nbins)
	real*8 lfac, qlike, chi2

cf2py intent(in) ndat, lndfac, phases, nbins, lfac
cf2py intent(out) qlike, chi2
cf2py intent(inplace) phibins, bins

c---  Locals:
	real*8 twopi
	parameter (twopi = 6.2831853071795862d0)
	real*8 del, arg, dphi, phi
cb	real*8 binf(nbmax), logS
	real*8 plimult, chisqr
c	logical show
c	integer i

c	if (w .gt. 0.) then
c	    qlike = w*exp(-(w-15.)**2/2.)
c	    return
c	endif

c---  Do the integral.
c	show = abs(w-.90162) .lt. .00005
	qlike = 0.
	chi2 = 0.
	phi = 0.
	npts = 0
	del = twopi / nbins
	sum = 0.
40	continue
	    call fbindata(ndat, phases, phi, nbins, phibins,
     *                    bins, dphi)
c	if (show) then
c	    write(*,'(2(1pg13.4),6i4,1pg12.4)') w, phi,(bins(i),i=1,6),
c     +               chisqr(ndat,nbins,bins)
c	endif
	    arg = plimult(ndat, nbins, bins) - lndfac + lfac
cc	write(*,'(4(1pg12.3),i4)') w, phi, dphi, bins(1)
cc	write(*,'(a,3(1pg13.4))') '---> ', arg+lndfac-lfac, lfac, arg
cb	    call bindata(tlo, thi, data, ndat, w, phi, nbins, bins, 
cb     +                   binf, logS, dphi)
cc	write(*,'(4(1pg13.4),i4)') w, phi, dphi, logS, bins(1)
cb	    arg = plimult(ndat, nbins, bins) - lndfac + logS + lfac
	    npts = npts + 1
	    if (phi+dphi .lt. del) then
	        qlike = qlike + exp(arg)*dphi
	        chi2 = chi2 + chisqr(ndat,nbins,bins)*dphi
	        phi = phi + dphi
	        goto 40
	    endif
	qlike = qlike + exp(arg)*(del-phi)
	qlike = qlike * nbins
	chi2 = chi2 + chisqr(ndat,nbins,bins)*(del-phi)
	chi2 = chi2 * nbins / twopi
cc	write(*,'(a,1pg13.4)') '===> ',qlike

	return
	end

c-----------------------------------------------------------------------
	function lodds (nbins, ndat, intimult, lfac)

c---	Log ODDS ratio
c
c	The log of the odds ratio in favor of an nbins model over 
c	a constant model.  
c
c	intimult is the integral of the (phase-averaged) inverse multiplicity
c	times the frequency parameter (w, wdot, wddot) prior over the
c	frequency parameter space.
c
c	lfac should be the same log factor passed to cp_w.

c---  Arguments:
	real*8 lodds, intimult, lfac
	integer nbins, ndat

cf2py intent(in) nbins, ndat, intimult, lfac

c---  Locals:
	real*8 twopi, arg, curf, factln
	integer curnb, curnd
	parameter (twopi = 6.2831853071795862d0)
	data curnb, curnd /-1, -1/
	save curf, curnb, curnd

	if (nbins .ne. curnb .or. ndat .ne. curnd) then
	    curnb = nbins
	    curnd = ndat
	    arg = nbins
	    curf = factln(nbins-1) + ndat * log(arg)
	    curf = curf - factln(ndat+nbins-1)
	    curf = curf + factln(ndat)
	endif

	lodds = curf + log(intimult/twopi) - lfac

	return
	end

c-----------------------------------------------------------------------
	subroutine rate_wb (ndat, lndfac, phases, nbins, phibins, 
     +                      bins, lfac, nt, rphases, r, r2)

c---	"RATE estimate conditional on W and Bin #"
c
c	The rate is estimated at the input phases in rphases;
c	r and r2 will contain the estimate and uncertainty.
c	The input phases should be in [0,2*pi].
c
c	lfac is the logarithm of an arbitrary factor by which the
c	result is scaled, to facilitate computations that might
c	otherwise be out of the range of real*8 arithmetic.
c
c	The phase integral is done exactly.
c
c	phibins and bins are used as storage and declared
c	"inplace" for Python to avoid repeated memory allocation.
c
c	The "cb" commented lines would implement Phil Gregory's bin 
c	factors, using a slower version of fbindata not included here.

c---  Arguments:
	integer ndat, nbins, bins(nbins), nt
	real*8 lndfac, phases(ndat), phibins(nbins), lfac
	real*8 rphases(nt), r(nt), r2(nt)

cf2py intent(in) ndat, lndfac, phases, nbins
cf2py intent(in) lfac, nt, rphases
cf2py intent(out) r, r2
cf2py intent(inplace) phibins, bins

c---  Locals:
	integer n, j
	real*8 del, arg, dphi, phi, norm, rc, r2c
cb	real*8 binf(nbmax), logS
	real*8 plimult, one, two, twopi
	parameter (twopi = 6.2831853071795862d0)
	parameter (one = 1.d0, two = 2.d0)

c---  Do the integral over phi exactly.
	norm = 0.
	do 10 n=1, nt
	    r(n) = 0.
	    r2(n) = 0.
10	continue
	phi = 0.
	npts = 0
	del = twopi / nbins
40	continue
cb	    call bindata(tlo, thi, data, ndat, w, phi, nbins, bins, 
cb     +                   binf, logS, dphi)
cb	    arg = plimult(ndat, nbins, bins) - lndfac + logS + lfac
	    call fbindata(ndat, phases, phi, nbins, phibins,
     *                    bins, dphi)
	    arg = plimult(ndat, nbins, bins) - lndfac + lfac
	    arg = exp(arg)
	    norm = norm + arg*dphi
	    npts = npts + 1
	    do 60 n=1, nt
	        j = one + nbins*mod(rphases(n) + phi, twopi)/twopi
	        rc = nbins * (bins(j)+one)/(ndat+nbins)
	        r2c = rc * nbins * (bins(j)+two)/(ndat+nbins+one)
	        if (phi+dphi .lt. del) then
	            r(n) = r(n) + rc*arg*dphi
	            r2(n) = r2(n) + r2c*arg*dphi
	        else
	            r(n) = r(n) + rc*arg*(del-phi)
	            r2(n) = r2(n) + r2c*arg*(del-phi)
	        endif
60	    continue
	    if (phi+dphi .lt. del) then
	        phi = phi + dphi
	        goto 40
	    endif
	do 100 n=1, nt
	    r(n) = r(n)/norm
	    r2(n) = r2(n)/norm
c	    write(*,'(i3,5(1pg13.4))') nbins,rphases(n),r(n),r2(n),norm,lfac
100	continue

	return
	end

c-----------------------------------------------------------------------
	function chisqr (ndat, nbins, bins)

c---  	"CHI SQuaRed"
c
c	Return the value of the traditional epoch-folding Pearson's
c	chi**2 value for the binned data.  This assumes no
c	"dead time" in any bins.

c+++  Arguments:
	real*8 chisqr
	integer ndat, nbins, bins(nbins)

cf2py intent(in) ndat, nbins, bins

c+++  Locals:
	integer i
	real*8 expect

	expect = ndat
	expect = expect / nbins
	chisqr = 0.
	do 20 i=1, nbins
	    chisqr = chisqr + (bins(i) - expect)**2
20	continue
	chisqr = chisqr / expect

	return
	end

c-----------------------------------------------------------------------
	function factln(n)

c	Logarithm of the factorial of n, adapted from Numerical Recipes.

	integer n
	real*8 factln

cf2py intent(in) n

	real*8 x, a(100), gammln
	data A/100*-1.d0/
	save a

	if (n .LT. 0) then
          pause 'negative factorial'
	else if (n .LE. 99) then
        if (a(n+1) .lt. 0.d0) then
            x = n + 1.d0
            a(n+1) = gammln(x)
        endif
        factln = a(n+1)
	else
	    x = n + 1.d0
            factln = gammln(x)
	endif

	return
	end

c-----------------------------------------------------------------------
        function gammln(xx)

c---    Adapted from 'Numerical Recipes,' this function calculates
c       the gamma function, used for calculating Gaussian deviates.

cf2py	intent(in) xx

        real*8 gammln
        real*8 xx, cof(6), stp, half, one, fpf, x, tmp, ser
        integer j
        data cof, stp /76.18009173d0, -86.50532033d0, 24.01409822d0,
     1      -1.231739516d0, 0.120858003d-2, -0.536382d-5,
     2      2.50662827465d0/
        data half, one, fpf /0.5d0, 1.0d0, 5.5d0/
        save cof, stp, half, one, fpf

        x = xx - one
        tmp = x + fpf
        tmp = (x + half)*log(tmp) - tmp
        ser = one
        do 10 j=1, 6
            x = x + one
            ser = ser + cof(j) / x
10      continue
        gammln = tmp + log(stp*ser)

        return
        end

