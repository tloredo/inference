c=======================================================================
c	Subprograms for output analysis of time series produced
c	by Markov Chain Monte Carlo (MCMC) methods.
c
c	June 1993  TJL
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	subroutine odiag(ndat, data, del, mean, arsig, spsig, ord, clen,
     *                  wrk1, wrk2, wrk3, ok)

c---	"Output DIAGnostics"
c
c	Given the time series of length ndat in data, return the
c	number of times, del, to delete so that the remaining
c	time series is stationary.  Return the mean and standard
c	deviation for that time series.  Two standard deviation
c	estimates are returned: arsig from an AR model, and
c	spsig from fitting the low frequency power spectrum.
c	Also return the correlation length from the AR model (clen)
c	and the order used (ord).
c
c	wrk2,3 are real*8 workspaces, at least nt elements 
c	in length.  wrk1 must be at least the power of 2 >= nt 
c	in length.
c
c	ok = .false. if the time series is not stationary over
c	at least half its length, in which case the mean and
c	std dev'ns are not estimated.
c
c	del is determined by standardizing parts of the time
c	series to an approximate Brownian bridge, and using the
c	Cramer vonMises statistic to test for deviations from
c	the bridge (which is the dist'n for a *stationary* series).
c	See Ripley (1987), Heidelberg & Welch (1983), and
c	Schruben (1983) for details.  We test at the 10% false
c	rejection level, as suggested by H&W.
c
c	The AR order for arsig is determined using Schwarz's BIC.
c	The AR coefficients are estimated with the Burg algorithm.
c	See Ripley (1987), Fishman (1973, 1978), and Numerical
c	Recipes for details.
c
c	The spectral standard deviation is found by fitting a 
c	quadratic to the log power at the first 25 frequencies
c	of the DFT.  See Heidelberg & Welch (1981) for details.

c---  Arguments:
	integer ndat, del, ord
	real*8 data(*), wrk1(*), wrk2(*), wrk3(*)
	real*8 mean, arsig, spsig, clen
	logical ok

c---  Locals:
	integer dtrials, maxo
	real*8 Ccrit
	parameter (maxo = 30, dtrials = 6, Ccrit = 0.347)
	real*8 tfac(dtrials)
	save tfac
	data tfac /0., .1, .2, .3, .4, .5/
	integer t, th, tlo, i, n, nt, del, nu
	real*8 var, dt, CvM, fac, sum

c---  First look for the initial transient.  Begin by estimating
c---  the variance from the second half of the series by
c---  subtracting off its mean and fitting an AR model.  Note that
c---  we want the variance of a single point, not the variance for
c---  the mean, so we multiply the latter by the length of the series.
	th = ndat / 2
	tlo = ndat - th + 1
	mean = 0.
	do 20 t=tlo, ndat
	    mean = mean + data(t)
20	continue
	mean = mean / th
	do 40 t=tlo, ndat
	    wrk1(t-tlo+1) = data(t) - mean
40	continue
	call bestar (th, wrk1, ord, var, clen, wrk2, wrk3)
	var = th * var

c===  Check for stationarity, deleting increasing fractions of the data.
	do 120 i=1, dtrials
	    del = tfac(i) * ndat
	    tlo = del + 1
	    n = ndat - del
	    fac = 1. / (var*n)
	    dt = 1. / n

c===  First find the mean of the undeleted part.
	    mean = 0.
	    do 60 t=tlo, ndat
	        mean = mean + data(t)
60	    continue
	    mean = mean / n

c===  Now the standardized time series, for t in [0,1], is 
c===  fac * (sum from n*t - n*t*mean).  Integrate its square
c===  to get the CvM statistic.  Note that it vanishes at
c===  t=0 and t=1, so we don't have to worry about the
c===  endpoints of the trapezoid rule.
	    CvM = 0.
	    do 100 nt=1, n-1
	        sum = 0.
	        do 80 t=1, nt
	            sum = sum + data(tlo+t-1)
80	        continue
	        CvM = CvM + fac*(sum - nt*mean)**2
100	    continue
	    CvM = dt * CvM
c	    write(*,'(a,i5,1pg12.4)') 'del, CvM: ',del, CvM

c===  If we pass the test, leave the loop; otherwise, try the
c===  next deletion.
	    if (CvM .lt. Ccrit) goto 200
120	continue

c---  If we made it here, the time series is not stationary over
c---  any appreciable part.
	ok = .false.
	del = ndat
	return

c---  Here's where we go if we've identified a stationary part.
c---  The mean has already been estimated; now estimate the
c---  std deviation and corr length.
200	continue

c---  Here's the spectral estimate.
	call specvar (n, data(tlo), spsig, nu, wrk1)
	spsig = sqrt(spsig)

c---  Now get the AR estimate, which requires we subtract the mean.
	do 220 t=1, n
	    wrk1(t) = data(t+del) - mean
220	continue
	call bestar (n, wrk1, ord, arsig, clen, wrk2, wrk3)
	arsig = sqrt(arsig)
	ok = .true.

	return
	end

c-----------------------------------------------------------------------
	subroutine bestar (ndat, data, ord, var, clen, wrk, wrks)

c---	"BEST AutoRegressive model"
c
c	Find the best AR model to fit the time series of length ndat
c	in data.  Return the order, estimated variance, and correlation
c	length.  wrk and wrks are workspaces at least (ndat-1) in
c	size.
c
c	The AIC is used to determine what order is "best".  I'd use
c	the BIC, except the samples are correlated, and we're using
c	this routine to measure the correlation!
c
c	Orders up to 30 are examined.

c---  Arguments:
	integer ndat, ord
	real*8 data(*), var, clen, wrk(*), wrks(*)

c---  Locals:
	integer maxo, csiz
	parameter (maxo = 30, csiz = (maxo*(maxo+1))/2)
	integer o, s, cur, off, omax
	real*8 ln, sum, AIC, AIClo, v, c
	real*8 vars(maxo), coef(csiz)

	ln = ndat
	ln = log(ln)
	omax = min(ndat-1,maxo)
	do 100 o=1, omax
	    cur = o - 1
	    call nxtarc (ndat, data, cur, vars, coef, wrk, wrks)
	    AIC = ndat*log(vars(o)) + 2.*o
c	    BIC = ndat*log(vars(o)) + o*ln
	    sum = 0.
	    off = (o*(o-1))/2
	    do 120 s=1, o
	        sum = sum + coef(off+s)
120	    continue
	    c = 1. / (1.-sum)**2
	    v = vars(o) / (ndat/c)
	    if (o .eq. 1) then
	        ord = o
	        AIClo = AIC
	        var = v
	        clen = c
	    else
	        if (AIC .lt. AIClo) then
	            AIClo =AIC
	            ord = o
	            var = v
	            clen = c
	        endif
	    endif
c	    write(*,80) o, sqrt(v), c, AIC
c80	    format('  Order: ',i4,'  sig = ',1pg12.4,
c     *        '  Clen = ',1pg12.4,'  AIC = ',1pg12.4)
100	continue

	return
	end

c-----------------------------------------------------------------------
	subroutine nxtarc (ndat, data, cur, var, coef, wrk, wrks)

c---	"NeXT Auto Regressive Coefficients"
c
c	Assuming the autoregressive coefficients have been
c	calculated for orders up to cur, as well as estimates
c	of the noise variance, calculate coefficients and
c	a noise variance estimate for the next order using
c	Burg's algorithm.
c
c	var(p) is the variance estimate for order p
c	coef(p*(p-1)/2+s) is the coefficient for lag s
c	         for order p
c	cwrk is a workspace for coefficients, at least cur+1 in size
c	wrk and wrks are workspaces at least (ndat-1) in size
c
c	NOTE:  Do not alter the workspaces between calls!! 
c	Successive calls use the results in vars, coefs, and
c	the workspaces.
c
c	The BIC for a particular order is 
c	          ndat*log(var(ord)) + log(ndat)*ord
c	A criterion for a good order is that which minimizes BIC.
c	There is a quasibayes justification for this.  However, it
c	isn't really valid here, because it assumes independent samples.
c
c	A heuristic alternative is the AIC,
c	          ndat*log(var(ord)) + 2.*ord
c	Ripley (1987) recommends this instead of the BIC.
c	Fishman (1973, 1978) offers significance tests for the order.
c
c	This routine was adapted from MEMCOF in Press et al,
c	Numerical Recipes (1993).

c---  Arguments:
	integer ndat, cur
	real*8 data(*), var(*), coef(*), wrk(*), wrks(*)

c---  Locals:
	integer t, ord, off, poff, s, r
	real*8 numer, denom

c===  Treat the cur=0 (order 1) case specially.
	if (cur .eq. 0) then

c---  Initialize workspaces for the recursion, and the var estimate.
	    var(1) = 0.
	    do 20 t=1, ndat
	        var(1) = var(1) + data(t)**2
20	    continue
	    var(1) = var(1) / ndat
	    wrk(1) = data(1)
	    wrks(ndat-1) = data(ndat)
	    do 40 t=2, ndat-1
	        wrk(t) = data(t)
	        wrks(t-1) = data(t)
40	    continue

c---  Get estimates for the first order
	    numer = 0.
	    denom = 0.
	    do 60 t=1, ndat-1
	        numer = numer + wrk(t)*wrks(t)
	        denom = denom + wrk(t)**2 + wrks(t)**2
60	    continue
	    coef(1) = 2. * numer / denom
	    var(1) = var(1) * (1. - coef(1)**2)
	    return
	endif

c===  For cur>0, use Burg's recursive procedure.
	ord = cur + 1
	off = (ord*(ord-1))/2
	poff = (cur*(cur-1)/2)

c---  Start with vars = previous vars, and the first cur
c---  coefficients for the next order = to those from previous order.
	var(ord) = var(cur)
	do 80 s=1, cur
	    coef(off+s) = coef(poff+s)
80	continue

c---  Now adjust them to get results for order (cur+1).
	numer = 0.
	denom = 0.
	do 100 t=1, ndat-cur-1
	    wrk(t) = wrk(t) - coef(poff+cur)*wrks(t)
	    wrks(t) = wrks(t+1) - coef(poff+cur)*wrk(t+1)
100	continue
	do 120 t=1, ndat-ord
	    numer = numer + wrk(t)*wrks(t)
	    denom = denom + wrk(t)**2 + wrks(t)**2
120	continue
	coef(off+s) = 2. * numer / denom
	var(ord) = var(ord) * (1. - coef(off+s)**2)
	do 140 r=1, s-1
	    coef(off+r) = coef(poff+r) - coef(off+s)*coef(poff+s-r)
140	continue

	return
	end

c-----------------------------------------------------------------------
	subroutine specvar (ndat, data, var, nu, wrk)

c---	"SPECtral estimate of VARiance"
c
c	Estimate the variance to associate with the sample mean
c	of the correlated time series of length ndat in data(*).
c	The estimate is returned as var, and the # of degrees
c	of freedom to associate with it as nu.
c
c	wrk is a real*8 workspace of length >= ndat that is a
c	power of 2.
c
c	We (somewhat arbitrarily) require ndat >= 200.  It should
c	probably be considerably larger!  If ndat < 200, return
c	var = nu = 0.
c
c	If the power at a low frequency is <0, return var=0, nu=-1.
c
c	The method is due to Heidelber & Welch (1981), and involves
c	fitting a quadratic to the log of the sum of the power
c	spectrum at two adjacent frequencies (the sum reduces
c	skewness), and correcting for bias.

c---  Arguments:
	integer ndat, nu
	real*8 data(*), var, wrk(*)

c---  Locals:
	real*8 log2, J0, psvar
	integer nf, nf2
	parameter (log2 = 0.6931471805599453d0, nf = 25, nf2 = 2.*nf)
	parameter (J0 = 0.270, psvar = 0.645)
	integer npow, nft, t, n, n1, flist(3)
	real*8 J(nf), f(nf), sig(nf), coef(3), cov(3,3)
	real*8 rn, chi, fac, mean

c---  Functions:
	external poly

c---  Make sure the data are long enough.
	if (ndat .lt. 4*nf) then
	    var = 0.
	    nu = 0
	    return
	endif

c---  Find how long the FFT will be.  Note that the realft routine 
c---  will need *half* the length.
	rn = ndat
	npow = log(rn) / log2
	nft = 2**npow
	if (nft .lt. ndat) nft = 2 * nft

c---  Copy the data to the workspace, where it will be padded and
c---  replaced with its DFT.  Find the mean on the way, for padding.
	mean = 0.
	do 20 t=1, ndat
	    wrk(t) = data(t)
	    mean = mean + data(t)
20	continue
	mean = mean / ndat

c---  Pad it with its mean.
	do 40 t=ndat+1, nft
	    wrk(t) = mean
40	continue

c---  Do an FFT of the prewhitened residuals.
	call realft(wrk, nft, 1)

c---  Collect the first 2*nf powers

c---  Collect nf logs of sums of powers at adjacent frequencies.  
c---  Note that the zero frequency DFT is not used.
c---  Along the way, build the nf-by-3 matrix for a least squares fit 
c---  to a quadratic for the first nf sums.
	do 60 n=1, nf
	    n1 = 4*(n-1) + 3
	    J(n) = 0.5 * (wrk(n1)**2   + wrk(n1+1)**2 + 
     +                    wrk(n1+2)**2 + wrk(n1+3)**2) / ndat
	    if (J(n) .le. 0.) then
	        var = 0.
	        nu = -1
	        return
	    endif
	    J(n) = log(J(n)) + J0
	    f(n) = 0.5 * (4.*n-1) / ndat
	    sig(n) = 1.
60	continue

c---  Do the fit.
	flist(1) = 1
	flist(2) = 2
	flist(3) = 3
	call lfit(f, J, sig, nf, coef, flist, 3, cov, 3, chi, poly)

c---  To get the estimate at f=0, we don't just use the intercept
c---  (1st coef), because we have fit the log.  We have to correct
c---  for the expectation of a log normal mean.
	fac = exp(-0.5*psvar*cov(1,1))
	var = fac * exp(coef(1)) / ndat
	nu = nint( 2./(1./fac - 1.) )

	return
	end

c-----------------------------------------------------------------------
	subroutine poly (f, vals, num)

c---	"POLYnomial"
c
c	Return values of the terms of the quadratic at f.

c---  Arguments:
	integer num
	real*8 f, vals(*)

	vals(1) = 1.
	vals(2) = f
	vals(3) = f*f

	return
	end

