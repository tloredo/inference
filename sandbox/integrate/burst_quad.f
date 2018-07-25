c-----------------------------------------------------------------------
	subroutine set_bqpts (nF, np, na)

c---	"SET Burst Quadrature PoinTS"
c
c	Sets the # of points for quadrature of the
c	event rate times the flux & position likelihoods.
c
c	Arguments:
c	nF is the # of points for the flux integral.
c	np is the # of points for the polar direction integral.
c	na is the # of points for the azimuthal direction integral.
c
c	Total # of points used per burst is nf*np*na.

c+++  Arguments:
	integer nF, np, na

c+++  Globals:
	real*8 zero, twopi
	parameter (zero = 0., twopi = 6.283185307179586d0)
	integer nqmax
	parameter (nqmax = 50)
	integer nFqpts
	real*8 Fabsc(nqmax), Fwts(nqmax)
	common /Fqcb/ Fabsc, Fwts, nFqpts
	integer nmu, nphi
	real*8 pabsc(nqmax), pwts(nqmax)
	common /dqcb/ pabsc, pwts, nmu, nphi
	include 'data.cb'

c---  Get GH points for F.
	if (nF .gt. nqmax) call quit('Too many F quad points!')
	nFqpts = nF
	if (nF .eq. 1) then
	    Fabsc(1) = 0.
	    Fwts(1) = 1.
	else
	    call gauher(Fabsc, Fwts, nFqpts)
	endif

c---  Get GL points for phi.
	if (na .gt. nqmax) call quit('Too many azimuth quad points!')
	nphi = na
	if (na .eq. 1) then
	    pabsc(1) = 0.
	    pwts(1) = twopi
	else
	    call gauleg(zero, twopi, pabsc, pwts, nphi)
	endif

c---  Set the # of polar points; the abscissas & weights must be
c---  set for each burst with set_bquad.
	if (np .gt. nqmax) call quit('Too many polar quad points!')
	nmu = np

c---  If there is data, set the mu quadrature rules.
	if (nb .gt. 0) then
	    call set_bqmu (nb, dsigs, kvals, absc, wts, mumax)
	endif

c---  Note that things must be recalculated.
	call dmchange

	return
	end

c-----------------------------------------------------------------------
	subroutine set_bqmu (nb, dsigs, kvals, absc, wts, npmax)

c---	"SET Burst Quadrature MU"
c
c	Sets the abscissas & weights for polar integrals over
c	the individual burst likelihood functions.
c
c	Arguments:
c	nb is the # of bursts
c	dsigs = "sigmas" reported in 1B catalog (68.3% confidence)
c	npmax is the 1st physical dimension of absc.
c
c	Returns:
c	kvals stores the kappa values for the Fisher dist'n
c	    modeling the direction uncertainty for burst i.
c	absc(m,i) is the m'th abscissa for the polar integral 
c	    for burst i.
c	wts(m,i) is the m'th weight for the polar integral for burst i.

c+++  Arguments:
	integer nb, npmax
	real*8 dsigs(*), kvals(*), absc(npmax,*), wts(npmax,*)

c+++  Globals:
	real*8 zero, twopi
	parameter (zero = 0., twopi = 6.283185307179586d0)
	integer nqmax
	parameter (nqmax = 50)
	integer nFqpts
	real*8 Fabsc(nqmax), Fwts(nqmax)
	common /Fqcb/ Fabsc, Fwts, nFqpts
	integer nmu, nphi
	real*8 pabsc(nqmax), pwts(nqmax)
	common /dqcb/ pabsc, pwts, nmu, nphi

c+++  Locals:
	real*8 one, minus, C, C1, eps
	integer nrmax
	parameter (one = 1., minus = -1., C = 0.683, C1 = 1.-C)
	parameter (eps = 1.e-5, nrmax = 25)
	integer i, j, iter
	real*8 gla(nqmax), glw(nqmax), lo, fac, off
c	real*8 kold, k, f, csig, ek, ekc, rek, rat, sh
	real*8 kold, k, f, csig, arg, ekc, ekk, num, den

c---  Get GL points for mu for each burst.  First get GL points over
c---  [-1,1]; then transform them for each burst.  To transform
c---  them, first find k for that burst using Newton-Raphson.  Use
c---  the Gaussian (small angle) approximation as a starting guess.
	if (nmu .gt. npmax) call quit('Too many burst mu points!')
	if (nmu .eq. 1) then
	    gla(1) = 0.
	    glw(1) = 1.
	else
	    call gauleg (minus, one, gla, glw, nmu)
	endif
	do 30 i=1, nb
	    kold = 2.297707 / dsigs(i)**2
	    csig = cos(dsigs(i))
	    iter = 0
10	    continue
	    if (iter .ge. nrmax) call quit('Too many iter for kappa!')

c---  This commented code fails for small dsigs due to the exp.
c	    ek = exp(kold)
c	    ekc = exp(kold*csig)
c	    rek = 1. / ek
c	    sh = ek - rek
c	    rat = (ek - ekc) / (ek - rek)
c	    f = rat - C
c	    k = kold - 
c     *          f / ((ek-csig*ekc)/sh - ((ek-ekc)/sh)*((ek+rek)/sh))

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
	    k = kold -
     *          f * den / (2.*num*ekk/den - arg*ekc)

	    iter = iter + 1
	    if (abs(k-kold)/k .gt. eps) then
	        kold = k
	        goto 10
	    endif
	    kvals(i) = k
	    if (k .lt. 250.) then
	        lo = exp(-2.*k)
	    else
	        lo = 0.
	    endif
            fac = 0.5 * (1. - lo)
            off = 0.5d0 * (1. + lo)
            do 20 j=1, nmu
                absc(j,i) = (log(fac*gla(j) + off) / k) + 1.
                wts(j,i) = 0.5 * glw(j) / twopi
20	    continue
30	continue

	return
	end

c-----------------------------------------------------------------------
	subroutine F_quad (F, sig, dR)

c---	"Flux QUADrature"
c
c	Returns the quadrature of dN/dFdO*Gaussian for a burst 
c	with detected peak photon # flux F, with uncertainty sig,
c	for an isotropic model.
c
c	The rate function used is dR_dF, which may interpolate
c	from a grid.

c+++  Arguments:
	real*8 F, sig, dR

c+++  Globals:
	integer nqmax
	parameter (nqmax = 50)
	integer nFqpts
	real*8 Fabsc(nqmax), Fwts(nqmax)
	common /Fqcb/ Fabsc, Fwts, nFqpts
	real*8 etabar
	common /Fqetab/ etabar

c+++  Locals:
	real*8 rtpi, rt2, pi4, fac
	parameter (rtpi = 1.772453850905516d0, rt2= 1.414213562373095d0)
	parameter (pi4 = 12.56637061435917d0, fac = rtpi*pi4)
	integer i
	real*8 FF

c+++  Functions:
	real*8 dR_dF

c---  Here we do the outermost integral, over F.
c---  For nFqpts = 1, it's easy!
	if (nFqpts .eq. 1) then
	    dR = dR_dF (F)
	    return
	endif

c---  Otherwise, it's still pretty easy...  just sum the quadrature,
c---  noting that the Gaussian is included in the wts.
c---  Omit points from negative F, if they arise.
	dR = 0.
	do 20 i=1, nFqpts
	    FF = sig*rt2*Fabsc(i) + F
	    if (FF .le. 0.) goto 20
	    dR = dR + Fwts(i)*dR_dF(FF)
20	continue

c---  We divide dR by rtpi from scaling the quadrature rule, and
c---  by an additional 4*pi since we integrated dR/dF, but need
c---  the integral of dR/dFdO.
	dR = dR / fac

	return
	end

c-----------------------------------------------------------------------
	subroutine d_quad (F, mu_b, l, k, absc, wts, dR)

c---	"Direction QUADrature"
c
c	Returns the direction part of the Fisher-weighted integral
c	of the burst rate.  Gauss-Legendre quadrature is used.

c+++  Arguments:
	real*8 F, mu_b, l, k, absc(*), wts(*), dR

c+++  Globals:
	integer nqmax
	parameter (nqmax = 50)
	integer nmu, nphi
	real*8 pabsc(nqmax), pwts(nqmax)
	common /dqcb/ pabsc, pwts, nmu, nphi

c+++  Locals:
	real*8 one, pi, rt2, C_F
	parameter (one = 1., pi = 3.141592653589793d0)
	parameter (rt2 = 1.414213562373095d0)
	parameter (C_F = 0.1125395395196383d0)
	integer i, j
	real*8 s_b, mu, sn, phi, mu_v, s_v, l_v, ddR

c+++  Function:
	real*8 dR_dFdO

c---  For nmu = 1, it's easy!
	if (nmu .eq. 1) then
	    dR = dR_dFdO(F, mu_b, l)
	    return
	endif

c---  Otherwise, it's a little harder.  We work in terms
c---  of a coordinate system with z axis along (l,b); the (x,y)
c---  axes are arbitrary, since the Fisher dist'n is azimuthally
c---  symmetric.
	dR = 0.
	s_b = sqrt(one - mu_b*mu_b)
	do 40 i=1, nmu
	    mu = absc(i)
	    sn = sqrt(one - mu*mu)

c---  Here is the azimuthal part, which is easy once we transform
c---  tau to mu and rotate to (l,b).  For the rotation, we take
c---  the GNP to be in the phi=0 plane.
	    ddR = 0.
	    do 20 j=1, nphi
	        phi = pabsc(j)
	        call sunrotate (mu_b, s_b, l, mu, sn, phi, 
     +                          mu_v, s_v, l_v)
	        ddR = ddR + pwts(j)*dR_dFdO(F, mu_v, l_v)
20	    continue
	    dR = dR + wts(i)*ddR
40	continue

	return
	end

c-----------------------------------------------------------------------
	subroutine Flb_quad (F, sig, l, b, k, absc, wts, dR)

c---	"Flux-l-b QUADrature"
c
c	Returns the quadrature of dN/dCdO*Gaussian for a burst 
c	with detected peak photon # flux F, with uncertainty sig,
c	and position (l,b), with uncertainty k.
c	Gauss-Hermite quadrature is used for the F quadrature;
c	Gauss-Legendre is used for the direction quadratures.

c+++  Arguments:
	real*8 F, sig, l, b, k, absc(*), wts(*), dR

c+++  Globals:
	integer nqmax
	parameter (nqmax = 50)
	integer nFqpts
	real*8 Fabsc(nqmax), Fwts(nqmax)
	common /Fqcb/ Fabsc, Fwts, nFqpts

c+++  Locals:
	real*8 rtpi, rt2
	parameter (rtpi = 1.772453850905516d0, rt2= 1.414213562373095d0)
	integer i
	real*8 FF, val, mu_b

c---  Here we do the outermost integral, over F.
c---  For nFqpts = 1, it's easy!
	mu_b = sin(b)
	if (nFqpts .eq. 1) then
	    call d_quad (FF, mu_b, l, k, absc, wts, dR)
	    return
	endif

c---  Otherwise, it's still pretty easy...  just sum the quadrature,
c---  noting that the Fisher dist'n factor is included in the wts.
c---  Omit points from negative F, if they arise.
	dR = 0.
	do 20 i=1, nFqpts
	    FF = sig*rt2*Fabsc(i) + F
	    if (FF .le. 0.) goto 20
	    call d_quad (FF, mu_b, l, k, absc, wts, val)
	    dR = dR + Fwts(i)*val
20	continue
	dR = dR / rtpi

	return
	end

