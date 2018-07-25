c-----------------------------------------------------------------------
c	Subprograms for calculating the Cmax dist'n of bursters
c	in a matter-dominated Friedman cosmology with constant
c	burst rate per unit (comoving) coordinate volume, and
c	power-law "luminosity" function and spectrum.
c
c	22 Apr 92 TJL
c-----------------------------------------------------------------------
	subroutine set_HO (hh, O)

c---	"SET H0 and Omega"
c
c	Set the current value of the Hubble constant (in units of
c	100 km / s / Mpc) and Omega.

c+++  Arguments:
	real*8 hh, O

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut
	data h/0./, Omega/0./, el/0./, eu/0./, e1/0./, e2/0./

	h = hh
	Omega = O
	q0 = O/2.

	return
	end

c-----------------------------------------------------------------------
	subroutine set_lfun (power, nu_u, ratio)

c---	"SET Luminosity FUNction"
c
c	nu_u = Nu * H0**2 * F_spec(0) / 4*pi*c**2
c
c	Note:  set_HO, set_nspec and set_det must be called before this!

c+++  Arguments:
	real*8 power, nu_u, ratio

c+++  Locals:
	real*8 c100, pi4, zero
	parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
	parameter (zero = 0.)
	real*8 F_spec

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

	if (h*Omega .eq. 0. .or. (el .eq. eu) .or. (e1 .eq. e2)) 
     +            pause 'Premature set_lfun!'
	p = power
	Nu = nu_u * pi4 * c100**2 / (h**2 * F_spec(zero))
	Nl = Nu / ratio
	if (p .ne. 1.) then
	    A = (p-1) / ( ratio**(p-1) - 1. )
	else
	    A = 1. / log(ratio)
	endif

	return
	end

c-----------------------------------------------------------------------
	subroutine set_blfun (sindx, nu_u, ratio)

c---	"SET Beamed Luminosity FUNction"
c
c	Set the luminosity function from beaming of a source with
c	spectrum ~ E**(-sindx).  This luminosity function is not
c	normalized.
c
c	nu_u = Nu * H0**2 * F_spec(0) / 4*pi*c**2
c
c	Note:  set_HO, set_nspec and set_det must be called before this!

c+++  Arguments:
	real*8 sindx, nu_u, ratio

c+++  Locals:
	real*8 c100, pi4, zero
	parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
	parameter (zero = 0.)
	real*8 F_spec

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

	if (h*Omega .eq. 0. .or. (el .eq. eu) .or. (e1 .eq. e2)) 
     +            pause 'Premature set_blfun!'
	p = (3.+sindx) / (2.+sindx)
	Nu = nu_u * pi4 * c100**2 / (h**2 * F_spec(zero))
	Nl = Nu / ratio
	A = 1.

	return
	end

c-----------------------------------------------------------------------
	subroutine set_dens (n_0, power)

c---	"SET number DENSity"
c
c	n_0 is in units of bursts / sec / Mpc**3

c+++  Arguments:
	real*8 n_0, power

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

	n0 = n_0
	beta = power

	return
	end

c-----------------------------------------------------------------------
	subroutine set_spec (power, lo, hi, dlo, dhi)

c---	"SET Number SPECtrum"
c
c	Set the spectrum of bursts, and the spectral range of the
c	detector.

c+++  Arguments:
	real*8 power, lo, hi, dlo, dhi

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

	e1 = dlo
	e2 = dhi
	alpha = power
	el = lo
	eu = hi
	if (alpha .ne. 1.) then
	    K = (e2/eu)**(1.-alpha) / ( (eu/el)**(alpha-1) - 1. )
	else
	    K = 1. / log(eu/el)
	endif
	if (e1 .gt. 0. .and. el .gt. e1) 
     *      write(*,'(a)') 'WARNING:  el > e1!'

	return
	end

c-----------------------------------------------------------------------
	subroutine set_zcut (zc)

c---	"SET Z CUToff"

c+++  Arguments:
	real*8 zc

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

	zcut = zc

	return
	end

c-----------------------------------------------------------------------
	function cdR_dOdC (C, tol)

c---	"Cosmological dR / dOmega dCmax"

c+++  Arguments:
	real*8 cdR_dOdC, C, tol

c+++  Locals:
	real*8 c100, pi4, zero
	parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
	parameter (zero = 0.)
	real*8 z1, z2, zm, z3, tz1, tz2, ratio, pow1, dR
	real*8 dR_intgd
	external dR_intgd

c+++  Globals:
	real*8 gC, pow, rpow1
	common /ccb/ gC, pow, rpow1
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

c---  Find the z range to integrate over.
	gC = C
	call z_range(C, z1, z2, zm)
	if (z2 .le. z1) then
	    cdR_dOdC = 0.
	    return
	endif
	if (zm .gt. z1 .and. zm .lt. z2) then
	    z3 = z2
	    z2 = zm
	else
	    z3 = 0.
	endif
c	write(9,'(a,3(1pg12.4))') 'z: ',z1,z2,z3

c---  Find the power law transformation for the z integral & do it.
	pow = 0.
	rpow1 = 1.
	ratio = dR_intgd(z2) / dR_intgd(z1)
	pow = - log(ratio) / log(z2/z1)
	if (pow .eq. 1.) pow = 0.
	pow1 = 1. - pow
	rpow1 =1. / pow1
	tz1 = z1**pow1
	tz2 = z2**pow1
	call qromb(dR_intgd, tz1, tz2, cdR_dOdC, tol)
	cdR_dOdC = cdR_dOdC / pow1
c	write(9,'(a,5(1pg12.4))') 'I1: ',C,dR_intgd(z1),
c     +      dR_intgd(z2),pow,cdR_dOdC

c---  Do any second integral, if required.
	if (z3 .gt. z2) then
	    pow = 0.
	    rpow1 = 1.
	    ratio = dR_intgd(z3) / dR_intgd(z2)
	    pow = - log(ratio) / log(z3/z2)
	    if (pow .eq. 1.) pow = 0.
	    pow1 = 1. - pow
	    rpow1 =1. / pow1
	    tz1 = z2**pow1
	    tz2 = z3**pow1
	    call qromb(dR_intgd, tz1, tz2, dR, tol)
	    dR = dR / pow1
c	write(9,'(a,5(1pg12.4))') 'I2: ',C,dR_intgd(z2),
c     +      dR_intgd(z3),pow,dR
	    cdR_dOdC = cdR_dOdC + dR
	endif

	cdR_dOdC = pi4 * n0 * cdR_dOdC * c100 / h

	return
	end

c-----------------------------------------------------------------------
	function dR_intgd (tz)

c---	"dR INTeGranD"
c
c	The integrand for the z integral, wrt transformed z.

c+++  Arguments:
	real*8 dR_intgd, tz

c+++  Locals:
	real*8 c100, pi4, zero
	parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
	parameter (zero = 0.)
	real*8 z, Ndot, d, d2, z1
	real*8 F_spec, d_prop

c+++  Globals:
	real*8 gC, pow, rpow1
	common /ccb/ gC, pow, rpow1
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

	z = tz**rpow1
	z1 = 1. + z
	d = d_prop(z)
	d2 = d*d
	Ndot = pi4 * d2 * z1 * gC / F_spec(z)
	Ndot = Ndot / Nu
	dR_intgd = z**pow * d2*d2 * z1**(-beta) * A * Ndot**(-p) /
     +         ( z1 * sqrt(1.+Omega*z) * F_spec(z) * Nu )
c	write(9,'(4(1pg12.4))') gC, z, Ndot, dR_intgd

	return
	end

c-----------------------------------------------------------------------
	function d_prop (z)

c---	"Distance (PROPer)"

c+++  Arguments:
	real*8 d_prop, z

c+++  Locals:
	real*8 c100
	parameter (c100 = 2997.9246)

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

	d_prop = c100 * (z*q0 + (1.-q0)*(1.-sqrt(1+Omega*z))) /
     +           (h * q0**2 * (1.+z))

	return
	end

c-----------------------------------------------------------------------
	subroutine z_range (C, z1, z2, zm)

c---	"Z RANGE for integrations"
c
c	Return the range of z such that the luminosity function
c	at Ndot corresponding to count rate C is nonzero.
c
c	zm is where the upper spectral cutoff just enters the
c	trigger range.

c+++  Arguments:
	real*8 C, z1, z2, zm

c+++  Locals:
	real*8 c100, pi4, tol, ddel, zero, Mpc
	integer itmax
	parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
	parameter (tol = 0.0001, ddel = 0.001, zero = 0., itmax = 25)
	parameter (Mpc = 3.08568d24)
	real*8 x, xx, qq, xcur, z, dx, zmax, del
	integer it
	real*8 F_spec, d_prop

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

c---  Note the maximum z (pushes the upper spectral cutoff below
c---  the lower detector limit).
	zmax = eu/e1 - 1.

c---  First go for z1.  Get a rough guess...
	x = Nl / (pi4 * C)
	xx = x * h**2 / c100**2
	qq = (1.-q0)
	if (qq .gt. 0.) then
	    qq = 2. * qq * (1 + 0.25*qq)
	    z1 = sqrt( (sqrt(1+2.*xx*F_spec(zero)*qq)-1.)/qq )
	else
	    z1 = 0.5 * (sqrt(1.+4.*xx*F_spec(zero)) - 1.)
c	    print *, 'Special z1: ',xx,F_spec(zero),z1
	endif

c---  Now refine it with Newton-Raphson.
	it = 0
	del = min(0.003*z1,ddel)
20	xcur = (1.+z1) * d_prop(z1)**2 / F_spec(z1)
c	write(*,'(a,i3,4g12.4)') 'Guess1: ',it,z1,xcur,x
	if (abs(x-xcur) .gt. tol*x) then
	    it = it + 1
	    if (it .gt. itmax) pause 'z_range:  Too many NRs for z1!'
	    z = z1 + del
	    dx = ((1.+z) * d_prop(z)**2 / F_spec(z)  -  xcur) / del
	    z1 = z1 - (xcur-x)/dx
	    if (z1 .gt. zmax) z1 = sqrt((z-del)*zmax)
	    if (z1 .lt. 0.) z1 = sqrt(z-del)
	    goto 20
	endif

c---  Now do z2.  Get a rough guess...
	x = Nu / (pi4 * C)
	xx = x * h**2 / c100**2
	qq = (1.-q0)
	if (qq .gt. 0.) then
	    qq = 2. * qq * (1 + 0.25*qq)
	    z2 = sqrt( (sqrt(1+2.*xx*F_spec(zero)*qq)-1.)/qq )
	else
	    z2 = 0.5 * (sqrt(1.+4.*xx*F_spec(zero)) - 1.)
c	    print *, 'Special z2: ',xx,F_spec(zero),z2
	endif

c---  Now refine it with Newton-Raphson.
	it = 0
	del = min(0.003*z2,ddel)
40	xcur = (1.+z2) * d_prop(z2)**2 / F_spec(z2)
c	write(*,'(a,i3,4g12.4)') 'Guess2: ',it,z2,xcur,x
	if (abs(x-xcur) .gt. tol*x) then
	    it = it + 1
	    if (it .gt. itmax) pause 'z_range:  Too many NRs for z2!'
	    z = z2 + del
	    dx = ((1.+z) * d_prop(z)**2 / F_spec(z)  -  xcur) / del
	    z2 = z2 - (xcur-x)/dx
	    if (z2 .gt. zmax) z2 = sqrt((z-del)*zmax)
	    if (z2 .lt. 0.) z2 = sqrt(z-del)
	    goto 40
	endif
	if (zcut .gt. 0. .and. z2 .gt. zcut) z2 = zcut

c---  Last, find zm.
	zm = (eu/e2) - 1.

	return
	end

c-----------------------------------------------------------------------
	function F_spec (z)

c---	"Fraction of SPECtrum"

c+++  Arguments:
	real*8 F_spec, z

c+++  Locals:
	real*8 r, f, x1, x2

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

	r = el / eu
	f = (1. + z)/eu
	x1 = max( e1*f, r )
	x2 = min( e2*f, 1. )
	if (x2 .le. x1) then
	    F_spec = 0.
	else
	    if (alpha .eq. 1.) then
	        F_spec = K * log(x2/x1)
	    else
	        F_spec = K * ((x2/x1)**(alpha-1.) - 1.) * 
     *                   (1.+z)**(1.-alpha)
	    endif
	endif

	return
	end

c-----------------------------------------------------------------------
	function cdR_dz (Cmin, z)

c---	"Cosmological dR/dz"
c
c	Return the burst rate per unit redshift from redshift z, for
c	bursts with observed C above Cmin.

c+++  Arguments:
	real*8 cdR_dz, Cmin, z

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

c+++  Locals:
	real*8 c100, pi4
	parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
	real*8 area, z1, NCmin, N1r
	real*8 d_prop, F_spec

c---  If z is out of range, return 0.
	if (zcut .gt. 0. .and. z .gt. zcut) then
	    cdR_dz = 0.
	    return
	endif

c---  Find the Ndot range to integrate over.
	area = pi4 * d_prop(z)**2
	z1 = 1. + z
	NCmin = area * z1 * Cmin / F_spec(z)
	if (NCmin .gt. Nu) then
	    cdR_dz = 0.
	    return
	endif
	N1r = max(Nl, NCmin) / Nu

c---  The integral is analytic, since it's just a powerlaw.
	cdR_dz = A * (N1r**(1.-p) - 1) / (p - 1.)

c---  Multiply by the rate * area-redshift factor.
	cdR_dz = cdR_dz * n0 * area * c100 / 
     +           (z1**(beta+2.) * sqrt(1.+Omega*z) * h)

	return
	end

c-----------------------------------------------------------------------
	function scdR_dOdC (C, nu_sc)

c---	"Standard Candle cosmological dR/dCd(Omega)"
c
c	Return the burst rate per unit redshift from redshift z, for
c	bursts with observed C above Cmin.
c
c	This uses globals set by the "set" subroutines above,
c	except that set_lfun need not be called.

c+++  Arguments:
	real*8 scdR_dOdC, C, nu_sc

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

c+++  Locals:
	real*8 c100, pi4, tol, ddel, zero, Mpc
	integer itmax
	parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
	parameter (tol = 0.0001, ddel = 0.001, zero = 0., itmax = 25)
	parameter (Mpc = 3.08568d24)
	real*8 Ndot, x, xx, qq, xcur, z, z2, dx, zmax, del, z1
	integer it
	real*8 dF, r, dd, q1, rt, dp, dCdz
	real*8 F_spec, d_prop

c---  Convert from nu to Ndot.
	Ndot = nu_sc * pi4 * c100**2 / (h**2 * F_spec(zero))

c---  Note the maximum z: the z that pushes the upper 
c---  spectral cutoff below the lower detector limit.
	zmax = eu/e1 - 1.

c---  First find the z of the bursts.  Start with a rough guess...
	x = Ndot / (pi4 * C)
	xx = x * h**2 / c100**2
	qq = (1.-q0)
	if (qq .gt. 0.) then
	    qq = 2. * qq * (1 + 0.25*qq)
	    z = sqrt( (sqrt(1+2.*xx*F_spec(zero)*qq)-1.)/qq )
	else
	    z = 0.5 * (sqrt(1.+4.*xx*F_spec(zero)) - 1.)
	endif
	if (F_spec(z) .eq. 0) pause 'SC: F_spec vanishes!'

c---  Now refine it with Newton-Raphson.
	it = 0
	del = min(0.003*z,ddel)
20	xcur = (1.+z) * d_prop(z)**2 / F_spec(z)
c	write(*,'(a,i3,6g12.4)') 'Guess: ',it,z,xcur,x,(xcur-x)/dx,zmax
	if (abs(x-xcur) .gt. tol*x) then
	    it = it + 1
	    if (it .gt. itmax) pause 'z_range:  Too many NRs for z!'
	    z2 = z + del
	    dx = ((1.+z2) * d_prop(z2)**2 / F_spec(z2)  -  xcur) / del
	    z = z - (xcur-x)/dx
	    if (z .gt. zmax) z = sqrt((z2-del)*zmax)
	    if (z .lt. 0.) z = 0.5 * (z2-del)
c	    if (z .lt. 0.) z = sqrt(z2-del)
	    goto 20
	endif
	z1 = 1. + z

c---  If beyond zcut, return 0.
	if (zcut .gt. 0. .and. z .gt. zcut) then
	    scdR_dOdC = 0.
	    return
	endif

c---  Calculate F0'/F0.
	if (z1 .gt. eu/e2) then
	    r = eu / (e1*z1)
	    if (alpha .ne. 1.) then
	        dF = (1.-alpha) * r**(alpha-2.) /
     +               (z1 * (r**(alpha-1.) - 1.))
	    else
	        dF = - 1. / (z1 * log(r))
	    endif
	else
	    if (alpha .ne. 1.) then
	        dF = (1. - alpha) / z1
	    else
	        dF = 0.
	    endif
	endif

c---  Calculate d'/d.
	dp = d_prop(z)
	q1 = 1. - q0
	rt = sqrt(Omega*z + 1.)
	dd = c100 / (h * q0**2 * z1)
	dd = dd * (q0 - q0*q1/rt - (z*q0 + q1*(1.-rt))/z1)
	dd = dd / d_prop(z)

c---  Use these to calculate the Jacobian to convert the delta-function
c---  from C to z.
	dCdz = Ndot * F_spec(z) / (pi4 * dp**2 * z1)
	dCdz = dCdz * abs(dF - 2.*dd - 1./z1)

c---  Finally, here's the rate.
	scdR_dOdC = (c100/h) * dp**2 * n0 /
     +           (dCdz * z1**(beta+2.) * sqrt(1.+Omega*z))

	return
	end

c-----------------------------------------------------------------------
	function scdR_dz (nu_sc, Cmin, z)

c---	"Standard Candle cosmological dR/dz"
c
c	Return the burst rate per unit redshift from redshift z, for
c	bursts with observed C above Cmin.

c+++  Arguments:
	real*8 scdR_dz, nu_sc, Cmin, z

c+++  Globals:
	real*8 h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K
	real*8 n0, beta, e1, e2, zcut
	common /coscb/ h, Omega, q0, p, Nl, Nu, A, alpha, el, eu, K, 
     +          n0, beta, e1, e2, zcut

c+++  Locals:
	real*8 c100, pi4, zero
	parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
	parameter (zero = 0.)
	real*8 Ndot, area, z1, Cz
	real*8 d_prop, F_spec

c---  If z=0, there is no volume!  Return to avoid division by 0.
	if (z .eq. 0.) then
	    scdR_dz = 0.
	    return
	endif

c---  If beyond zcut, return 0.
	if (zcut .gt. 0. .and. z .gt. zcut) then
	    scdR_dz = 0.
	    return
	endif

c---  Check that there are bursts available at the required C values.
	Ndot = nu_sc * pi4 * c100**2 / (h**2 * F_spec(zero))
	area = pi4 * d_prop(z)**2
	z1 = 1. + z
	Cz = Ndot * F_spec(z) / (area * z1)
	if (Cz .lt. Cmin) then
	    scdR_dz = 0.
	    return
	endif

c---  Here it is...
	scdR_dz = n0 * area * c100 / 
     +           (z1**(beta+2.) * sqrt(1.+Omega*z) * h)

	return
	end

c-----------------------------------------------------------------------
	

c-----------------------------------------------------------------------



