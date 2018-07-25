c=======================================================================
c	Subroutines for modeling the spatial distribution of GRBs.
c
c	30 Oct 91 TJL
c	   Feb 92 rewritten
c	09 Mar 92 changed variables in integrals to x=E^-alpha
c	10 Apr 92 further mods to speed integration; they assume
c		  powerlaw E dist'n
c=======================================================================
c	This first bunch sets the models for spatial & energy dist'ns.
c=======================================================================
	subroutine set_R0 (R_0)

c	"SET R0"
c
c	Set radial distance to Sol along GC x-axis.

c+++  Arguments:
	real*8 R_0

c+++  Globals:
	real*8 R0
	common /Rcb/ R0
	data R0/0./

	R0 = R_0

	return
	end

c-----------------------------------------------------------------------
	subroutine set_n_mod (model, np)

c---	"SET Number density MODel"
c
c	Specify what burst rate per unit volume model to use.
c	Returns the number of parameters required in np.

c+++  Arguments:
	integer model, np

c+++  Globals:
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel
	data ndpar/4,5/, dmodel/0/

	if (dmodel .gt. dmmax) pause 'set_n_mod:  No such model!'
	dmodel = model
	np = ndpar(model)

	return
	end

c-----------------------------------------------------------------------
	subroutine set_n_par (params)

c---	"SET Number density PARameters"
c
c	Specify parameters for the burst number density model.

c+++  Arguments:
	real*8 params(*)

c+++  Locals:
	integer i

c+++  Globals:
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel

c---  Power-law, bounded, spherical model (C,p,Rl,Ru).  n=Cr^-p
	if (dmodel .eq. 1) then
	    do 20 i=1, ndpar(1)
	        dpar(i) = params(i)
20	    continue

c---  Oblate spheroid model, 1/r**2 with core (C,e,rc,Rl,Ru).
c---  Store (minor axis)**2 rather than e, and rc**2 rather than rc.
	else if (dmodel .eq. 2) then
	    do 40 i=1, ndpar(2)
	        dpar(i) = params(i)
40	    continue
	    if (dpar(4) .eq. 0. .and. dpar(3) .eq. 0.) 
     +         pause 'Cannot have zero core and no inner boundary!'
	    if (dpar(2) .ge. 1. .or. dpar(2) .lt. 0.)
     +         pause 'Illegal eccentricity!'
	    dpar(2) = 1. - dpar(2)**2
	    dpar(3) = dpar(3)**2
	    
	else 
	    pause 'set_n:  No such model!'
	endif

	return
	end

c-----------------------------------------------------------------------
	subroutine set_E_mod(model, np)

c---	"SET Energy distribution"
c
c	Specify the burst energy distribution model.
c	Return # of parameters required in np.

c+++  Arguments:
	integer model, np

c+++  Globals:
	integer emmax, epmax
	parameter (emmax = 3, epmax = 5)
	integer emodel, nepar(emmax)
	real*8 epar(epmax)
	common /emodcb/ epar, nepar, emodel
	data nepar/3, 3, 3/, emodel/0/

	if (model .gt. emmax) pause 'set_E_mod:  No such model!'
	emodel = model
	np = nepar(model)

	return
	end

c-----------------------------------------------------------------------
	subroutine set_E_par (params)

c---	"SET Energy dist'n PARameters"
c
c	Specify parameters for the burst energy dist'n.

c+++  Arguments:
	real*8 params(*)

c+++  Locals:
	integer i
	real*8 pi4
	parameter (pi4 = 4.d0*3.141592653589793d0)

c+++  Globals:
	integer emmax, epmax
	parameter (emmax = 3, epmax = 5)
	integer emodel, nepar(emmax)
	real*8 epar(epmax)
	common /emodcb/ epar, nepar, emodel
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel
	real*8 R0
	common /Rcb/ R0

c---  Bounded power-law model (p,El,Eu,[A]).  A is calculated here.
	if (emodel .eq. 1) then
	    do 20 i=1, nepar(1)
	        epar(i) = params(i)
20	    continue
	    if (epar(1) .eq. 1.) then
	        epar(4) = 1. / log(epar(3)/epar(2))
	    else
	        epar(4) = (epar(1)-1.) / 
     +                    ((epar(3)/epar(2))**(epar(1)-1.) - 1.)
	        epar(4) = epar(4) * epar(3)**epar(1)
	    endif

c---  Bounded power-law, with different parameters (p, r, Sc).
c---  r = Eu/El, Sc = Eu/4pi(R0**2+Rc**2).  This is only valid
c---  if the current density model is model 2.  Store
c---  El and Eu rather than r and Sc.  Also, calculate & save A.
	else if (emodel .eq. 2) then
	    if (dmodel .ne. 2) pause 'E model 2 requires halo model!'
	    epar(1) = params(1)
	    epar(3) = params(3) * pi4 * (r0**2 + dpar(3))
	    epar(2) = epar(3) / params(2)
	    if (epar(1) .eq. 1.) then
	        epar(4) = 1. / log(epar(3)/epar(2))
	    else
	        epar(4) = (epar(1)-1.) / 
     +                    ((epar(3)/epar(2))**(epar(1)-1.) - 1.)
	        epar(4) = epar(4) * epar(3)**epar(1)
	    endif

c---  Bounded power-law, with different parameters (p, r, r1).
c---  r = Eu/El, Eu = 4pi*r1**2 (assumes Smin=1).  Store
c---  El and Eu rather than r and rmax.  Also, calculate & save A.
	else if (emodel .eq. 3) then
	    epar(1) = params(1)
	    epar(3) = pi4 * params(3)**2
	    epar(2) = epar(3) / params(2)
	    write(*,'(a,2g13.4)') 'E range: ',epar(2),epar(3)
	    if (epar(1) .eq. 1.) then
	        epar(4) = 1. / log(epar(3)/epar(2))
	    else
	        epar(4) = (epar(1)-1.) / 
     +                    ((epar(3)/epar(2))**(epar(1)-1.) - 1.)
	        epar(4) = epar(4) * epar(3)**epar(1)
	    endif

	else 
	    pause 'set_E_par:  No such model!'
	endif

	return
	end

c-----------------------------------------------------------------------
c=======================================================================
c	Here are the models; spatial coords are wrt the GC, with
c	z towards NGP and x towards Sol.
c=======================================================================
	function dn_dV (x, y, z)

c---	"dn / dV"
c
c	Return the burst rate per unit volume as a function of
c	GC coordinates.

c+++  Arguments:
	real*8 dn_dV, x, y, z

c+++  Locals:
	real*8 r
	real*8 del
	parameter (del = 1. + 1.d-10)

c+++  Globals:
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel

c---  Bounded spherical distribution.
	if (dmodel .eq. 1) then
	    r = sqrt(x*x + y*y + z*z)
	    if (r .ge. dpar(3)/del .and. r .le. dpar(4)*del) then
	        dn_dV = dpar(1)
	        if (dpar(2) .ne. 0.) dn_dV = dn_dV / r**dpar(2)
	    else
	        dn_dV = 0.
	    endif

c---  Bounded 1/r**2 oblate spheroid distribution with core.
	else if (dmodel .eq. 2) then
	    r = sqrt(x*x + y*y + z*z/dpar(2))
	    if (r .ge. dpar(4)/del .and. r .le. dpar(5)*del) then
	        dn_dV = dpar(1) / (dpar(3) + r**2)
	    else
	        write(*,'(a,4(1pg13.4))') 'Outside -- r= ',r,x,y,z
	        dn_dV = 0.
	    endif

	endif

	return
	end

c-----------------------------------------------------------------------
	function E_dist (E)

c---	"Energy DISTribution"
c
c	Return the probability density for a burst having an energy E.

c+++  Arguments:
	real*8 E_dist, E

c+++  Locals:
	real*8 del
	parameter (del = 1. + 1.d-10)

c+++  Globals:
	integer emmax, epmax
	parameter (emmax = 3, epmax = 5)
	integer emodel, nepar(emmax)
	real*8 epar(epmax)
	common /emodcb/ epar, nepar, emodel

c---  Power-law distribution:
	if (emodel .le. 3) then
	    if (E .ge. epar(2)/del .and. E .le. epar(3)*del) then
	        E_dist = epar(4) / E**epar(1)
	    else
	        E_dist = 0.
	        write(*,'(a,1pg13.4)') 'Outside -- E= ',E
	    endif
	endif

	return
	end

c=======================================================================
c	This bunch returns ranges of the coordinates & energy.
c	Here, spatial coordinate origin is at Sol, with z towards GC
c	and x towards NGP.
c=======================================================================
	subroutine r_range_1 (mu, phi, rl, ru)

c	"Radial RANGE 1"
c
c	Range of r inside the spatial dist'n in the specified drxn,
c	for the 1st intersection of the l.o.s. with the distribution.

c+++  Arguments:
	real*8 mu, phi, rl, ru

c+++  Locals:
	real*8 one
	parameter (one = 1.d0)
	real*8 muc, outer, inner
	real*8 mu2, C, din, dout, rm, rm2, disc1, disc2

c+++  Globals:
	real*8 R0
	common /Rcb/ R0
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel

c---  Bounded spherical power-law model:
	if (dmodel .eq. 1) then
	    inner = dpar(3)
	    outer = dpar(4)
c>>>  R0 = 0 is a simple special case.
	    if (R0 .eq. 0.) then
	        rl = inner
	        ru = outer
c>>>  If inside the inner boundary, there is one simple intsxn.
	    else if (R0 .le. inner) then
	        rl = R0 * (sqrt(mu*mu-one+(inner/R0)**2) + mu)
	        ru = R0 * (sqrt(mu*mu-one+(outer/R0)**2) + mu)
c>>>  If between the boundaries, rl=0, and ru depends on whether
c>>>  the los hits the inner boundary or not.
	    else if (R0 .gt. inner .and. R0 .lt. outer) then
	        rl = 0.
	        muc = sqrt(one - (inner/R0)**2)
	        if (mu .gt. muc) then
	            ru = R0 * (mu - sqrt(mu*mu-one+(inner/R0)**2))
	        else
	            ru = R0 * (sqrt(mu*mu-one+(outer/R0)**2) + mu)
	        endif
c>>>  Finally, if outside the outer boundary, there is no intsxn if
c>>>  mu is too small.  Otherwise, rl is set by the first intsxn with
c>>>  the outer boundary, and ru depends on whether the los hits
c>>>  the inner boundary.
	    else
	        muc = sqrt(one - (outer/R0)**2)
	        if (mu .gt. muc) then
	            rl = R0 * (mu - sqrt(mu*mu-one+(outer/R0)**2))
	            muc = sqrt(one - (inner/R0)**2)
	            if (mu .gt. muc) then
	                ru = R0 * (mu - sqrt(mu*mu-one+(inner/R0)**2))
	            else
	                ru = R0 * (mu + sqrt(mu*mu-one+(outer/R0)**2))
	            endif
	        else
	            rl = 0.
	            ru = 0.
	        endif
	    endif

c---  Bounded oblate spheroid.
	else if (dmodel .eq. 2) then
	    inner = dpar(4)
	    outer = dpar(5)
	    mu2 = mu*mu
	    C = sqrt(mu2 + (1.-mu2)*(sin(phi)**2+cos(phi)**2/dpar(2)))
	    rm = R0*mu
	    rm2 = rm*rm
	    din = R0*R0 - inner*inner
	    dout = R0*R0 - outer*outer
c>>>  R0 = 0 is a simpler special case.
	    if (R0 .eq. 0.) then
	        rl = inner * C
	        ru = outer * C
c>>>  If inside the inner boundary, there is one simple intsxn.
	    else if (R0 .le. inner) then
	        rl = (rm + sqrt(rm2 - din*C))/C
	        ru = (rm + sqrt(rm2 - dout*C))/C
c>>>  If between the boundaries, rl=0, and ru depends on whether
c>>>  the los hits the inner boundary or not, which is determined
c>>>  by the sign of the discriminant (rm2 - din*C).
	    else if (R0 .gt. inner .and. R0 .lt. outer) then
	        rl = 0.
	        disc1 = rm2 - din*C
	        if (disc1 .gt. 0.) then
	            ru = (rm - sqrt(disc1))/C
	        else
	            ru = (rm + sqrt(rm2 - dout*C))/C
	        endif
c>>>  Finally, if outside the outer boundary, there is no intsxn if
c>>>  mu is too small.  Otherwise, rl is set by the first intsxn with
c>>>  the outer boundary, and ru depends on whether the los hits
c>>>  the inner boundary.  The 1st cond'n depends on the sign of
c>>>  (rm2 - dout*C), the second on (rm2 - din*C).
	    else
	        disc1 = rm2 - dout*C
	        if (disc1 .gt. 0.) then
	            rl = (rm - sqrt(disc1))/C
	            disc2 = rm2 - din*C
	            if (disc2 .gt. 0.) then
	                ru = (rm - sqrt(disc2))/C
	            else
	                ru = (rm + sqrt(disc1))/C
	            endif
	        else
	            rl = 0.
	            ru = 0.
	        endif
	    endif

	endif

c	write(9,'(a,4(1pg13.4))') 'r1: ',mu,phi,rl,ru
c	call r_range(mu,phi)
	return
	end

c-----------------------------------------------------------------------
	subroutine r_range_2 (mu, phi, rl, ru)

c	"Radial RANGE 2"
c
c	Range of r inside the spatial dist'n in the specified drxn,
c	for a possible 2nd intersection of the l.o.s. with the 
c	distribution.

c+++  Arguments:
	real*8 mu, phi, rl, ru

c+++  Locals:
	real*8 one
	parameter (one = 1.d0)
	real*8 muc, inner, outer
	real*8 mu2, C, din, dout, rm, rm2, disc1, disc2

c+++  Globals:
	real*8 R0
	common /Rcb/ R0
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel

c---  Bounded spherical power-law model:
	if (dmodel .eq. 1) then
	    inner = dpar(3)
	    outer = dpar(4)
c>>>  If R0=0, or if inside the inner boundary, there was only 1 intsxn.
	    if (R0 .eq. 0. .or. R0 .le. inner) then
	        rl = 0.
	        ru = 0.
c>>>  If between the boundaries, there is a 2nd intsxn only if
c>>>  the los hits the inner boundary.
	    else if (R0 .gt. inner .and. R0 .lt. outer) then
	        muc = sqrt(one - (inner/R0)**2)
	        if (mu .gt. muc) then
	            rl = R0 * (mu + sqrt(mu*mu-one+(inner/R0)**2))
	            ru = R0 * (mu + sqrt(mu*mu-one+(outer/R0)**2))
	        else
	            rl = 0.
	            ru = 0.
	        endif
c>>>  If outside the outer boundary, there is a 2nd intsxn only if
c>>>  the los hits both the outer and inner boundaries.
	    else
	        muc = sqrt(one - (outer/R0)**2)
	        if (mu .gt. muc) then
	            muc = sqrt(one - (inner/R0)**2)
	            if (mu .gt. muc) then
	                rl = R0 * (mu + sqrt(mu*mu-one+(inner/R0)**2))
	                ru = R0 * (mu + sqrt(mu*mu-one+(outer/R0)**2))
	            else
	                rl = 0.
	                ru = 0.
	            endif
	        else
	            rl = 0.
	            ru = 0.
	        endif
	    endif

c---  Bounded oblate spheroid.
	else if (dmodel .eq. 2) then
c>>>  If R0=0, or if inside the inner boundary, there was only 1 intsxn.
	    if (R0 .eq. 0. .or. R0 .le. inner) then
	        rl = 0.
	        ru = 0.
	        return
	    endif
	    inner = dpar(4)
	    outer = dpar(5)
	    mu2 = mu*mu
	    C = sqrt(mu2 + (1.-mu2)*(sin(phi)**2+cos(phi)**2/dpar(2)))
	    rm = R0*mu
	    rm2 = rm*rm
	    din = R0*R0 - inner*inner
	    dout = R0*R0 - outer*outer
c>>>  If between the boundaries, there is a 2nd intsxn only if
c>>>  the los hits the inner boundary.
	    if (R0 .gt. inner .and. R0 .lt. outer) then
	        disc1 = rm2 - din*C
	        if (disc1 .gt. 0.) then
	            rl = (rm + sqrt(disc1))/C
	            ru = (rm + sqrt(rm2 - dout*C))/C
	        else
	            rl = 0.
	            ru = 0.
	        endif
c>>>  If outside the outer boundary, there is a 2nd intsxn only if
c>>>  the los hits both the outer and inner boundaries.
	    else
	        disc1 = rm2 - dout*C
	        disc2 = rm2 - din*C
	        if (disc1 .gt. 0. .and. disc2 .gt. 0.) then
	            rl = (rm + sqrt(disc2))/C
	            ru = (rm + sqrt(disc1))/C
	        else
	            rl = 0.
	            ru = 0.
	        endif
	    endif

	endif

c	write(9,'(a,4(1pg13.4))') 'r2: ',mu,phi,rl,ru
	return
	end

c-----------------------------------------------------------------------
	subroutine E_range (el, eh)

c---	"Energy RANGE"
c
c	This returns the range of energies covered by the current
c	energy distribution.

c+++  Arguments:
	real*8 el, eh

c+++  Globals:
	integer emmax, epmax
	parameter (emmax = 3, epmax = 5)
	integer emodel, nepar(emmax)
	real*8 epar(epmax)
	common /emodcb/ epar, nepar, emodel

	if (emodel .eq. 1 .or. emodel .eq. 2 .or. emodel .eq. 3) then
	    el = epar(2)
	    eh = epar(3)
	endif

	return
	end

c-----------------------------------------------------------------------
	subroutine r_min (rmin)

c---	"Radial coordinate Minimum"
c
c	This returns the minimum radius from Sol that lies 
c       inside the dist'n.

c+++  Arguments:
	real*8 rmin

c+++  Locals:
	real*8 zero
	parameter (zero = 0.)

c+++  Globals:
	real*8 R0
	common /Rcb/ R0
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel

c---  Bounded sphere:
	if (dmodel .eq. 1) then
	    if (R0 .lt. dpar(3)) then
	        rmin = dpar(3) - R0
	    else if (R0 .ge. dpar(3) .and. R0 .le. dpar(4)) then
	        rmin = 0.
	    else
	        rmin = R0 - dpar(4)
	    endif

c---  Bounded oblate spheroid:
c---  If we are close enough to the origin, the nearest point
c---  moves away from (z=0) towards (x=0).
	else if (dmodel .eq. 2) then
	    if (R0 .lt. dpar(4)) then
	        if (R0 .ge. dpar(4)*(1.-dpar(2))) then
	            rmin = dpar(4) - R0
	        else
	            rmin = sqrt(dpar(2)*dpar(4)**2 - 
     +                          R0**2*dpar(2)/(1.-dpar(2)))
	        endif
	    else if (R0 .ge. dpar(4) .and. R0 .le. dpar(5)) then
	        rmin = 0.
	    else
	        if (R0 .ge. dpar(5)*(1.-dpar(2))) then
	             rmin = R0 - dpar(5)
	        else
	            rmin = sqrt(dpar(2)*dpar(5)**2 - 
     +                          R0**2*dpar(2)/(1.-dpar(2)))
	        endif
	    endif

	endif

	return
	end

c-----------------------------------------------------------------------
	subroutine mu_range (mulo, muhi, f)

c	"MU RANGE"
c
c	Range of polar angle cosine (mu = 1 toward GC) that lies
c	inside the dist'n.  f is a factor by which integrals over
c	mu should be multiplied, to account for possible symmetry.

c+++  Arguments:
	real*8 mulo, muhi, f

c+++  Locals:
	real*8 one
	parameter (one = 1.d0)

c+++  Globals:
	real*8 R0
	common /Rcb/ R0
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel

c---  Bounded spherical power-law model:
c---  The range is -1 to 1, unless R0 is outside the dist'n.
	if (dmodel .eq. 1) then
	    if (R0 .ge. dpar(4)) then
	        mulo = sqrt(one - (dpar(4)/R0)**2)
	    else
	        mulo = - one
	    endif
	    muhi = one
	    f = one

c---  Bounded oblate spheroid:
c---  The range is -1 to 1, unless R0 is outside the dist'n, in
c---  which case mulo is set by the widest part (phi = pi/2).
c---  This gives the same value of mulo as for the spherical case.
	else if (dmodel .eq. 2) then
	    if (R0 .ge. dpar(5)) then
	        mulo = sqrt(one - (dpar(5)/R0)**2)
	    else
	        mulo = - one
	    endif
	    muhi = one
	    f = one

	endif

	return
	end

c-----------------------------------------------------------------------
	subroutine phi_range (philo, phihi, f)

c	"PHI RANGE"
c
c	Range of azimuth to integrate.  f is a factor by which 
c	integrals over mu should be multiplied, to account for 
c	possible symmetry.


c+++  Arguments:
	real*8 philo, phihi, f

c+++  Locals:
	real*8 pi
	parameter (pi = 3.141592653589793d0)

c+++  Globals:
	real*8 R0
	common /Rcb/ R0
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel

c---  Bounded spherical power-law model:
c---  Integrate over 0 to pi/2, and multiply by 4 for symmetry.
c---  (Actually, any single phi point is enough!)
	if (dmodel .eq. 1) then
	    philo = 0.
	    phihi = pi / 2.
	    f = 4.

c---  Bounded oblate spheroid.
c---  Integrate over 0 to pi/2, and multiply by 4 for symmetry.
	else if (dmodel .eq. 2) then
	    philo = 0.
	    phihi = pi / 2.
	    f = 4.

	endif

	return
	end

c=======================================================================
c	This bunch calculates the observed distributions.
c=======================================================================
	function dn_dSdO (S, mu, phi, eps)

c---	"dn / d S d Omega"
c
c	Return the burst rate per unit size per unit steradian 
c	in the specified direction at the specified size.
c
c	eps sets the accuracy of the Romberg quadratures used.
c
c	Note on power-law transformations:
c
c	We transform the dependent variables of integrals to new
c	params related by a powerlaw to ease computing of the
c	integrals.  This integral behaves roughly as follows:
c	If f(E) ~ E^-p
c	and n(r) ~ r^-s, then alpha = p + s/2 - 3/2, and
c	the change of variables for E is to x = E^(1-alpha),
c	with dE =  E^+alpha dx/(1-alpha).  We store 1-alpha.

c+++  Arguments:
	real*8 dn_dSdO, S, mu, phi, eps

c+++  Locals:
	real*8 pi32, pi4, one
	parameter (pi32 = 5.56832800d0, pi4 = 4.d0*3.14159265358d0)
	parameter (one = 1.d0)
	real*8 rl, ru, El, Eu, E1, E2, dn2, x1, x2, ratio, p
	real*8 i_dSdO
	external i_dSdO

c+++  Globals:
	real*8 gS, gmu, gphi
	common /dsdocb/ gS, gmu, gphi
	real*8 alpha1, ra1
	common /transcb/ alpha1, ra1

	gS = S
	gmu = mu
	gphi = phi

c--1  The integral may have two parts.  Here's the first part:
	call e_range (El, Eu)
	call r_range_1 (mu, phi, rl, ru)
	if (ru .gt. rl) then
	    E1 = max(El, pi4*S*rl*rl)
	    E2 = min(Eu, pi4*S*ru*ru)
	    if (E2 .gt. E1) then

c--1  Determine the power law to use for transforming the integral.
	        alpha1 = 1.
	        ra1 = 1.
	        ratio = i_dSdO(E2) / i_dSdO(E1)
	        p = - log(ratio) / log(E2/E1)
	        if (p .lt. one) p = p + .5
	        if (p .eq. one) then
	            alpha1 = one
	        else
	            alpha1 = one - p
	        endif
	        ra1 = one / alpha1
	        x1 = E1**alpha1
	        x2 = E2**alpha1
	        call romb(i_dSdO, x1, x2, dn_dSdO, eps)
	    else
	        dn_dSdO = 0.
	    endif
c	    write(9,'(12x,5(1pg12.4))') E1,E2,El,Eu,dn_dSdO

c--2  Here's the possible second part:
	    call r_range_2 (mu, phi, rl, ru)
	    if (ru .gt. rl) then
	        E1 = max(El, pi4*S*rl*rl)
	        E2 = min(Eu, pi4*S*ru*ru)
	        if (E2 .gt. E1) then

c--2  Determine the power law to use for transforming the integral.
	            alpha1 = 1.
	            ra1 = 1.
	            ratio = i_dSdO(E2) / i_dSdO(E1)
	            p = - log(ratio) / log(E2/E1)
	            if (p .lt. one) p = p + .5
	            if (p .eq. one) then
	                alpha1 = one
	            else
	                alpha1 = one - p
	            endif
	            ra1 = one / alpha1
	            x1 = E1**alpha1
	            x2 = E2**alpha1
	            call romb(i_dSdO, x1, x2, dn2, eps)
	            dn_dSdO = dn_dSdO + dn2
	        endif
	    endif
	    dn_dSdO = dn_dSdO / (S**2.5d0 * 16.d0 * pi32)

c--0  If there isn't even a first part, it's really simple!
	else
	    dn_dSdO = 0.
	endif

	return
	end

c-----------------------------------------------------------------------
	function i_dSdO (xx)

c---	"Integrand for dn_dSdO"
c
c	Evaluate the integrand for dn_dSdO, for integration
c	over burst energy (transformed to xx=E^(1-alpha)).

c+++  Arguments:
	real*8 i_dSdO, xx

c+++  Locals:
	real*8 one, pi4
	parameter (one = 1.d0, pi4 = 4.d0*3.14159265358d0)
	real*8 E, r, s, x, y, z, pow
	real*8 dn_dV, E_dist

c+++  Globals:
	real*8 gS, gmu, gphi
	common /dsdocb/ gS, gmu, gphi
	real*8 R0
	common /Rcb/ R0
	real*8 alpha1, ra1
	common /transcb/ alpha1, ra1

c---  Transform from x to E.
	E = xx**ra1

c---  Calculate GC coordinates for the relevant point.
	r = sqrt(E/(pi4*gS))
	s = sqrt(one - gmu*gmu)
	x = R0 - r*gmu
	y = r * s * sin(gphi)
	z = r * s * cos(gphi)
	pow = 2.5 - alpha1

c---  Here's the integrand.
c	i_dSdO = E_dist(E) * E * sqrt(E) * dn_dV(x,y,z) <--- wrt E
	i_dSdO = E_dist(E) * E**pow * dn_dV(x,y,z) / alpha1
	if (i_dSdO .eq. 0.) write(*,'(7g11.3)')
     +  gS,gmu,gphi,E,r,sqrt(r*r+R0**2-2.*R0*r*gmu)

	return
	end

c-----------------------------------------------------------------------
	function dn_dO (Sl, mu, phi, eps)

c---	"dn / d Omega"
c
c	Return the burst rate per unit steradian in the specified
c	direction, for bursts with size > Sl.
c
c	eps sets the accuracy of the Romberg quadratures used.

c+++  Arguments:
	real*8 dn_dO, Sl, mu, phi, eps

c+++  Locals:
	real*8 pi4, one
	parameter (pi4 = 4.d0*3.14159265358d0, one = 1.d0)
        real*8 rl, ru, S_min, S_max, S1, El, Eu, E1, E2, Ea
        real*8 p1, i1, i2
        real*8 ratio, p
        real*8 inner, oromb, E_dist

	real*8 zero
        PARAMETER(zero=0.)

c+++  Globals:
	real*8 R0
	common /Rcb/ R0
	real*8 gSl, gmu, gphi, geps
	common /docb/ gSl, gmu, gphi, geps

c---  Make sure Sl is not larger than the largest possible S.
	call r_range_1 (mu, phi, rl, ru)
	call E_range (El, Eu)
	if (rl .gt. 0.) then
	    S_max = Eu / (pi4*rl*rl)
	    if (Sl .ge. S_max) then
	        dn_dO = 0.
	        return
	    endif
	endif

c---  Set global variables for communicating with innermost integral.
	gSl = Sl
	gmu = mu
	gphi = phi
	geps = eps

c---  Find range for outermost integration over E for the 1st r range.
	S_min = El / (pi4*ru*ru)
	S1 = max( Sl, S_min)
	E1 = max(El, pi4*S1*rl**2)
	E2 = Eu
	if (E1 .ge. E2) then
	    dn_dO = 0.
	    return
	endif

c---  Determine the power law to use for transforming the integral,
c---  and then do it.
	if (rl .gt. R0/5.) then
	    p = - log(E_dist(E2)/E_dist(E1)) / log(E2/E1)
	else
	    Ea = E1
	    i1 = inner(Ea,rl,ru,S1)
	    if (i1 .eq. 0.) then
	       Ea = min(2.*E1,sqrt(E1*E2))
	       i1 = inner(Ea,rl,ru,S1)
	    endif
	    i2 = inner(E2,rl,ru,S1)
	    if (i1*i2 .eq. 0.) then
	        p = 0.
	    else
	        ratio = i2 / i1
	        p = - log(ratio) / log(E2/Ea)
	    endif
	endif
	if (p .eq. one) then
	    p1 = one
	else
	    p1 = one - p
	endif
c	write(9,'(a,4g12.4)') 'r1: ',rl,ru,E1,E2
c	write(9,'(4x,4g12.4)') Ea,E2,i1,i2
	dn_dO = oromb(p1, E1, E2, rl, ru, S1)

c---  Now repeat for the 2nd r range.
	call r_range_2 (mu, phi, rl, ru)
	if (ru .le. rl) return

c---  First, find the E range to integrate over.
	S_min = El / (pi4*ru*ru)
	S1 = max( Sl, S_min)
	E1 = max(El, pi4*S1*rl**2)
	if (E1 .ge. E2) return

c---  Determine the power law to use for transforming the integral,
c---  and then do it.
	if (rl .gt. R0/5.) then
	    p = - log(E_dist(E2)/E_dist(E1)) / log(E2/E1)
	else
	    Ea = E1
	    i1 = inner(Ea,rl,ru,S1)
	    if (i1 .eq. 0.) then
	       Ea = min(2.*E1,sqrt(E1*E2))
	       i1 = inner(Ea,rl,ru,S1)
	    endif
	    if (i1*i2 .eq. 0.) then
	    p = 0.
	    else
	        ratio = i2 / i1
	        p = - log(ratio) / log(E2/Ea)
	    endif
	endif
	if (p .eq. one) then
	    p1 = one
	else
	    p1 = one - p
	endif
c	write(9,'(a,4g12.4)') 'r2: ',rl,ru,E1,E2
c	write(9,'(4x,4g12.4)') Ea,E2,i1,i2
	dn_dO = dn_dO + oromb(p1, E1, E2, rl, ru, S1)

	return
	end

c-----------------------------------------------------------------------
	function oromb (p1, E1, E2, rl, ru, S1)

c---	"Outer ROMBerg"
c
c	Perform the outermost integral over E for dn/dO.

c+++  Arguments:
	real*8 oromb, p1, E1, E2, rl, ru, S1

c+++  Locals:
	real*8 pi4, one
	parameter (pi4 = 4.d0*3.14159265358d0, one = 1.d0)
        real*8 E, y1, y2, py, y_rng, dy, y
        real*8 inner

	real*8 zero
	integer JMAX, JMAXP, K, KM
        PARAMETER(zero=0.,JMAX=20,JMAXP=JMAX+1,K=5,KM=4)
        integer j, l, i, it
        real*8 s(JMAXP), h(JMAXP), sum, tnm, del

c+++  Globals:
	real*8 gSl, gmu, gphi, geps
	common /docb/ gSl, gmu, gphi, geps

	y1 = E1**p1
	y2 = E2**p1
	py = one/p1
c	write(9,'(a,6g10.3)') 'mu,p1,E,i: ',gmu,p1,E1,E2

c---  Do the outer integral by Romberg with the trap. rule.  Do it
c---  explicitly here, so that we can call the qromb subroutine
c---  for the inner integral, and not have any "crosstalk".
c	call qromb(inner, y1, y2, oromb, eps)
	y_rng = y2 - y1
	h(1) = 1.
	do 40 j=1, jmax
	    if (j .eq. 1) then
	        s(j) = E1**(one-p1) * inner(E1,rl,ru,S1)
c	if (gmu.gt..23) write(9,'(a,i3,2g12.4)') '--> ',j,E,s(j)
	        s(j) = s(j) + E2**(one-p1) * inner(E2,rl,ru,S1)
c	if (gmu.gt..23) write(9,'(a,i3,2g12.4)') '--> ',j,E,s(j)
	        s(j) = 0.5 * y_rng * s(j) / p1
	        it = 1
	    else
	        tnm = it
	        dy = y_rng / tnm
	        y = y1 + 0.5*dy
	        sum = 0.
	        do 20 i=1, it
	            E = y**py
	            sum = sum + E**(one-p1) * inner(E,rl,ru,S1)
c	if (gmu.gt..23) write(9,'(a,i3,5g12.4)') 
c     .      '--> ',j,E,rl,min( ru, sqrt(E/(pi4*S1)) ),
c     .      inner(E,rl,ru,S1),E**(one-p1) *inner(E,rl,ru,S1)
	            y = y + dy
20	        continue
	        sum = sum / p1
	        s(j) = 0.5 * (s(j) + y_rng*sum/tnm)
	        it = 2*it
	    endif
	    if (j .ge. k) then
	        l = j - km
	        call polint(h(l), s(l), k, zero, oromb, del)
	        if (j .ge. 10) write(*,'(i5,g13.5)') j, oromb
c	        if (j .gt. 5) write(9,'(i5,g13.4)') j, oromb
	        if (abs(del) .lt. geps*abs(oromb)) return
	        if (oromb .eq. zero) return
	    endif
	    s(j+1) = s(j)
	    h(j+1) = 0.25 * h(j)
40	continue	
        print *, 'Error:  ',oromb,del
        PAUSE 'oromb:  Too many outer steps.'

	end

c-----------------------------------------------------------------------
	function inner (E, rl, ru, S1)

c---	"INNERmost integral"
c
c	Return the value of the inner integral over r, as a function
c	of E and global variables Sl, mu, and phi.  This inner 
c	integral may have two parts.

c+++  Arguments:
	real*8 inner, E, rl, ru, S1

c+++  Locals:
	real*8 pi4
	parameter (pi4 = 4.d0*3.14159265358d0)
	real*8 r1, r2
	real*8 E_dist, i_inner
	external i_inner

c+++  Globals:
	real*8 gSl, gmu, gphi, geps
	common /docb/ gSl, gmu, gphi, geps

	r1 = rl
	r2 = min( ru, sqrt(E/(pi4*S1)) )
	if (r2 .gt. r1) then
	    call romb(i_inner, r1, r2, inner, geps)
	else
	    inner = 0.
	endif

c	if (gmu.gt..23) write(9,'(a,5g12.4)') 
c     .       '--> ',E,r1,r2,inner,E_dist(E)*inner
	inner = E_dist(E) * inner

	return
	end

c-----------------------------------------------------------------------
	function i_inner (r)

c---	"Integrand for INNERmost integral"
c
c	Return the innermost integral, over r, required for dn_dO.

c+++  Arguments:
	real*8 i_inner, r

c+++  Locals:
	real*8 one
	parameter (one = 1.)
	real*8 s, x, y, z
	real*8 dn_dV

c+++  Globals:
	real*8 gSl, gmu, gphi, geps
	common /docb/ gSl, gmu, gphi, geps
	real*8 R0
	common /Rcb/ R0

c---  Calculate GC coordinates for the relevant point.
	s = sqrt(one - gmu*gmu)
	x = R0 - r*gmu
	y = r * s * sin(gphi)
	z = r * s * cos(gphi)

c---  Here's the integrand: just the number in dr.
	i_inner = r**2 * dn_dV(x,y,z)
c	if (i_inner .eq. 0.) write(*,'(7g11.3)') gSl,gmu,gphi,
c     +       sqrt(r*r+R0**2-2.*R0*r*gmu),dn_dV(x,y,z)
	return
	end

c=======================================================================
	subroutine r_range (mu, phi)

c	"Radial RANGE 1"
c
c	Range of r inside the spatial dist'n in the specified drxn,
c	for the 1st intersection of the l.o.s. with the distribution.

c+++  Arguments:
	real*8 mu, phi, rl, ru

c+++  Locals:
	real*8 one
	parameter (one = 1.d0)
	real*8 muc, outer, inner

c+++  Globals:
	real*8 R0
	common /Rcb/ R0
	integer dmmax, dpmax
	parameter (dmmax = 2, dpmax = 5)
	integer dmodel, ndpar(dmmax)
	real*8 dpar(dpmax)
	common /dmodcb/ dpar, ndpar, dmodel

	    inner = dpar(4)
	    outer = dpar(5)
c>>>  R0 = 0 is a simple special case.
	    if (R0 .eq. 0.) then
	        rl = inner
	        ru = outer
c>>>  If inside the inner boundary, there is one simple intsxn.
	    else if (R0 .le. inner) then
	        rl = R0 * (sqrt(mu*mu-one+(inner/R0)**2) + mu)
	        ru = R0 * (sqrt(mu*mu-one+(outer/R0)**2) + mu)
c>>>  If between the boundaries, rl=0, and ru depends on whether
c>>>  the los hits the inner boundary or not.
	    else if (R0 .gt. inner .and. R0 .lt. outer) then
	        rl = 0.
	        muc = sqrt(one - (inner/R0)**2)
	        if (mu .gt. muc) then
	            ru = R0 * (mu - sqrt(mu*mu-one+(inner/R0)**2))
	        else
	            ru = R0 * (sqrt(mu*mu-one+(outer/R0)**2) + mu)
	        endif
c>>>  Finally, if outside the outer boundary, there is no intsxn if
c>>>  mu is too small.  Otherwise, rl is set by the first intsxn with
c>>>  the outer boundary, and ru depends on whether the los hits
c>>>  the inner boundary.
	    else
	        muc = sqrt(one - (outer/R0)**2)
	        if (mu .gt. muc) then
	            rl = R0 * (mu - sqrt(mu*mu-one+(outer/R0)**2))
	            muc = sqrt(one - (inner/R0)**2)
	            if (mu .gt. muc) then
	                ru = R0 * (mu - sqrt(mu*mu-one+(inner/R0)**2))
	            else
	                ru = R0 * (mu + sqrt(mu*mu-one+(outer/R0)**2))
	            endif
	        else
	            rl = 0.
	            ru = 0.
	        endif
	    endif


	write(9,'(a,4(1pg13.4))') 'spherical r1: ',mu,phi,rl,ru
	return
	end

c***********************************************************************
      SUBROUTINE ROMB(FUNC,A,B,SS,eps)

        external FUNC
        real*8 func, a, b, ss, dss, eps
        integer j, l, jmax, jmaxp, k, km
      PARAMETER(JMAX=20,JMAXP=JMAX+1,K=3,KM=2)
      real*8 h(jmaxp), s(jmaxp)

        if(a.eq.b) then
           ss=0.
           return
        endif

      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          L=J-KM
          CALL POLINT(H(L),S(L),K,0.,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
          if (SS .eq. 0.) return
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
11    CONTINUE
      PAUSE 'ROMB:  Too many steps.'
      END
 
