c-----------------------------------------------------------------------
c	This files contains subprograms used to evaluate the
c	likelihood function at one point in parameter space.
c
c	01 Dec 93  TJL  (Adapted from earlier files)
c	Last Mod:  7 Feb 95 --- modified to better handle SSC model
c-----------------------------------------------------------------------
	subroutine newpar

c---	"NEW PARameters"
c
c	Reset global variables to indicate that the likelihood must
c	be calculated because parameters have changed.
c
c	If in addition the data or detector charactersitics have
c	changed, call lkreset instead.

c+++  Globals:
	include 'modchar.cb'

c---  Note that R_eff and the dR/dF grid must be calculated.
	didR = .false.

c---  Note that L itself must also be recalculated.
	didl = .false.

	return
	end

c-----------------------------------------------------------------------
	subroutine lkreset

c---	"LiKelihood RESET"
c
c	Reset global variables to indicate that the likelihood must
c	be calculated from scratch.

c+++  Globals:
	include 'modchar.cb'

c---  Note that R_eff and the dR/dF grid must be calculated.
	didR = .false.

c---  Note that L itself must also be recalculated.
	didl = .false.

c---  Set the "previous" values of nu in the compound model to
c---  ensure it will be recalculated.
	pnu_h = -1.
	pnu_c = -1.

	return
	end

c-----------------------------------------------------------------------
	subroutine set_ltype ( marg )

c---	"SET Likelihood TYPE"
c
c	If marg = .true., we calculate the marginal likelihood for
c	shape parameters.  If marg = .false., we calculate the
c	full likelihood.

c+++  Arguments:
	logical marg

c+++  Globals:
	include 'modchar.cb'

	margl = marg
	call lkreset

	return
	end

c-----------------------------------------------------------------------
	function ljlike ()

c---	"Log Joint quasi LIKElihood"

c+++  Arguments:
	real*8 ljlike

c+++  Globals:
	include 'modchar.cb'
	include 'data.cb'

c+++  Locals:
	real*8 tiny, day
	parameter (tiny = 1.d-200, day = 86400.)
	integer i
	real*8 dR

c+++  Functions:
	real*8 T_obs

c---  Don't waste any time repeating a calculation!
	if (didl) then
	    ljlike = llcur
	    return
	endif

c---  We haven't yet coded the full likelihood for superposed models.
	if (.not. margl .and. curtype .eq. 3) 
     *      call quit('Cannot do full L for superposed models!')

c---  Set the parameters for all models, and calculate normalizing
c---  integrals.
	call setmod
	call lksetup

c---  Calculate the log likelihood.
	ljlike = 0.
	do 500 i=1, nb
	    if (useb(i)) then
	        call bquad(i, dR)
	        if (dR .gt. tiny) then
	            ljlike = ljlike + log(dR)
	        else
	            ljlike = ljlike - 500.
	        endif
	    endif
500	continue
	ll_prod = ljlike
	ll_exp = T_obs()*day*R_eff(1)
	if (.not. margl) ljlike = ljlike - ll_exp

c---  Adjust some globals to reflect having done the calculation.
	llcur = ljlike
	didl = .true.

	return
	end

c-----------------------------------------------------------------------
	subroutine bquad (i, dR)

c---	"Burst QUADrature"
c
c	Return the quadrature of the total differential rate with
c	the individual burst likelihood function with the 
c	characteristics of burst i.

c+++  Arguments:
	integer i
	real*8 dR

c+++  Globals:
	include 'modchar.cb'
	include 'data.cb'

c+++  Locals:
	real*8 one
	parameter (one = 1.)
	real*8 fac, R
	real*8 T_obs

c---  For most models, this is straightforward:  Call the appropriate
c---  subroutine for isotropic or anisotropic models.
	if (curtype .eq. 1) then
	    call F_quad(Fvals(i), Fsigs(i), dR)
	    if (margl) dR = dR / R_eff(1)
	else if (curtype .eq. 2) then
	    call Flb_quad (Fvals(i), Fsigs(i), lvals(i), bvals(i), 
     *               kvals(i), absc(1,i), wts(1,i), dR)
	    if (margl) dR = dR / R_eff(1)
	else
	    if (.not. margl) 
     *          call quit('Can only calculate marginal for ssc!')

c---  But for the compound SSC model, it is more complicated.
c---  Calculate the cosmo & halo parts separately, since the
c---  required integrals are of different dimensions.  Don't
c---  waste time on an integral if it is not needed (f_h = 0 or 1).
c---  Also, the f_h=1 case is special because of how we
c---  parameterize the weights.
	  if (f_h .eq. 1.) then
	    if (.not. didh) then
	        call Flb_quad (Fvals(i), Fsigs(i), lvals(i), bvals(i), 
     *               kvals(i), absc(1,i), wts(1,i), dR_h(i))
	    endif
	    dR = dR_h(i) / R_eff(2)
	  else
	    if (.not. didc) then
	        call F_quad (Fvals(i), Fsigs(i), dR_c(i))
	    endif
	    if (f_h .eq. 0.) then
	        dR = dR_c(i) / R_eff(1)
	    else
	        fac = f_h / (one - f_h)
	        if (.not. didh) then
	          call Flb_quad (Fvals(i), Fsigs(i), lvals(i), bvals(i),
     *                 kvals(i), absc(1,i), wts(1,i),  dR_h(i))
	        endif
	        R = R_eff(1) / (one - f_h)
	        dR = (dR_c(i) + fac*R_eff(1)*dR_h(i)/R_eff(2)) / R
	    endif
	  endif
	endif

c---  Include the observing time factor.
	if (margl) dR = dR / T_obs()

	return
	end

c-----------------------------------------------------------------------
	subroutine set_Fgrid (F1, F2, F3, n1, n2, n3, nu)

c---	"SET F GRID"
c
c	Set the grid in flux used to calculate normalizing integrals
c	for the burst rate.
c
c	The grid uses n1 pts from F1 to the threshold of eta
c	  (these are used to interpolate dR/dF for isotropic models),
c	n2 points from the threshold to F2
c	and n3 points from F2 to F3.
c	An additional segment from F3 to infinity can be added
c	with set_infF.
c
c	On return, nu is the # of bursts with peak flux estimates
c	in [Fth,Fmax], where Fmax = F3 or infinity.

c+++  Arguments:
	real*8 F1, F2, F3
	integer n1, n2, n3, nu

c+++  Globals:
	include 'modchar.cb'
	include 'data.cb'

c+++  Locals:
	real*8 w, dF
	logical Flog
	integer i, nn, indx

c---  Set the grid.  First check the # of points requested.
	nFg = abs(n1) + abs(n2) + abs(n3)
	nF1 = abs(n1)
	if (nFg .gt. nFmax) call quit('Too many points in F grid!')

c---  Now do up to the threshold.
	if (F1 .ge. F_th) call quit('Lo F limit must be < F_th!')
	F_lo = F1
	F_hi = F3
	call setstp (F_lo, F_th, n1, Flog, dF)
	do 20 i=1, abs(n1)-1
	    call stpinf(i, F_lo, dF, Flog, Fgrid(i), w)
20	continue

c---  Now do from the threshold to F2.
	if (F2 .le. F_th) call quit('F2 must be > F_th!')
	Fgrid(abs(n1)) = F_th
	if (n2 .lt. 0) then
	    nn = n2 - 1
	else
	    nn = n2 + 1
	endif
	indx = abs(n1) - 1
	call setstp (F_th, F2, nn, Flog, dF)
	do 30 i=2, abs(n2)+1
	    call stpinf(i, F_th, dF, Flog, Fgrid(indx+i), w)
30	continue

c---  Now do from F2 to F3.
	if (n3 .lt. 0) then
	     nn = n3 -1
	else
	    nn = n3 + 1
	endif
	indx = abs(n1) + abs(n2) - 1
	call setstp (F2, F_hi, nn, Flog, dF)
	do 50 i=2, abs(n3)+1
	    call stpinf(i, F2, dF, Flog, Fgrid(indx+i), w)
50	continue

c---  Find out how many bursts lie above the cutoff.
	call flagb(nu)

	return
	end

c-----------------------------------------------------------------------
	subroutine set_infF (ch, nu)

c---	"SET INFinite F flag"
c
c	Set a flag indicating whether to tag a powerlaw piece onto
c	the F grid that goes up to infinite flux.

c+++  Arguments:
	character*(*) ch
	integer nu

c+++  Globals:
	include 'modchar.cb'

	infF = (ch(1:1) .eq. 't') .or. (ch(1:1) .eq. 'T')
	call flagb(nu)

	return
	end

c-----------------------------------------------------------------------
	subroutine set_monF (i)

c---	"SET MONitor interval for F grid"
c
c	Set how often during every computation of an F grid to
c	let the user know how things are going.

c+++  Arguments:
	integer i

c+++  Globals:
	include 'modchar.cb'

	monF = i

	return
	end

c-----------------------------------------------------------------------
	function dR_dF (F)

c---	"d(Rate)/d(Flux)"
c
c	Returns the differential burst rate per unit Flux by
c	interpolating the Fgrid.  This is only good for
c	isotropic models.
c
c	For compound models, this returns only the first
c	component.

c+++  Arguments:
	real*8 dR_dF, F

c+++  Globals:
	include 'modchar.cb'

c+++  Locals:
	integer nlo, nhi, sign, new
	real*8 slope

	if (F .lt. F_lo) then
	    dR_dF = 0.
	    return
	endif

c---  Interpolate/extrapolate the cosmo part with a powerlaw.
c---  First, locate F in Fgrid.
        nlo = 1
        nhi = nFg
        sign = (Fgrid(nFg)-Fgrid(1)) / abs(Fgrid(nFg)-Fgrid(1))
20      if (nhi - nlo .gt. 1) then
            new = (nhi + nlo) / 2
            if ((Fgrid(new)-F)/sign .gt. 0.) then   
                nhi = new
            else
                nlo = new
            endif
            go to 20
        endif

c---  Find the powerlaw slope.
	slope = log(dRdF(nhi,1)/dRdF(nlo,1)) /
     +          log(Fgrid(nhi)/Fgrid(nlo))

c---  Interpolate/extrapolate with this slope.
	dR_dF = dRdF(nlo,1) * (F/Fgrid(nlo))**slope

	return
	end

c-----------------------------------------------------------------------
	function dR_dF_2 (F)

c---	"d(Rate)/d(Flux) for component 2"
c
c	Returns the differential burst rate per unit Flux by
c	interpolating the Fgrid for the *second* component of
c	compound models.

c+++  Arguments:
	real*8 dR_dF_2, F

c+++  Globals:
	include 'modchar.cb'

c+++  Locals:
	integer nlo, nhi, sign, new
	real*8 slope

	if (F .lt. F_lo) then
	    dR_dF_2 = 0.
	    return
	endif

c---  Interpolate/extrapolate the cosmo part with a powerlaw.
c---  First, locate F in Fgrid.
        nlo = 1
        nhi = nFg
        sign = (Fgrid(nFg)-Fgrid(1)) / abs(Fgrid(nFg)-Fgrid(1))
20      if (nhi - nlo .gt. 1) then
            new = (nhi + nlo) / 2
            if ((Fgrid(new)-F)/sign .gt. 0.) then   
                nhi = new
            else
                nlo = new
            endif
            go to 20
        endif

c---  Find the powerlaw slope.
	slope = log(dRdF(nhi,2)/dRdF(nlo,2)) /
     +          log(Fgrid(nhi)/Fgrid(nlo))

c---  Interpolate/extrapolate with this slope.
	dR_dF_2 = dRdF(nlo,2) * (F/Fgrid(nlo))**slope

	return
	end

c-----------------------------------------------------------------------
	subroutine lksetup

c---	"LiKelihood SETUP"
c
c	Setup global values needed for evaluation of the likelihood
c	function.  This involves calculating the normalization for
c	the burst rate and making a grid to interpolate dR/dF for 
c	isotropic models.
c
c	For halo models, the normalizing constant can optionally be
c	interpolated from a grid read by the set.halo command.
c
c	For the ssc model, we do not waste time recalculating a
c	model component that has not changed.

c+++  Globals:
	include 'main.cb'
	include 'modchar.cb'

c+++  Locals:
	real*8 zero
	parameter (zero = 0.)
	integer n, nF, nelem, nrd
	real*8 ls, vec(nFmax)
	logical ok
	real*8 dR_dFdO, eta_F, R_sch, eta_max
	logical use_hg

c---  Don't waste any time repeating a calculation!
	if (didR) return

c---  Complain if an Fgrid hasn't been specified.
	if (nFg .eq. 0) call quit('Must specify F grid!')

c---  If the model is anisotropic, complain if angular quadrature
c---  points have not been specified.
	if (curtype .eq. 2) call chk_sky
	if (curtype .eq. 3) call chk_sky

c---  See if this is a compound model.
	if (curtype .eq. 3) then
	    nelem = 2
	else
	    nelem = 1
	endif

c---  Loop over the elements of the model.
	do 200 n=1, nelem
	    comp = n

c===  The ssc compound model is treated specially.
c===  We don't recalculated a component if it was previously
c===  calculated and its parameters have not changed.
c===  Also, we don't calculate a component if it is not needed
c===  (ie, if its fraction is 0).  Finally, we don't do the
c===  halo part from scratch if we can get it from a grid.
	    if (curmod .eq. 8) then
	        if (n .eq. 1 .and. f_h .eq. 1.) then
	            R_eff(n) = 0.
	            pnu_c = -1.
                    goto 200
	        else if (n .eq. 2 .and. f_h .eq. 0.) then
	            R_eff(n) = 0.
	            pnu_h = -1.
                    goto 200
	        else if (n .eq. 2 .and. use_hg()) then
	            R_eff(n) = R_sch(nu_h)
	            goto 200
	        else if (n .eq. 1 .and. didc .and. nu_c .eq. pnu_c)then
	            goto 200
	        else if (n .eq. 2 .and. didh .and. nu_h .eq. pnu_h)then
	            goto 200
	        else if (n .eq. 1) then
	            pnu_c = nu_c
	        else if (n .eq. 2) then
	            pnu_h = nu_h
	        endif
	    endif

c===  If this is the file spec model, read in the Fgrid and the rates.
	if (curmod .eq. 12) then
	    call rdcol (flun2, 'iso_rate.dat','iso', '<', '>', nFmax,
     *             1, Fgrid, nrd, ok)
	    if (.not. ok .or. nrd .ne. nFg) then
	        print *, 'Data read from iso_rate.dat: ',nrd
	        print *, Fgrid(1), nrd
	        call quit('lkcalc: ISO_RATE grid mismatch!')
	    endif
	    call rdcol (flun2, 'iso_rate.dat','iso', '<', '>', nFmax,
     *             4, vec, nrd, ok)
	    do 20 nF=1, nFg
	        dRdF(nF,n) = vec(nF)
20	    continue
	    call cpar
	    goto 85
	endif

c===  Otherwise, evaluate dR/dF on the grid.  If the model is 
c===  anisotropic, this requires angular integrals performed by 
c===  sky_int.
c===  NOTE:  Since sky_int includes the angular part of the
c===  1B eta function, we must divide by eta_max.
c===  For isotropic models, we do a very low flux part of the grid
c===  to let us later interpolate dR/dF.
	    do 80 nF=1, nFg
	        if (monF .gt. 0 .and. mod(nF,monF) .eq. 0) then
	          print *,'  Evaluating F grid pt ',nF,' of ',nFg
	        endif
	        if (curtype .eq. 1) then
	            dRdF(nF,n) = pi4 * dR_dFdO(Fgrid(nF),zero,zero)
	        else if (curtype .eq. 2) then
c	            if (nF .lt. nF1)
c     *               print *, '**** Wasting time below nF1!!!!!! ****'
c	            if (nF .lt. 0) then
	            if (nF .lt. nF1) then
	                dRdF(nF,n) = 0.
	            else
	                call sky_int(Fgrid(nF), dRdF(nF,n),n)
	                dRdF(nF,n) = dRdF(nF,n) / eta_max()
	            endif
	        else if (curtype .eq. 3) then
	            if (n .eq. 1) then
	                dRdF(nF,n) = pi4 * dR_dFdO(Fgrid(nF),zero,zero)
	            else
	                if (nF .lt. nF1) then
	                    dRdF(nF,n) = 0.
	                else
	                    call sky_int(Fgrid(nF), dRdF(nF,n),n)
	                    dRdF(nF,n) = dRdF(nF,n) / eta_max()
	                endif
	            endif
	        endif
80	    continue

c===  Now that we have the grid, use it to calculate integrals
c===  assuming powerlaw behavior between grid points.  But since
c===  the eta_F factor may make eta_F*dR vanish at some points, 
c===  some eta_F*dR intervals must use the trapezoid rule instead.
c===  Start the calculation at the threshold grid point.
85	    continue
	    R_eff(n) = 0.
	    do 100 nF=nF1, nFg-1
	        if (eta_F(Fgrid(nF+1))*dRdF(nF+1,n)*
     *              eta_F(Fgrid(nF))*dRdF(nF,n) .gt. 0.) then
	            ls = ( log(eta_F(Fgrid(nF+1))*dRdF(nF+1,n)/
     +                         (eta_F(Fgrid(nF))*dRdF(nF,n))) / 
     +                     log(Fgrid(nF+1)/Fgrid(nF)) )  +  1.
	            if (ls .ne. 0.) then
	                R_eff(n) = R_eff(n) + 
     +                           eta_F(Fgrid(nF))*dRdF(nF,n)*Fgrid(nF)*
     +                           ((Fgrid(nF+1)/Fgrid(nF))**ls - 1.) / ls
	            else
	                R_eff(n) = R_eff(n) + 
     +                           eta_F(Fgrid(nF))*dRdF(nF,n)*Fgrid(nF)*
     +                           log(Fgrid(nF+1)/Fgrid(nF))
	            endif
c	        write(9,'(2i3,4(1pg12.4))') n,nF,Fgrid(nF),
c     *               eta_F(Fgrid(nF))*dRdF(nF,n),
c     *               ls, R_eff(n)
	        else
	            R_eff(n) = R_eff(n) + 
     +                         0.5*(eta_F(Fgrid(nF))*dRdF(nF,n)
     +                         +  eta_F(Fgrid(nF+1))*dRdF(nF+1,n)) *
     +                         (Fgrid(nF+1)-Fgrid(nF))
	        endif
100	    continue

c===  Alert stdout if we are not near the homogeneous part at
c===  the end of the grid, unless dRdF vanishes, or unless
c===  we are using a phenomenological model.
	    if (dRdF(nFg,n)*dRdF(nFg-1,n) .gt. 0.) then
	        ls = ( log(dRdF(nFg,n)/dRdF(nFg-1,n)) / 
     +                 log(Fgrid(nFg)/Fgrid(nFg-1)) )
	        if (abs(ls+2.5) .gt. 0.02 .and. 
     +              (curmod .gt. 3 .and. curmod .ne. 13 .and.
     +               curmod .ne. 9)) then
	            write(*,'(a,1pg12.4)') 'ls = ', ls
	            print *, 'Note:  Grid ends before homogeneity...'
	        endif
	    endif

c===  If requested, add on what you get if you extrapolate to
c===  infinite F (but not if dRdF bottomed out!).
	    if (infF .and. (dRdF(nFg,n)*dRdF(nFg-1,n) .gt. 0.)) then
	        ls = ( log(dRdF(nFg,n)/dRdF(nFg-1,n)) / 
     +                 log(Fgrid(nFg)/Fgrid(nFg-1)) ) + 1.
	        if (ls .ge. 0.) then
	            write(*,'(a,1pg12.4)') 'ls = ',ls-1.
	            print *, Fgrid(nFg-1), Fgrid(nFg)
	            print *, dRdF(nFg-1,n), dRdF(nFg,n)
	            call quit('R_tot integral diverges!')
	        endif
	        R_tot(n) = R_tot(n) - dRdF(nFg,n)*Fgrid(nFg) / ls
	        ls = ( log(eta_F(Fgrid(nFg))*dRdF(nFg,n)/
     +                     (eta_F(Fgrid(nFg-1))*dRdF(nFg-1,n))) / 
     +                 log(Fgrid(nFg)/Fgrid(nFg-1)) ) + 1.
	        if (ls .ge. 0.) then
	            write(*,'(a,1pg12.4)') 'ls = ',ls-1.
	            call quit('R_eff integral diverges!')
	        endif
	        R_eff(n) = R_eff(n) - 
     +                     eta_F(Fgrid(nFg))*dRdF(nFg,n)*Fgrid(nFg)**ls
     +                     / ls
	    endif

c===  Go for the next model element, if compound.
200	continue
c	write(*,'(a,2(1pg14.6))') 'R1, 2 = ',R_eff(1),R_eff(2)

c---  Note that the total rates are calculated, and that's it.
	didR = .true.

	return
	end

c-----------------------------------------------------------------------
	subroutine globalq (n, R, Re)

c---	"GLOBAL Quantities"
c
c	Return the values of integrals of the rate based on model n.
c	If n=0, return grand totals.

c+++  Arguments:
	real*8 R, Re
	integer n

c+++  Globals:
	include 'modchar.cb'

	if (n .lt. 0 .or. n .gt. nmod) 
     +      call quit('globalq:  Illegal model #!')
	R = R_tot(n)
	Re = R_eff(n)

	return
	end

c-----------------------------------------------------------------------
	subroutine wglobalq (lun)

c---	"Write GLOBAL Quantities"
c
c	Write global quantities for each model element and for the
c	sum to lun.  If lun=-1, use stdout.

c+++  Arguments:
	integer lun

c+++  Globals:
	include 'modchar.cb'

c+++  Locals:
	integer n

	if (lun .ge. 0) then
	    write(lun,20) R_tot(0), R_eff(0)
	else
	    write(*,20) R_tot(0), R_eff(0)
	endif
	if (nmod .gt. 1) then
	    if (lun .ge. 0) then
	        write(lun,30) nmod
	        do 10 n=1, nmod
	            write(lun,40) n,R_tot(n),R_eff(n)
10	        continue
	        write(lun,'(a)')
	    else
	        write(*,30) nmod
	        do 15 n=1, nmod
	            write(*,40) n,R_tot(n),R_eff(n)
15	        continue
	        write(*,'(a)')
	    endif
	endif
20	format('Global quantities for current model:',/,
     *         '  R_tot = ',1pg12.4,/,
     *         '  R_eff = ',1pg12.4,/)
30	format('Contributions from ',i2,
     *         ' model elements (#, R_tot, R_eff):')
40	format('  element ',i2,5(1pg12.4))
	return
	end

c-----------------------------------------------------------------------
	subroutine moments (lun, fname)

c	"MOMENTS"
c
c	Find angular moments vs. flux, and total moments, for
c	the current model, and write them out to the specified file.

c+++  Arguments:
	integer lun
	character*(*) fname

c+++  Globals:
	include 'main.cb'
	include 'modchar.cb'

c+++  Locals:
	real*8 zero, third
	parameter (zero = 0., third = 1.d0/3.d0)
	integer nF
	real*8 dD(nFmax), dQ(nFmax), gamma(nFmax), Dtot, Qtot, gammab
	real*8 rat, frat, lfrat, edR, edR1, mrat, dF, e, e1
	real*8 ls, R, dR
	real*8 eta_F, eta_max

c---  Complain if an Fgrid hasn't been specified.
	if (nFg .eq. 0) call quit('Must specify F grid!')

c---  If the model is anisotropic, complain if angular quadrature
c---  points have not been specified.
	call setmod
	if (curtype .eq. 2) call chk_sky
	if (curtype .eq. 3) call chk_sky

c---  Note we're messing with the flux grid.
	didR = .false.

c---  See if this is a compound model or isotropic model.
	if (curtype .eq. 3) then
	    call quit('Cannot do moments for compound models!')
	else if (curtype .eq. 1) then
	    call quit('Uhhh... this model is isotropic!')
	endif


c===  Evaluate moments on the flux grid.  
c===  NOTE:  Since sky_int includes the angular part of the
c===  1B eta function, we must divide by eta_max, because we
c===  will multiply by eta_F, which includes the average exposure
c===  factor.
	    do 80 nF=nF1, nFg
	        if (monF .gt. 0 .and. mod(nF,monF) .eq. 0) then
	          print *,'  Evaluating F grid pt ',nF,' of ',nFg
	        endif
	        call sky_mom(Fgrid(nF), dRdF(nF,1), dD(nF), dQ(nF))
	        dRdF(nF,1) = dRdF(nF,1) / eta_max()
	        dD(nF) = dD(nF) / eta_max()
	        dQ(nF) = dQ(nF) / eta_max()
c	        dD(nF) = -dRdF(nF,1)*0.01266
c	        dQ(nF) = dRdF(nF,1)*0.32848333
80	    continue

c===  Now that we have the grid, use it to calculate integrals
c===  assuming powerlaw behavior between grid points.  But since
c===  the eta_F factor may make eta_F*dR vanish at some points, 
c===  some eta_F*dR intervals must use the trapezoid rule instead.
c===  Start the calculation at the threshold grid point.
85	    continue
	    R_eff(1) = 0.
	    Dtot = 0.
	    Qtot = 0.
	    do 100 nF=nF1, nFg-1
	        e = eta_F(Fgrid(nF))
	        e1 = eta_F(Fgrid(nF+1))
	        edR = e*dRdF(nF,1)
	        edR1 = e1*dRdF(nF+1,1)
	        dF = Fgrid(nF+1)-Fgrid(nF)
	        if (edR1*edR .gt. 0.) then
	            rat = edR1 / edR
	            frat = Fgrid(nF+1) / Fgrid(nF)
	            lfrat = log(frat)
	            ls = log(rat)/lfrat  +  1.
	            if (ls .ne. 0.) then
	                R_eff(1) = R_eff(1) + edR*Fgrid(nF)*
     +                           (frat**ls - 1.) / ls
	            else
	                R_eff(1) = R_eff(1) + edR*Fgrid(nF)*log(frat)
	            endif
	            if (dD(nF)*dD(nF+1) .le. 0.) then
	                Dtot = Dtot + 0.5*(e*dD(nF)+e1*dD(nF+1))*dF
	            else
	                mrat =  e1*dD(nF+1) / (e*dD(nF))
	                ls = log(mrat)/lfrat + 1.
	                if (ls .ne. 0.) then
	                    Dtot = Dtot + e*dD(nF)*Fgrid(nF)*
     +                             (frat**ls - 1.) / ls
	                else
	                    Dtot = Dtot + e*dD(nF)*Fgrid(nF)*log(frat)
	                endif
	            endif
	            if (dQ(nF) .le. 0. .or. dQ(nF+1) .le. 0.) then
	                Qtot = Qtot + 0.5*(e*dQ(nF)+e1*dQ(NF+1))*dF
	            else
	                mrat =  e1*dQ(nF+1) / (e*dQ(nF))
	                ls = log(mrat)/lfrat + 1.
	                if (ls .ne. 0.) then
	                    Qtot = Qtot + e*dQ(nF)*Fgrid(nF)*
     +                               (frat**ls - 1.) / ls
	                else
	                    Qtot = Qtot + e*dQ(nF)*Fgrid(nF)*log(frat)
	                endif
	            endif
	        else
	            R_eff(1) = R_eff(1) + 0.5*(edR + edR1) * dF
	            Dtot = Dtot + 0.5*(e*dD(nF) + e1*dD(NF+1))*dF
	            Qtot = Qtot + 0.5*(e*dQ(nF) + e1*dQ(NF+1))*dF
	        endif
100	    continue

c===  Alert stdout if we are not near the homogeneous part at
c===  the end of the grid, unless dRdF vanishes, or unless
c===  we are using a phenomenological model.
	    if (dRdF(nFg,1)*dRdF(nFg-1,1) .gt. 0.) then
	        ls = ( log(dRdF(nFg,1)/dRdF(nFg-1,1)) / 
     +                 log(Fgrid(nFg)/Fgrid(nFg-1)) )
	        if (abs(ls+2.5) .gt. 0.02 .and. 
     +              (curmod .gt. 3 .and. curmod .ne. 13 .and.
     +               curmod .ne. 9)) then
	            write(*,'(a,1pg12.4)') 'ls = ', ls
	            print *, 'Note:  Grid ends before homogeneity...'
	        endif
	    endif

c===  If requested, add on what you get if you extrapolate to
c===  infinite F (but not if dRdF bottomed out!).
	    if (infF .and. (dRdF(nFg,1)*dRdF(nFg-1,1) .gt. 0.)) then
	        e = eta_F(Fgrid(nFg))
	        e1 = eta_F(Fgrid(nFg-1))
	        edR = e*dRdF(nFg,1)
	        edR1 = e1*dRdF(nFg-1,1)
	        rat = edR / edR1
	        frat = Fgrid(nFg)/Fgrid(nFg-1)
	        lfrat = log(frat)
	        ls = ( log(rat) / lfrat ) + 1.
	        if (ls .ge. 0.) then
	            write(*,'(a,1pg12.4)') 'ls = ',ls-1.
	            write(*,'(4(1pg12.4))') edR, edR1, dRdF(nFg,1),
     *                 dRdF(nFg-1,1)
	            call quit('R_eff integral diverges!')
	        endif
	        R_eff(1) = R_eff(1) - edr*Fgrid(nFg)**ls / ls
	        mrat = e*dD(nFg) / (e1*dD(nFg-1))
	        ls = log(mrat) / lfrat + 1.
	        if (ls .ge. 0.) then
	            write(*,'(a,1pg12.4)') 'ls = ',ls-1.
	            call quit('Dtot integral diverges!')
	        endif
	        Dtot = Dtot - e*dD(nFg)*Fgrid(nFg)**ls / ls
	        mrat = e*dQ(nFg) / (e1*dQ(nFg-1))
	        ls = log(mrat) / lfrat + 1.
	        if (ls .ge. 0.) then
	            write(*,'(a,1pg12.4)') 'ls = ',ls-1.
	            call quit('Qtot integral diverges!')
	        endif
	        Qtot = Qtot - e*dQ(nFg)*Fgrid(nFg)**ls / ls
	    endif

c---  Now find the log slope of the cumulative.
	R = R_eff(1)
	gammab = 0.
	do 200 nF=nF1, nFg-1
	    if (R .eq. 0.) then
	        gamma(nF) = gamma(nF-1)
	        goto 200
	    endif
	    edR = eta_F(Fgrid(nF))*dRdF(nF,1)
	    edR1 = eta_F(Fgrid(nF+1))*dRdF(nF+1,1)
	    dF = Fgrid(nF+1)-Fgrid(nF)
	    if (edR1*edR .gt. 0.) then
	        rat = edR1 / edR
	        frat = Fgrid(nF+1) / Fgrid(nF)
	        lfrat = log(frat)
	        ls = log(rat)/lfrat  +  1.
	        if (ls .ne. 0.) then
	            dR = edR*Fgrid(nF)*(frat**ls - 1.) / ls
	        else
	            dR = edR*Fgrid(nF)*log(frat)
	        endif
	    else
	        dR = 0.5*(edR + edR1) * dF
	    endif
	    if (dR .ge. R) then
	        R = 0.
	        gamma(nF) = gamma(nF-1)
	        goto 200
	    endif
	    gamma(nF) = - log(1.-dR/R)/log(Fgrid(nF+1)/Fgrid(nF))
	    gammab = gammab + 0.5*gamma(nF)*(edR + edR1)*dF
	    R = R - dR
200	continue
	gamma(nFg) = gamma(nFg-1)

c---  Normalize the totals.
	Dtot = Dtot / R_eff(1)
	Qtot = (Qtot / R_eff(1)) - third
	gammab = gammab / R_eff(1)

c---  Open the file.
	call openf(lun, fname)

c---  Write the totals, then the running moments.
	write(lun, 500) R_eff(1), Dtot, Qtot, gammab
500	format('R_eff  ', 1pg14.4,/,
     *         'D.Q    ',2(1pg14.4),/,
     *         'gammab ',1pg14.4,/)

	do 520 nF=nF1, nFg
	    write(lun,540) Fgrid(nF), eta_F(Fgrid(nF)), dRdF(nF,1), 
     *         dD(nF)/dRdF(nF,1), (dQ(nF)/dRdF(nF,1))-third, gamma(nF)
520	continue
540	format(6(1pg12.4))
	close(lun)

	return
	end

c-----------------------------------------------------------------------
        subroutine modout(lun, termdoc)

c---    This subroutine writes out the parameter and symbol values
c       for all the formulae making up the current model.  The
c       results are written to the file attached to lun, and to
c       the terminal, if requested.

c+++  Arguments:
        integer lun
        logical termdoc

c+++  Globals:
	include 'data.par'
	include 'data.cb'
        include 'modsyms.cb'
        include 'modchar.cb'

c+++  Locals:
	real*8 day, year, kpc3, Mpc2
	parameter (day = 86400., year = 3.1557d7)
	parameter (kpc3 = 2.938000d64, Mpc2 = 9.521409d48)
	integer i, l
	real*8 value, tval, scale
	character*80 line
	integer nwlen

c+++  Functions:
	real*8 T_obs

c---  Write out the formulae and their parameters and symbols.
        do 60 i=1, nform
            call forout(lun, i, termdoc, .false.)
60      continue

c---  Write out the values of any symbols.  If transformations have
c---  been specified, write transformed values, too.
        do 100 i=1, nsym
            value = symval(i)
            if (tsym(i)) then
                call trans(i, value, tval)
                write(line, 120) symbls(i), value, tval
            else
                write(line, 140) symbls(i), value
            endif
            l = nwlen(line)
            write(lun, '(a)') line(1:l)
            if (termdoc) write(*, '(a)') line(1:l)
100     continue
120     format(a8,' = ',1pg14.6,'  ('' = ',1pg14.6,')')
140     format(a8,' = ',1pg14.6)

c---  When we calculate the marginal likelihood for 
c---  shape parameters, write out the mode for the scale factor
c---  (using a Jeffreys prior).  We actually write exp(s), where
c---  s is the mode of p(log scale).
	if (margl) then
	    if (curtype .eq. 1) then
	        scale = nuse / (T_obs() * day * R_eff(1))
	    else if (curtype .eq. 2) then
	        scale = nuse / (T_obs() * day * R_eff(1))
	    else if (curtype .eq. 3) then
	        scale = -1.
	    endif
c	    write(lun,200) scale, scale*year*kpc3, 1.e9*scale*year
c	    write(lun,200) scale, scale*year
	    write(lun,200) scale
	    if (termdoc) write(*,200) scale
	endif
c200	format('scale    = ',1pg14.6,'  (= ',1pg14.6,' /(yr*kpc**3))',/,
c     *         '           ',14x,'  (= ',1pg14.6,' /(yr*Gpc**3))')
200	format('scale    = ',1pg14.6)

        write(lun,*)
        if (termdoc) write(*,*)

        return
        end

c-----------------------------------------------------------------------
        subroutine wmodel (lun, name)

c---	"Write MODEL"
c
c---    Write out the cosmo part in rdpar format, including 
c 	globals and the values in the current F grid.

c+++  Arguments:
	integer lun
	character*(*) name

c+++  Globals:
        include 'modchar.cb'
        include 'data.par'
        include 'modsyms.cb'
        include 'data.cb'

c+++  Locals:
        integer i, j, n, nF
        real*8 cum, c1, c2, l1, l2, g, ls, F, Fn, Fp, cl1, cl2, cg
        real*8 dR(nFmax), R(nFmax), ljl, fac, Rtot, arg1, arg2, B
        real*8 tiny
        parameter (tiny = 1.e-30)

c+++  Functions:
        real*8 eta_F, eta_max, ljlike
        integer nwlen
        logical use_hg

c---  We don't have all the available info for the SSC model if we've
c---  been using a grid for the halo R_eff values.
	if (curmod .eq. 8 .and. use_hg()) then
	    call quit('WMODEL:  Cannot write SSC model if using grid!')
	endif

c---  Also, parts of the "do 15" loop presume only 1 model element.
	if (nform .gt. 1) call quit('WMODEL:  nform too large!')

c---  Make sure the current model has been set up.  Initialize the
c---  various output arrays.
c	call setmod
c	call lksetup
	if (.not. didl) ljl = ljlike()
	do 10 nF=1, nFg
	    dR(nF) = 0.
	    R(nF) = 0.
10	continue

c---  Open the file.
	call openf(lun, name)

c---  First write a header identifying the model elements.  Write
c---  out globals for each element where appropriate.  Accumulate
c---  the total rate as we go.
	write(lun,'(a,i3)') 'elements  ',nform
	n = 0
	if (curmod .eq. 8) then
	    if (f_h .ne. 0. .and. f_h .ne. 1.) then
	        fac = f_h / (1. - f_h)
	    endif
	endif
        do 200 i=1, nform
            write(lun,'(a)') fornam(fornum(i))
                n = n + 1
                do 15 nF=1, nFg
                    if (curmod .ne. 8) then
                        dR(nF) = dR(nF) + dRdF(nF,n)
                    else if (f_h .eq. 0.) then
                        dR(nF) = dRdF(nF,1)
                    else 
                        if (nF .lt. nF1) then
	                    call sky_int(Fgrid(nF), dRdF(nF,2),2)
	                    dRdF(nF,2) = dRdF(nF,2) / eta_max()
                        endif
                        if (f_h .eq. 1.) then
                            dR(nF) = dRdF(nF,2)
                        else
                            dR(nF) = (dRdF(nF,1) + 
     *                               fac*dRdF(nF,2)*R_eff(1)/R_eff(2)) 
                        endif
                    endif
15	        continue
	        Rtot = 0.
	        do 18 nF = nF1, nFg-1
                    if (eta_F(Fgrid(nF))*dR(nF) .eq. 0.) then
                        R(nF) = R(nF) +
     +                          0.5*eta_F(Fgrid(nF+1))*dR(nF+1) *
     +                          (Fgrid(nF+1) - Fgrid(nF))
                    else
	              ls = ( log(eta_F(Fgrid(nF+1))*dR(nF+1)/
     +                           (eta_F(Fgrid(nF))*dR(nF))) / 
     +                       log(Fgrid(nF+1)/Fgrid(nF)) ) + 1.
	              if (ls .ne. 0.) then
	                  R(nF) = R(nF) + 
     +                            eta_F(Fgrid(nF))*dR(nF)*Fgrid(nF)*
     +                            ((Fgrid(nF+1)/Fgrid(nF))**ls - 1.) /ls
	              else
	                  R(nF) = R(nF) + 
     +                            eta_F(Fgrid(nF))*dR(nF)*Fgrid(nF)*
     +                            log(Fgrid(nF+1)/Fgrid(nF))
	              endif
	            endif
	            Rtot = Rtot + R(nF)
18	        continue
	        if (curmod .ne. 8) then
                    write(lun,180) R_eff(n)
                else
                    write(lun,180) Rtot
                endif

c===  Write power law model params.
            if (fornum(i) .eq. 1) then
                write(lun,20) (params(j,i),j=1,1)

c===  Write rat'l broken power law model params.
            else if (fornum(i) .eq. 2) then
                write(lun,30) (params(j,i),j=1,3)

c===  Write broken power law model params.
            else if (fornum(i) .eq. 3) then
                write(lun,40) (params(j,i),j=1,4)

c===  Write st'd candle cosmo model params.
            else if (fornum(i) .eq. 4) then
                write(lun,50) (params(j,i),j=1,7)

c===  Write power-law luminosity func'n cosmo model params.
            else if (fornum(i) .eq. 5) then
                write(lun,60) (params(j,i),j=1,10)

c===  Write st'd candle halo model params.
            else if (fornum(i) .eq. 6) then
                write(lun,70) (params(j,i),j=1,5)

c===  Write ower-law luminosity func'n halo model params.
            else if (fornum(i) .eq. 7) then
                write(lun,80) (params(j,i),j=1,6)

c===  Write superposed st'd candle model params.
            else if (fornum(i) .eq. 8) then
                write(lun,90) (params(j,i),j=1,8)

c===  Write peak duration broken power law model params.
            else if (fornum(i) .eq. 9) then
                write(lun,100) (params(j,i),j=1,3), F_b, g2

c===  Write peak dur'n st'd candle cosmo model params.
            else if (fornum(i) .eq. 10) then
                write(lun,120) (params(j,i),j=1,7), 1./(1.+params(7,i))

c===  Write beamed cosmo model params.
            else if (fornum(i) .eq. 11) then
                write(lun,130) (params(j,i),j=1,8)

c===  File spec model has no params!
            else if (fornum(i) .eq. 12) then
                continue

c===  Write rat'l broken power law (theta) model params.
            else if (fornum(i) .eq. 13) then
                write(lun,140) params(1,i), tan(params(2,i)),
     *             params(2,i), params(3,i)

            endif

c===  Go get next ingredient.
200     continue

20	format('gamma   ',1pg12.4)
30	format('g1.g2   ',2(1pg12.4),/,
     .         'Fb      ',1pg12.4)
40	format('g1.g2   ',2(1pg12.4),/,
     .         'Fb.l    ',2(1pg12.4))
50	format('spectrum   ',3(1pg12.4),/,
     .         'Omega    ',1pg13.4,/,
     .         'n0       ',1pg13.4,/,
     .         'beta     ',1pg13.4,/,
     .         'nu       ',1pg13.4)
60	format('spectrum   ',3(1pg12.4),/,
     .         'Omega      ',1pg13.4,/,
     .         'n0         ',1pg13.4,/,
     .         'beta.zcut  ',2(1pg13.4),/,
     .         'p.nu_u.rho ',3(1pg13.4))
70	format('amp     ',1pg13.4,/,
     .         'e.rhoc  ',2(1pg13.4),/,
     .         'nu      ',1pg13.4,/,
     .         'f_A     ',1pg13.4)
80	format('n0.Rc.fA ',3(1pg13.4),/,
     .         'p.nu_u.rho ',3(1pg13.4))
90	format('spectrum     ',3(1pg13.4),/,
     .         'Omega.nu_c ',2(1pg13.4),/,
     .         'rhoc.nu_h  ',2(1pg13.4),/,
     .         'f_h     ',1pg13.4)
100	format('g1     ',1pg12.4,/,
     .         'sigma  ',1pg12.4,/,
     .         'tau_0  ',1pg12.4,/,
     .         'F_tau  ',1pg12.4,/,
     .         'g2     ',1pg12.4)
120	format('spectrum   ',3(1pg12.4),/,
     .         'Omega    ',1pg13.4,/,
     .         'beta     ',1pg13.4,/,
     .         'nu       ',1pg13.4,/,
     .         'zcrit    ',1pg13.4,/,
     .         'tfac     ',1pg13.4)
130	format('spectrum   ',3(1pg12.4),/,
     .         'Omega    ',1pg13.4,/,
     .         'beta     ',1pg13.4,/,
     .         'nu       ',1pg13.4,/,
     .         'D_l      ',1pg13.4,/,
     .         'tau      ',1pg13.4)
140	format('g1.g2   ',2(1pg12.4),/,
     .         'theta   ',1pg12.4,/,
     .         'Fb      ',1pg12.4)

180	format('R.eff   ',1pg12.4)

c---  Write the efficiency type and various global info.
	write(lun,'(a,i5)') 'eta.type  ',neta
        if (didl) then
            write(lun,'(a,a)') 'dfile  ',dfile(1:nwlen(dfile))
            if (dslxn .eq. 0) then
                write(lun,'(a)') 'dslxn  All'
            else if (dslxn .eq. 1) then
                write(lun,'(a)') 'dslxn  T90'
            else if (dslxn .eq. 2) then
                write(lun,'(a)') 'dslxn  T50'
            endif
            write(lun,'(a,2i5)') 'nb.nuse ',nb,nuse
            write(lun,'(a,1pg14.6)') 'llike  ',llcur
        endif
        write(lun,*)

c---  Now write out the F grid, including the cumulative dist'n,
c---  the logarithmic slope, and eta_F.
	cum = 1.
	do 400 i=nF1, nFg
	    F = Fgrid(i)
	    if (i .eq. nF1) then
	        Fn = Fgrid(i+1)
	        c1 = F
	        c2 = Fn
	        l1 = log(max(dR(i),tiny))
	        l2 = log(max(dR(i+1),tiny))
	        cl1 = log(cum)
	        cl2 = log(cum - R(i)/Rtot)
	    else if (i .eq. nFg) then
	        Fp = Fgrid(i-1)
	        c1 = Fp
	        c2 = F
	        l1 = log(max(dR(i-1),tiny))
	        l2 = log(max(dR(i),tiny))
	        arg1 = cum + R(i-1)/Rtot
	        if (arg1 .gt. 0. .and. cum .gt. 0.) then
	            cl1 = log(arg1)
	            cl2 = log(cum)
	        endif
	    else
	        Fp = Fgrid(i-1)
	        Fn = Fgrid(i+1)
	        c1 = Fp
	        c2 = Fn
	        l1 = log(max(dR(i-1),tiny))
	        l2 = log(max(dR(i+1),tiny))
	        arg1 = cum + R(i-1)/Rtot
	        arg2 = cum - R(i)/Rtot
	        if (arg1 .gt. 0. .and. arg2 .gt. 0.) then
	            cl1 = log(arg1)
	            cl2 = log(arg2)
	        endif
	    endif
	    g = - (l2 - l1) / log(c2/c1)
	    cg = - (cl2 - cl1) / log(c2/c1)
	    write(lun,'(6(1pg12.4))') F,dR(i),g,cum,cg,eta_F(F)
	    if (i .lt. nFg) cum = cum - R(i)/Rtot
400	continue

c---  If this is the SSC model, write out the separate parts, too,
c---  and the Bayes factor and odds for each burst.
	write(lun,*)
	write(lun,420) R_eff(1), R_eff(2)
420	format('Reff.c.h  ',2(1pg14.6))
	do 500 i=1, nFg
	    write(lun,'(a,4(1pg12.4))') 'iso  ',Fgrid(i),
     *         dRdF(i,1), dRdF(i,2), dR(i)
500	continue

	write(lun,*)
	do 600 i=1, nb
	    if (useb(i)) then
	       B = (dR_h(i)/R_eff(2)) / (dR_c(i)/R_eff(1))
	       write(lun,620) tnums(i), B, log10(B),
     *              fac*B, log10(fac*B)
	    endif
600	continue
620	format('burst  ',i5,4(1pg12.4))

c---  That's it!
	close(lun)

        return
        end

c-----------------------------------------------------------------------
        subroutine wdrate (lun, name, Flo, Fhi, nstp)

c---	"Write Differential RATE"
c
c---    Write out the current differential rate.

c+++  Arguments:
	integer lun
	character*(*) name
	real*8 Flo, Fhi
	integer nstp

c+++  Globals:
        include 'modchar.cb'
        include 'data.par'
        include 'modsyms.cb'
        include 'data.cb'

c+++  Locals:
	real*8 tiny, zero, year
	parameter (zero = 0., tiny = 1.e-30, year = 3.1557d7)
        integer nF, i, j
        real*8 dR1, dR2, F1, F2, ls, del, dF, s
        real*8 fac
        logical logstp

c+++  Functions:
        real*8 dR_dFdO

c---  Make sure the current model params have been set.
	call setmod

c---  If this is the file spec model, quit!
	if (curmod .eq. 12) then
	    call quit('WDRATE:  Cannot do file spec models!')
	endif

c---  Right now, we are not setup for superposed models.
	if (curtype .ne. 1 .and. curtype .ne. 2) 
     *      call quit('WDRATE:  Bad model type!')

c---  Open the file.
	call openf(lun, name)

c---  First write a header identifying the model elements.
c---  Also, set the scale factor that ensures the rate is
c---  bursts per year per unit flux.
	write(lun,'(a,i3)') 'elements  ',nform
	if (curmod .eq. 8) then
	    if (f_h .ne. 0. .and. f_h .ne. 1.) then
	        fac = f_h / (1. - f_h)
	    endif
	endif
	s = 1.
        do 200 i=1, nform
            write(lun,'(a)') fornam(fornum(i))

c===  Write power law model params.
            if (fornum(i) .eq. 1) then
                write(lun,20) (params(j,i),j=1,1)

c===  Write rat'l broken power law model params.
            else if (fornum(i) .eq. 2) then
                write(lun,30) (params(j,i),j=1,3)

c===  Write broken power law model params.
            else if (fornum(i) .eq. 3) then
                write(lun,40) (params(j,i),j=1,4)

c===  Write st'd candle cosmo model params.
            else if (fornum(i) .eq. 4) then
                write(lun,50) (params(j,i),j=1,7)
                s = year

c===  Write power-law luminosity func'n cosmo model params.
            else if (fornum(i) .eq. 5) then
                write(lun,60) (params(j,i),j=1,10)
                s = year

c===  Write st'd candle halo model params.
            else if (fornum(i) .eq. 6) then
                write(lun,70) (params(j,i),j=1,5)

c===  Write power-law luminosity func'n halo model params.
            else if (fornum(i) .eq. 7) then
                write(lun,80) (params(j,i),j=1,6)

c===  Write superposed st'd candle model params.
            else if (fornum(i) .eq. 8) then
                write(lun,90) (params(j,i),j=1,8)
                s = year

c===  Write peak duration broken power law model params.
            else if (fornum(i) .eq. 9) then
                write(lun,100) (params(j,i),j=1,3), F_b, g2

c===  Write peak dur'n st'd candle cosmo model params.
            else if (fornum(i) .eq. 10) then
                write(lun,120) (params(j,i),j=1,7), 1./(1.+params(7,i))
                s = year

c===  Write beamed cosmo model params.
            else if (fornum(i) .eq. 11) then
                write(lun,130) (params(j,i),j=1,8)
                s = year

c===  File spec model has no params!
            else if (fornum(i) .eq. 12) then
                continue

c===  Write rat'l broken power law (theta) model params.
            else if (fornum(i) .eq. 13) then
                write(lun,140) params(1,i), tan(params(2,i)),
     *             params(2,i), params(3,i)

            endif

c===  Go get next ingredient.
200     continue

20	format('gamma   ',1pg12.4)
30	format('g1.g2   ',2(1pg12.4),/,
     .         'Fb      ',1pg12.4)
40	format('g1.g2   ',2(1pg12.4),/,
     .         'Fb.l    ',2(1pg12.4))
50	format('spectrum   ',3(1pg12.4),/,
     .         'Omega    ',1pg13.4,/,
     .         'n0       ',1pg13.4,/,
     .         'beta     ',1pg13.4,/,
     .         'nu       ',1pg13.4)
60	format('spectrum   ',3(1pg12.4),/,
     .         'Omega      ',1pg13.4,/,
     .         'n0         ',1pg13.4,/,
     .         'beta.zcut  ',2(1pg13.4),/,
     .         'p.nu_u.rho ',3(1pg13.4))
70	format('amp     ',1pg13.4,/,
     .         'e.rhoc  ',2(1pg13.4),/,
     .         'nu      ',1pg13.4,/,
     .         'f_A     ',1pg13.4)
80	format('n0.Rc.fA ',3(1pg13.4),/,
     .         'p.nu_u.rho ',3(1pg13.4))
90	format('spectrum     ',3(1pg13.4),/,
     .         'Omega.nu_c ',2(1pg13.4),/,
     .         'rhoc.nu_h  ',2(1pg13.4),/,
     .         'f_h     ',1pg13.4)
100	format('g1     ',1pg12.4,/,
     .         'sigma  ',1pg12.4,/,
     .         'tau_0  ',1pg12.4,/,
     .         'F_tau  ',1pg12.4,/,
     .         'g2     ',1pg12.4)
120	format('spectrum   ',3(1pg12.4),/,
     .         'Omega    ',1pg13.4,/,
     .         'beta     ',1pg13.4,/,
     .         'nu       ',1pg13.4,/,
     .         'zcrit    ',1pg13.4,/,
     .         'tfac     ',1pg13.4)
130	format('spectrum   ',3(1pg12.4),/,
     .         'Omega    ',1pg13.4,/,
     .         'beta     ',1pg13.4,/,
     .         'nu       ',1pg13.4,/,
     .         'D_l      ',1pg13.4,/,
     .         'tau      ',1pg13.4)
140	format('g1.g2   ',2(1pg12.4),/,
     .         'theta   ',1pg12.4,/,
     .         'Fb      ',1pg12.4)


        write(lun,*)

c---  Now write out the rate.
	call setstp (Flo, Fhi, nstp, logstp, del)
	do 300 nF=1, abs(nstp)
	    call stpinf (nF, Flo, del, logstp, F1, dF)
	    F2 = F1 + 0.1*dF
	    if (curtype .eq. 1) then
	        dR1 = pi4 * dR_dFdO(F1,zero,zero)
	        dR2 = pi4 * dR_dFdO(F2,zero,zero)
	    else if (curtype .eq. 2) then
	        call sky_int0(F1, dR1)
	        call sky_int0(F2, dR2)
	    endif
	    ls = -(log(max(dR2,tiny)) - log(max(dR1,tiny))) / log(F2/F1)
	    write(lun,'(3(1pg12.4))') F1, s*dR1, ls
300	continue
	close(lun)
	
	return
	end
