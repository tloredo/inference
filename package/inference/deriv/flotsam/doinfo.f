c-----------------------------------------------------------------------
        subroutine doinfo (del, nrstp, ok, posev, ldet, lockham, cnum, 
     *     flun, ofile)

c---    "DO INFOrmation matrix calculation"
c
c       Calculate the information matrix of the posterior probability
c       about the current parameter variables, assumed to be the mode.
c	Also calculate its determinant and inverse, the Ockham
c	factor (up to the prior factor), and the condition number.
c
c       del is the fractional change in the posterior density used to 
c       determine the parameter changes used to calculate derivatives.

c+++  Argument:
	real*8 del, ldet, lockham, cnum
	integer nrstp, flun
	character*(*) ofile
	logical ok, posev

c+++  Global variables:
	include 'info.cb'

c+++  Local variables:
        integer npar
        real*8 twopi
        parameter (twopi = 2.d0*3.141592653589793d0)
        real*8 L(npmax,npmax), vec(npmax), w1(npmax), w2(npmax)
        real*8 nfunc
        external nfunc

        real*8 frac, tol, f, diff, mfrac
        integer maxnr, i, j, nit
        real*8 lpp_tv, dp(npmax), xval(npmax), mode(npmax), tval(npmax)
        real*8 lpjmax, lpjc, lpjnew, dlpj, t0, lpp
        external lpp_tv
        real*8 vlims(npmax,2)
        character*9 vsyms(npmax)
        parameter (mfrac = 0.01, tol = 0.03, maxnr = 15)
c        real*8 copy(npmax,npmax)
c        integer k

        if (del .eq. 0.) call quit('Bad scale for info calc!')

c---  Use current values of varying symbols for the mode.  Convert it
c---  to transformed coordinates.
        call varinf(.false., npar, vsyms, xval, vlims, dp, npmax)
        if (npar .gt. npmax) call quit('Too many params for I matrix!')
        call xttran(xval, mode, 'f')

c---  First set scales for the derivative calculations.
c---  For each variable, find the change in the variable value
c---  so that log(pp_tv) = log(pjmax) - del.  Use Newton-Raphson.
        lpjmax = lpp_tv(mode)
        lpjc = lpjmax - del
        call xttran(xval, tval, 'f')
        do 40 i=1, npar


c===  For first guess, use a quadratic approximation.
            f = 1.
            t0 = tval(i)
            if (t0 .gt. 0.) then
                frac = mfrac * t0
            else
                frac = mfrac
            endif
10          tval(i) = t0 + f*frac
            lpjnew = lpp_tv(tval)
            if (lpjnew .gt. lpjmax) call quit('Offset from mode!')
            diff = lpjmax - lpjnew
            if (diff .le. 0.) then
                f = 1.3 * f
                go to 10
            endif
            tval(i) = t0 + f*frac * sqrt(del/diff)

c===  Now do the NR iterations.
            nit = 0
20          lpjnew = lpp_tv(tval)
            if (lpjnew .gt. lpjmax) call quit('Offset from mode!')
            nit = nit + 1
            if (nit .gt. maxnr) call quit('doinfo:  Too many NR iters!')
            if (abs(lpjnew-lpjc) .gt. tol) then
                tval(i) = tval(i) + f*frac
                lpp = lpp_tv(tval)
                if (lpp .gt. lpjmax) call quit('Offset from mode!')
                dlpj = (lpp - lpjnew) / (f*frac)
                tval(i) = tval(i) - f*frac - (lpjnew-lpjc)/dlpj
                go to 20
            endif
            dp(i) = tval(i) - t0
            tval(i) = t0
40      continue

c---  Now actually calculate the Hessian.  Use nrstp Richardson
c---  extrapolation steps.
        call ninfol(lpp_tv, npar, mode, lpjmax, dp, nrstp, 
     *               info, ierr, npmax)

c---  Find its Cholesky decomposition and determinant. 
	call choldet(npar, info, L, npmax, ldet, ok)

c---  Find the eigens.  This call destroys the superdiagonal part
c---  of info, so restore it.
	call eigens(info, evecs, evals, order, npar, npmax, w1, w2)
	do 70 i=1, npar-1
	    do 60 j=i+1, npar
	        info(i,j) = info(j,i)
60	    continue
70	continue

c---  Find the determinant of V, and the Ockham factor.  If choldet
c---  worked, this is easy; otherwise, use eigens.
	if (ok) then
	    ldet = - ldet
	else
	    ldet = 0.
	    posev = .true.
	    do 75 i=1, npar
	        if (evals(i) .le. 0.) then
	            posev = .false.
	            goto 75
	        endif
	        ldet = ldet - log(evals(i))
75	    continue
	endif
	lockham = 0.5*npar*log(twopi) + 0.5*ldet

c---  Write out the mode, Hessian, & determinant.
        call openf(flun, ofile)
        write(flun,'(a,i3,/)')     '#.params ',npar
        do 80 i=1, npar
            write(flun,'(a,i3,3x,a)') 'param    ',i,vsyms(i)
80	continue
	write(flun,*)
        if (ok) then
            write(flun,'(a,1pg12.4)') 'log.det.V   ',ldet
            write(flun,'(a,1pg12.4)') 'log.Ockham  ',lockham
        else
            write(flun,'(a)') 'singular'
        endif
        write(flun,*)
        call fwrows('mode  ',xval,npar,5,flun,'(1pg12.4)')
        call fwrows('tmode ',xval,npar,5,flun,'(1pg12.4)')
        write(flun,*)
        write(flun,*)
        do 100 i=1, npar
            do 90 j=1, npar
                vec(j) = info(i,j)
90	    continue
	    call fwrows ('I.row ', vec, npar, 5, flun, '(1pg12.4)')
	    write(flun,*)
100	continue
	write(flun,*)
        do 110 i=1, npar
            do 105 j=1, npar
                vec(j) = ierr(i,j)
105	    continue
	    call fwrows ('I.err ', vec, npar, 5, flun, '(1pg12.4)')
	    write(flun,*)
110	continue
	write(flun,*)

c...  Copy info to test its inverse later.
c	do 105 i=1, npar
c	  do 105 j=1, npar
c	     copy(i,j) = info(i,j)
c105	continue

c---  Write out eigens.
	if (abs(evals(order(npar))) .gt. 0.) then
	    cnum = abs( evals(order(1)) / evals(order(npar)) )
	else
	    cnum = -1.
	endif
	do 120 i=1, npar
	    vec(i) = evals(order(i))
120	continue
	call fwrows('eigens  ', vec, npar, 5, flun, '(1pg12.4)')
	write(flun,*)
	do 160 i=1, npar
	    do 140 j=1, npar
	        vec(j) = evecs(i,order(j))
140	    continue
	    call fwrows ('O.row ', vec, npar, 5, flun, '(1pg12.4)')
	    write(flun,*)
160	continue
	write(flun,*)

c---  Calculate the inverse & write it out, if choldet worked.
	if (ok) then
	call cholinv(npar, L, info, npmax)        
	    do 240 i=1, npar
                do 220 j=1, npar
                    vec(j) = info(i,j)
220	        continue
	        call fwrows ('V.row ', vec, npar, 5, flun, '(1pg12.4)')
	        write(flun,*)
240	    continue
	endif
	close(flun)

c...  Test the inverse.  Multiply by the original, & check that
c...  its determinant is the inverse of the original.
c	do 200 i=1, npar
c	  do 190 j=1, npar
c	    f = 0.
c	    do 180 k=1, npar
c	       f = f + info(i,k)*copy(k,j)
c180	    continue
c	    print *, i, j, f
c190	  continue
c200	continue
c	call choldet(npar, info, L, npmax, ldet)
c	print *, 'invdet = ',ldet
	    
        return
        end

c-----------------------------------------------------------------------
	subroutine choldet (n, matrix, L, ndim, ldet, ok)

c---	"CHOLesky decomposition and DETerminant"
c
c	The Cholesky decomposition of the symmetric positive-definite 
c	n x n MATRIX is returned as the lower triangular matrix L; the
c	log determinant is returned as ldet.  The physical
c	dimensions of MATRIX and L are ndim x ndim.
c
c       The algorithm is from LINPACK, as published (in Algol) in
c       Num. Mat. 7, p. 362, 1965 (choldet1).  It is also described
c       (more clearly?) in R. Y. Rubinstein's book, "Simulation and
c       the Monte Carlo Method," p. 68.  The method of storage for 
c       L has been altered from choldet1.
c
c       The determinant code allows calculation of larger determinants
c       by writing  det = d1 * 2**d2, with d1 in the interval [1/16,1).
c	This is as in choldet1, except here we return the log det,
c	rather than d1 and d2.
c
c	For direct calculation of the determinant, comment the
c	two block-ifs after "d1 = d1 * x", and return det = d1.

c+++  Arguments:
	integer n, ndim
	real*8 matrix(ndim,ndim), L(ndim,ndim), ldet
	logical ok

c+++  Locals:
	real*8 d1, d2, x
	integer i, j, k

c---  Do the Cholesky decomposition recursively.
        d1 = 1.
        d2 = 0.
        do 100 i=1, n
            do 80 j=i, n
                x = matrix(i,j)
                if (i .gt. 1) then
                    do 35 k=i-1, 1, -1
                        x = x - L(j,k) * L(i,k)
35                  continue
                endif

c>>>  Here we get the determinant along the way as the product of
c>>>  squared diagonal elements.
                if (i .eq. j) then
                    d1 = d1 * x
                    if (x .le. 0) then
                        ldet = 0.
                        ok = .false.
                        return
                    endif
40                  if (abs(d1) .ge. 1.) then
                        d1 = d1 * 0.0625
                        d2 = d2 + 4
                        go to 40
                    endif
60                  if (abs(d1) .lt. 0.0625) then
                        d1 = d1 * 16.
                        d2 = d2 - 4
                        go to 60
                    endif
                    L(i,i) = sqrt(x)
                else
                    L(j,i) = x / L(i,i)
                    L(i,j) = 0.
                endif
80          continue
100     continue
        ldet = log(d1) + d2*log(2.)
	ok = .true.

	return
	end

c-----------------------------------------------------------------------
	subroutine cholinv (n, L, inv, ndim)

c---	"CHOLesky INVersion"
c
c	The inverse of a symmetric positive definite matrix is
c	calculated from its Cholesky decomposition, L, and returned
c	as INV.  The algorithm is from LINPACK, as described in
c	the reference cited in CHOLDET, above.

c+++  Arguments:
	integer n, ndim
	real*8 L(ndim,ndim), inv(ndim,ndim)

c+++  Locals:
	real*8 z
	integer i, j, k

c---  Copy L to the lower triangle of inv, inverting the diagonal.
	do 40 i=1, n
	    do 20 j=i, n
	        if (i .eq. j) then
	            inv(j,i) = 1. / L(j,i)
	        else
	            inv(j,i) = L(j,i)
	        endif
20	    continue
40	continue

c---  First, invert L; the inverse of L is also lower-triangular.
c---  Compute this in inv; it will be replaced by the inverse.
	do 100 i=1, n
	    do 80 j=i+1, n
	        z = 0.
	        do 60 k=j-1, i, -1
	            z = z - inv(j,k)*inv(k,i)
60	        continue
	        inv(j,i) = z * inv(j,j)
80	    continue
100	continue

c---  Then calculate the matrix inverse, which is symmetric.
	do 160 i=1, n
	    do 140 j=i, n
	        z = 0.
	        do 120 k=j, n
	            z = z + inv(k,j)*inv(k,i)
120	        continue
	        inv(j,i) = z
	        if (i .ne. j) inv(i,j) = inv(j,i)
140	    continue
160	continue

	return
	end
