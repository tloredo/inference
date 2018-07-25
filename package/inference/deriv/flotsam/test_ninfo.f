	program test_info

c	Test info matrix routines.
c
c	6 May 97:  TJL (adapted from earlier codes)

        integer npar
        parameter (npar = 4)
        real*8 mode(npar)
        real*8 dp(npar), info(npar,npar), ierr(npar,npar), pmax
        integer i, j
        real*8 mu(npar), h(npar,npar)
        common /parcb/ h, mu

	real*8 func
	external func

c---  Define the q-form.
        mu(1) = 1.0
        mu(2) = 3.0
        mu(3) = -2.5
        mu(4) = 2.
        h(1,1) = 1.
        h(2,1) = -.5
        h(2,2) = 1.
        h(3,1) = -.2
        h(3,2) = .5
        h(3,3) = 1.
        h(4,1) = .3
        h(4,2) = -.2
        h(4,3) = .1
        h(4,4) = 1.
        do 10 i=1, npar-1
            do 5 j=i+1, npar
                h(i,j) = h(j,i)
5	    continue
10	continue
        do 13 i=1, npar
            do 12 j=1, npar
                h(i,j) = 0.1*h(i,j)
12	    continue
13	continue

c---  Set steps.
        mode(1) = 1.
        mode(2) = 3.
        mode(3) = -2.5
        mode(4) = 2.

        dp(1) = 0.01
        dp(2) = 0.01
        dp(3) = 0.01
        dp(4) = 0.01
        pmax = func(mode)

        call ninfo(func, npar, mode, pmax, dp, 5, info, ierr, npar)
        write(*,*) 'Hess from ninfo:'
        write(*,'(4(1pg12.4))') ((info(i,j), j=1, npar), i=1, npar)
        write(*,*)
        write(*,*) 'Errors:'
        write(*,'(4(1pg12.4))') ((ierr(i,j), j=1, npar), i=1, npar)
        write(*,*)

        call ninfo(func, npar, mode, pmax, dp, 1, info, ierr, npar)
        write(*,*) 'Hess from ninfo, 1 iter:'
        write(*,'(4(1pg12.4))') ((info(i,j), j=1, npar), i=1, npar)
        write(*,*)
        write(*,*) 'Errors:'
        write(*,'(4(1pg12.4))') ((ierr(i,j), j=1, npar), i=1, npar)
        write(*,*)

        call ninfol(func, npar, mode, pmax, dp, 5, info, ierr, npar)
        write(*,*) 'Hess from ninfol:'
        write(*,'(4(1pg12.4))') ((info(i,j), j=1, npar), i=1, npar)
        write(*,*)
        write(*,*) 'Errors:'
        write(*,'(4(1pg12.4))') ((ierr(i,j), j=1, npar), i=1, npar)
        write(*,*)

        call ninfol(func, npar, mode, pmax, dp, 1, info, ierr, npar)
        write(*,*) 'Hess from ninfol, 1 iter:'
        write(*,'(4(1pg12.4))') ((info(i,j), j=1, npar), i=1, npar)
        write(*,*)
        write(*,*) 'Errors:'
        write(*,'(4(1pg12.4))') ((ierr(i,j), j=1, npar), i=1, npar)
        write(*,*)

        end

c-----------------------------------------------------------------------
        real*8 function func(x)

        integer npar
        parameter (npar = 4)
        real*8 x(*), q
        integer i, j
        real*8 mu(npar), h(npar,npar)
        common /parcb/ h, mu

        q = 0.
        do 20 i=1, npar
            do 10 j=1, i
                q = q + (x(i)-mu(i))*h(i,j)*(x(j)-mu(j))
                if (i .ne. j) q = q + (x(i)-mu(i))*h(i,j)*(x(j)-mu(j))
10          continue
20      continue
c        q = q + q*q

        func = exp(-q/2.)

        return
        end

c-----------------------------------------------------------------------
       subroutine quit(char)

        character*(*) char

        write(*,'(a)') char
        stop

        end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c	Routines for calculating numerical derivatives.
c
c        subroutine ninfo(f, N, mu, fmax, dp, iter, info, ierr, ndim)
c        subroutine ninfol(lf, N, mu, fmax, dp, iter, info, ierr, ndim)
c        subroutine rid_ex(h, rat, y, n, result, err)
c
c	nhess, ninfo, ninfol, and rid_ex are from an older collection.
c
c	Created 19 May 97 by Tom Loredo
c	
c-----------------------------------------------------------------------
	subroutine pk_scales (n, f, x, fmax, df, dx, ndim)

c---	"PeaK SCALES"
c
c	Find the scales dx(i) for variables x(i) that lead to a change
c	of df in the function f(x) about its current location, presumed
c	to be a local mode.
c
c	Note that the convergence criterion is a difference, not
c	a ratio, and is thus most appropriate when f is the logarithm
c	of a quantity of interest.
c
c	n	-->	Number of variables.
c	f	-->	The function name.
c	x	-->	The location of the mode.
c	fmax	-->	The value of f at the mode (to save a possibly
c			redundant function call).
c	df	-->	The change in f setting the scales.
c	dx	<->	Initially a guess at the scales; on return
c			the estimated scales.
c	ndim	-->	Physical dimension for x, dx.

c+++  Arguments:
	integer n, ndim
	real*8 f, x(ndim), fmax, df, dx(ndim)
	external f

c+++  Locals:
	integer maxnr
	real*8 tol
	parameter (maxnr = 20, tol = 0.01)
	integer i, nit
	real*8 f_c, s, x0, fnew, del, fp, dfdx

c---  Set the critical value.
	f_c = fmax - df
c	print *, 'fvals: ', fmax, f_c
c	write(9,'(a,3(1pg12.4))') 'fvals: ',x(1), fmax, f_c

c---  Loop thru the parameters.
        do 40 i=1, n

c===  For first guess, use a quadratic approximation.
            s = 1.
            x0 = x(i)
10          x(i) = x0 + s*dx(i)
            fnew = f(x)
            if (fnew .gt. fmax) call quit('Offset from mode!')
            del = fmax - fnew
            if (del .le. 0.) then
                s = 1.3 * s
                go to 10
            endif
            x(i) = x0 + s*dx(i) * sqrt(df/del)
c            print *, 'guess: ', i, x0, x(i)

c===  Now do the NR iterations.
            nit = 0
20          fnew = f(x)
c	print *, nit, x(i), fnew
c	write(9,'(i5,2(1pg12.4))') nit, x(i), fnew
            if (fnew .gt. fmax) call quit('Offset from mode!')
            nit = nit + 1
            if (nit .gt. maxnr)
     *          call quit('pk_scales:  Too many NR iters!')
            if (abs(fnew-f_c) .gt. tol) then
                x(i) = x(i) + s*dx(i)
                fp = f(x)
                if (fp .gt. fmax) call quit('Offset from mode!')
                dfdx = (fp - fnew) / (s*dx(i))
                x(i) = x(i) - s*dx(i) - (fnew-f_c)/dfdx
                go to 20
            endif
            dx(i) = x(i) - x0
            x(i) = x0
40      continue

	return
	end

c-----------------------------------------------------------------------
        subroutine ninfo(f, N, mu, fmax, dp, iter, info, ierr, ndim)

c---    Numerically calculate the information matrix of
c       the function f (of N parameters) about the point x = mu.
c       It is assumed that x = mu is the maximum of f, so all first
c       derivatives vanish there.
c
c       The information matrix = - (d/dx_i)(d/dx_j) log f, i.e.,
c       the negative Hessian, divided by f at the mode.
c
c       Input:
c
c           f       A real*8 function of a vector argument of length N.
c
c           N       The number of parameters of f.  (max = 15)
c
c           mu(ndim)    The parameters at the maximum of f.
c
c           fmax    The value of f at mu, passed as an argument
c                   to save a possibly redundant call to f.
c
c           dp(ndim)    An array of step sizes for the N parameters,
c                       used to calculate numerical derivatives of f.
c
c           iter    The number of iterations to use in the
c                   Richardson extrapolation.  Each interation
c                   kills off one more even power in the order of
c                   the error.
c
c           ndim    The dimension of the above arrays.
c
c       Returned:
c
c           info(ndim,ndim) The NxN negative Hessian matrix,
c                           divided by the function value at x = mu.
c
c	    ierr(ndim,ndim) Estimated errors in each term of info.
c			    If iter=1, ierr=0.
c
c       Algorithm:
c
c       Derivatives are calculated by iterated central differencing
c       with Richardson extrapolation.  The diagonal terms are
c       calculated first, and the resulting "sigmas" used to rescale
c       the parameters for calculation of off-diagonal terms.
c
c       Central differencing with step size h gives f'' to order h**2.
c       One Richardson step subtracts off the h**2 term, leaving
c       an error of order h**4, etc.  The Richardson extrapolation
c       is equivalent to finding N central differences with different
c       step sizes, and extrapolating to zero step size with a
c       Laplace expansion of degree N.  This is essentially
c       "Romberg differentiation".
c
c       This method is due to R.E. Kass, "Computing Observed Information
c       by Finite Differences", Commun. Stat. Sim. v. 16, 587-599 
c	(1987), with some errors in that paper here corrected, and the 
c       equivalence with Laplace/Romberg here identified.

        integer npmax, itmax
        real*8 ratio, one
        parameter (npmax = 15, itmax = 10, ratio = 0.707, one = 1.d0)

	integer ndim
        real*8 f, fmax, lfmax, mu1, mu2
        real*8 mu(ndim), dp(ndim), info(ndim,ndim), ierr(ndim,ndim)
        integer N, iter, i, j, k
        real*8 h(itmax), d2f(itmax), di

c---  Check N and iter.
        if (N .gt. npmax) pause 'ninfo:  N too large!'
        if (iter .gt. itmax) pause 'ninfo:  iter too large!'
        lfmax = log(fmax)

c---  First find the derivatives along the diagonal (unmixed).
        do 100 i=1, N
            mu1 = mu(i)

c>>>  Begin by finding iter second differences.
            do 20 j=1, iter
                h(j) = dp(i) * ratio**(iter-j)
                mu(i) = mu1 + h(j)
                d2f(j) = log(f(mu))
                mu(i) = mu1 - h(j)
                d2f(j) = (d2f(j) + log(f(mu)) - 2.*lfmax)/h(j)**2
20          continue

c>>>  Now extrapolate them to h = 0.
	    if (iter .gt. 1) then
                call rid_ex(dp(i), ratio, d2f, iter, info(i,i), di)
                info(i,i) = -info(i,i)
                ierr(i,i) = abs(di)
            else
                info(i,i) = -d2f(1)
                ierr(i,i) = 0.
            endif

c>>>  End loop over the diagonal.
            mu(i) = mu1
100     continue

c---  Now find the mixed derivatives
        do 200 i=1, N
            mu1 = mu(i)
            do 180 j=i+1, N
                mu2 = mu(j)

c>>>  Find iter second differences of f in terms of a
c>>>  new "diagonal" variable.
                do 120 k=1, iter
                    h(k) = ratio**(iter-k)
                    mu(i) = mu1 + h(k)/sqrt(info(i,i))
                    mu(j) = mu2 + h(k)/sqrt(info(j,j))
                    d2f(k) = log(f(mu))
                    mu(i) = mu1 - h(k)/sqrt(info(i,i))
                    mu(j) = mu2 - h(k)/sqrt(info(j,j))
                    d2f(k) = (d2f(k) + log(f(mu)) - 2.*lfmax)/h(k)**2
120             continue

c>>>  Now extrapolate them to h = 0, and eliminate the
c>>>  effects of the transformation.
		if (iter .gt. 1) then
                    call rid_ex(one, ratio, d2f, iter, info(i,j), di)
                    info(i,j) = -(info(i,j) + 2.) * 
     *                           sqrt(info(i,i)*info(j,j)) / 2.
	            ierr(i,j) = abs(di)*sqrt(info(i,i)*info(j,j)) / 2.
                    info(j,i) = info(i,j)
                    ierr(j,i) = ierr(i,j)
                else
                    info(i,j) = -(d2f(1) + 2.) * 
     *                           sqrt(info(i,i)*info(j,j)) / 2.
	            ierr(i,j) = 0.
                    info(j,i) = info(i,j)
                    ierr(j,i) = 0.
                endif

c>>>  End the loops.
                mu(j) = mu2
180         continue
            mu(i) = mu1
200     continue

        return
        end

c-----------------------------------------------------------------------
        subroutine ninfol(lf, N, mu, lfmax, dp, iter, info, ierr, ndim)

c---	"Numerical INFOrmation matrix from Logarithm"
c
c	This is identical to ninfo, above, except that it uses
c	the *logarithm* of the function, rather than the function
c	itself.  This saves a few calls to log(), but more
c	importantly, it works for functions whose logs are easily
c	calculated and too large (positive or negative) for
c	straightforward exponentiation.

c+++  Arguments:
	integer ndim
	real*8 lf, mu(ndim), lfmax, dp(ndim)
	real*8 info(ndim,ndim), ierr(ndim,ndim)
	integer N, iter

c+++  Locals:
        integer npmax, itmax
        real*8 ratio, one
        parameter (npmax = 15, itmax = 10, ratio = 0.707, one = 1.d0)

        real*8 mu1, mu2
        integer i, j, k
        real*8 h(itmax), d2f(itmax), di

c---  Check N and iter.
        if (N .gt. npmax) pause 'ninfo:  N too large!'
        if (iter .gt. itmax) pause 'ninfo:  iter too large!'

c---  First find the derivatives along the diagonal (unmixed).
        do 100 i=1, N
            mu1 = mu(i)

c>>>  Begin by finding iter second differences.
            do 20 j=1, iter
                h(j) = dp(i) * ratio**(iter-j)
                mu(i) = mu1 + h(j)
                d2f(j) = lf(mu)
                mu(i) = mu1 - h(j)
                d2f(j) = (d2f(j) + lf(mu) - 2.*lfmax)/h(j)**2
20          continue

c>>>  Now extrapolate them to h = 0.
	    if (iter .gt. 1) then
                call rid_ex(dp(i), ratio, d2f, iter, info(i,i), di)
                info(i,i) = -info(i,i)
                ierr(i,i) = abs(di)
            else
                info(i,i) = -d2f(1)
                ierr(i,i) = 0.
            endif

c>>>  End loop over the diagonal.
            mu(i) = mu1
100     continue

c---  Now find the mixed derivatives
        do 200 i=1, N
            mu1 = mu(i)
            do 180 j=i+1, N
                mu2 = mu(j)

c>>>  Find iter second differences of f in terms of a
c>>>  new "diagonal" variable.
                do 120 k=1, iter
                    h(k) = ratio**(iter-k)
                    mu(i) = mu1 + h(k)/sqrt(info(i,i))
                    mu(j) = mu2 + h(k)/sqrt(info(j,j))
                    d2f(k) = lf(mu)
                    mu(i) = mu1 - h(k)/sqrt(info(i,i))
                    mu(j) = mu2 - h(k)/sqrt(info(j,j))
                    d2f(k) = (d2f(k) + lf(mu) - 2.*lfmax)/h(k)**2
120             continue

c>>>  Now extrapolate them to h = 0, and eliminate the
c>>>  effects of the transformation.
	        if (iter .gt. 1) then
                    call rid_ex(one, ratio, d2f, iter, info(i,j), di)
                    info(i,j) = -(info(i,j) + 2.) * 
     *                           sqrt(info(i,i)*info(j,j)) / 2.
	            ierr(i,j) = abs(di)*sqrt(info(i,i)*info(j,j)) / 2.
                    info(j,i) = info(i,j)
                    ierr(j,i) = ierr(i,j)
                else
                    info(i,j) = -(d2f(1) + 2.) * 
     *                           sqrt(info(i,i)*info(j,j)) / 2.
	            ierr(i,j) = 0.
                    info(j,i) = info(i,j)
                    ierr(j,i) = 0.
                endif

c>>>  End the loops.
                mu(j) = mu2
180         continue
            mu(i) = mu1
200     continue

        return
        end

c-----------------------------------------------------------------------
        subroutine rid_ex(h, rat, y, n, result, err)

c---	"RIDder's EXtrapolation"
c
c       Extrapolate the derivative estimates in y to smaller
c	stepsize and/or higher order.  
c
c       Adapted from Numerical Recipes' DFRIDR routine.

        integer nmax
        real*8 big
        parameter (nmax = 10, big = 1.e30)

        real*8 h, rat, y(*), result, err
        real*8 con, con2, errt, fac, hh, a(nmax,nmax)
        integer i, j, n

        if (n .gt. nmax) pause 'zextrap2:  n too big!'
        con = 1. / rat
        con2 = con * con
        hh = h
        a(1,1) = y(n)
        err = big

        do 60 i=2, n
            hh = hh / con
            a(1,i) = y(n-i+1)
            fac = con2
            do 40 j=2, i
                a(j,i) = (a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
                fac = con2 * fac
                errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                if (errt .le. err) then
                    err = errt
                    result = a(j,i)
                endif
40          continue
	    if (abs(a(i,i)-a(i-1,i-1)) .ge. 2.*err) return
60      continue

        return
        end

