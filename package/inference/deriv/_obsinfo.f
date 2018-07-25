c-----------------------------------------------------------------------
c
c       Routines for calculating numerical 2nd derivatives for finding
c       the observed Fisher information matrix from a log likelihood.
c
c       Only pk_scales and ninfol are exposed to Python.
c
c       subroutine pk_scales (n, f, x, fmax, df, dx, ndim)
c       subroutine ninfol(lf, N, mu, fmax, dp, iter, info, ierr, ndim)
c       subroutine rid_ex(h, rat, y, n, result, err)
c
c       These are adapted from an older collection.
c
c       Created 19 May 97 by Tom Loredo
c       Altered for Python 5 Jan 2006
c
c-----------------------------------------------------------------------
        subroutine pk_scales (n, f, x, fmax, df, dx, maxnr, tol, status)

c---    "PeaK SCALES"
c
c       Find the scales dx(i) for variables x(i) that lead to a change
c       of df in the function f(x) about its current location, presumed
c       to be a local mode.  f(x) will decrease by df if any
c       *individual* x(i) is incremented by dx(i).
c
c       Note that the convergence criterion is a difference, not
c       a ratio, and is thus most appropriate when f is the logarithm
c       of a quantity of interest.
c
c       n       -->     Number of variables.
c       f       -->     The function name.
c       x       <->     The location of the mode.  Note this is
c                       altered in the subroutine, but restored.
c       fmax    -->     The value of f at the mode (to save a possibly
c                       redundant function call).
c       df      -->     The change in f setting the scales.
c       dx      <->     Initially a guess at the scales; on return
c                       the estimated scales.
c       maxnr   -->     Max # of Newton-Raphson iters to solve for dx
c       tol     -->     Tolerance for achieving df
c       status  -->     1 = success
c                       2 = input x found not to be at max
c                       3 = exceeded maxnr
c
c       In typical use maxnr=20 and tol=0.01 proves adequate.

c+++  Arguments:
        integer n, maxnr, status
        real*8 f, x(n), fmax, df, dx(n), tol
        external f

cf2py intent(in) n, f, x, fmax, df
cf2py intent(in,out) dx
cf2py intent(out) status

c+++  Locals:
        integer i, nit
        real*8 f_c, s, x0, fnew, del, fp, dfdx

        status = 1

c---  Set the critical value.
        f_c = fmax - df
C        print *, 'fvals: ', fmax, f_c
c       write(9,'(a,3(1pg12.4))') 'fvals: ',x(1), fmax, f_c

c---  Loop thru the parameters.
        do 40 i=1, n

c===  For first guess, use a quadratic approximation.
            s = 1.
            x0 = x(i)
10          x(i) = x0 + s*dx(i)
            fnew = f(n, x)
c            print *, 'guess: ', i, x0, x(i), fnew
            if (fnew .gt. fmax) then
                status = 2
                return
            endif
            del = fmax - fnew
            if (del .le. 0.) then
                s = 1.3 * s
                go to 10
            endif
            x(i) = x0 + s*dx(i) * sqrt(df/del)

c===  Now do the NR iterations.
            nit = 0
20          fnew = f(n, x)
c       print *, nit, x(i), fnew
c       write(9,'(i5,2(1pg12.4))') nit, x(i), fnew
            if (fnew .gt. fmax) then
                status = 2
                return
            endif
            nit = nit + 1
            if (nit .gt. maxnr) then
                status = 3
                return
            endif
            if (abs(fnew-f_c) .gt. tol) then
                x(i) = x(i) + s*dx(i)
                fp = f(n, x)
                if (fp .gt. fmax) then
                    status = 2
                    return
                endif
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
        subroutine ninfol(lf, n, mu, lfmax, dp, iter, info, ierr,
     *                    h, d2f, a)

c---    "Numerical INFOrmation matrix from Logarithm"
c
c       Numerically calculate the information matrix of the
c       function f=exp(lf) (of N parameters) about the point x = mu.
c       It is assumed that x = mu is the maximum of f, so all first
c       derivatives vanish there.
c
c       The information matrix = - (d/dx_i)(d/dx_j) log f, i.e.,
c       the negative Hessian, divided by f at the mode.
c
c               Note that lf is not exponentiated anywhere in this algorithm,
c       so it is safe to use functions with very large or small
c       logarithms.
c
c       Input:
c
c           lf      A real*8 function of a vector argument of length N.
c
c           n       The number of parameters of f.
c
c           mu(n)    The parameters at the maximum of f.
c
c           lfmax   The value of lf at mu, passed as an argument
c                   to save a possibly redundant call to lf.
c
c           dp(n)   An array of step sizes for the N parameters,
c                   used to calculate numerical derivatives of lf.
c                   pk_scales (above) can be used to set this.
c
c           iter    The number of iterations to use in the
c                   Richardson extrapolation.  Each interation
c                   kills off one more even power in the order of
c                   the error.
c
c           h, d2f  Work arrays (real*8) of size iter.
c
c           a       Work array (real*8) of size (iter,iter).
c
c       Returned:
c
c           info(n,n) The n by n negative Hessian matrix,
c                     divided by the function value at x = mu.
c
c           ierr(n,n) Estimated errors in each term of info.
c                     If iter=1, ierr=0.
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
c       "Romberg differentiation."
c
c       This method is due to R.E. Kass, "Computing Observed Information
c       by Finite Differences", Commun. Stat. Sim. v. 16, 587-599
c       (1987), with some errors in that paper here corrected, and the
c       connection with Laplace/Romberg here identified.


c+++  Arguments:
        integer n, iter
        real*8 lf, mu(n), lfmax, dp(n)
        real*8 info(n,n), ierr(n,n)
        real*8 h(iter), d2f(iter), a(iter,iter)
        external lf

cf2py  intent(in) lf, n, mu, lfmax, dp, iter
cf2py  intent(out) info, ierr
cf2py  intent(cache,hide) h, d2f, a

c+++  Locals:
        real*8 ratio, one
        parameter (ratio = 0.707, one = 1.d0)

        real*8 mu1, mu2, sig1, sig2, di, extrap
        integer i, j, k

c---  First find the derivatives along the diagonal (unmixed).
        do 100 i=1, n
            mu1 = mu(i)

c>>>  Begin by finding iter second differences.
            do 20 j=1, iter
                h(j) = dp(i) * ratio**(iter-j)
                mu(i) = mu1 + h(j)
                d2f(j) = lf(n, mu)
                mu(i) = mu1 - h(j)
                d2f(j) = (d2f(j) + lf(n, mu) - 2.*lfmax)/h(j)**2
20          continue

c>>>  Now extrapolate them to h = 0.
            if (iter .gt. 1) then
                call rid_ex(dp(i), ratio, d2f, iter, extrap, di, a)
                info(i,i) = -extrap
                ierr(i,i) = abs(di)
            else
                info(i,i) = -d2f(1)
                ierr(i,i) = 0.
            endif

c>>>  End loop over the diagonal.
            mu(i) = mu1
100     continue

c---  Now find the mixed derivatives
        do 200 i=1, n-1
            mu1 = mu(i)
            sig1 = 1./sqrt(info(i,i))
            do 180 j=i+1, n
                mu2 = mu(j)
                sig2 = 1./sqrt(info(j,j))

c>>>  Find iter second differences of f in terms of a
c>>>  new "diagonal" variable.
                do 120 k=1, iter
                    h(k) = ratio**(iter-k)
                    mu(i) = mu1 + h(k)*sig1
                    mu(j) = mu2 + h(k)*sig2
                    d2f(k) = lf(n, mu)
                    mu(i) = mu1 - h(k)*sig1
                    mu(j) = mu2 - h(k)*sig2
                    d2f(k) = (d2f(k) + lf(n, mu) - 2.*lfmax)/h(k)**2
120             continue

c>>>  Now extrapolate them to h = 0, and eliminate the
c>>>  effects of the transformation.
                if (iter .gt. 1) then
                    call rid_ex(one, ratio, d2f, iter, extrap, di, a)
                    info(i,j) = -0.5 * (extrap + 2.) / (sig1*sig2)
                    ierr(i,j) = 0.5 *abs(di) / (sig1*sig2)
                    info(j,i) = info(i,j)
                    ierr(j,i) = ierr(i,j)
                else
                    info(i,j) = -0.5 * (d2f(1) + 2.) / (sig1*sig2)
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
        subroutine rid_ex(h, rat, y, n, result, err, a)

c---    "RIDder's EXtrapolation"
c
c       Extrapolate the derivative estimates in y to smaller
c       stepsize and/or higher order.
c
c       a is a workspace array of size (n,n)
c
c       Adapted from Numerical Recipes' DFRIDR algorithm.

        real*8 big
        parameter (big = 1.e30)

        real*8 h, rat, y(*), result, err
        real*8 con, con2, errt, fac, hh, a(n,n)
        integer i, j, n

C        if (n .gt. nmax) pause 'zextrap2:  n too big!'
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

