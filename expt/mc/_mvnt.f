c-----------------------------------------------------------------------
c
c       Multivariate random number generators and associated
c       probability density functions, for wrapping with f2py.
c
c       All routines are real*8.
c
c
c       Created by Tom Loredo, 2006-03-07 (from a late '80s codebase)
c
c       Copyright (c) 2006, Tom Loredo
c
c-----------------------------------------------------------------------
        subroutine mvnini(n, S, L, norm, inv, h, err)

c---    This calculates the matrix L needed to sample a multivariate 
c       normal of n parameters with covariance matrix S (or its
c       inverse).   It just performs a Cholesky decomposition of C
c       to find the matrix L such that
c                   S = L L^
c       where L^ is the transpose of L. 
c
c       (Multivariate normal deviates can also be calculated using
c       eigenvectors of S, but the Cholesky decomposition method is
c       somewhat faster, since L is lower-triangular.)
c
c       Input:
c               n       Length of vector of deviates needed
c               S       The n by n covariance matrix of the MVN,
c                       or its inverse (**altered if inverse!**)
c               inv     A switch; =0 if S is the covariance matrix,
c                       =1 if it is the inverse
c               h       Vector workspace: real*8 h(n)
c
c       Output:
c               L       The n by n lower triangular Cholesky decomp.
c                       of S.
c               norm    The normalization constant for the mvg,
c                       norm = (2pi)^(-n/2) / sqrt(det S)
c               S       If the input S was the inverse of the covar.
c                       matrix, it is replace with the covar. matrix
c               err     Error flag:  0 = OK, 1 = Cholesky failed,
c                       2 = S inversion failed.
c
c       Algorithm:
c       The algorithm is from LINPACK, as published (in Algol) in
c       Num. Mat. 7, p. 362, 1965 (choldet1).  It is also described
c       (more clearly?) in R. Y. Rubinstein's book, "Simulation and
c       the Monte Carlo Method," p. 68.  The method of storage for 
c       L has been altered from choldet1.
c
c       The commented code allows calculation of larger determinants
c       by writing  det = d1 * 2**d2, with d1 in the range 1/16 to 1.
c
c       The inversion of S (if inv=1) is accomplished with the
c       Gauss-Jordan method, taking advantage of the PDS nature of S.
c       See Wilkinson and Reinsch, "Linear Algebra" (LINPACK gjdef1).
c       Note that inverse covar = - (Hessian at mode)/(func at mode).

        integer n
        real*8 S(n,n), L(n,n), h(n), norm
        integer inv, err, i, j, k
        real*8 d1, x, p, q, det, pi
        parameter (pi = 3.141592653589793d0)
c        integer d2

cf2py  intent(in) n, inv
cf2py  intent(in,out,copy) S
cf2py  intent(out) L, norm, err
cf2py  intent(cache,hide) h

        err = 0

c---  If S is the INVERSE covariance matrix (i.e. the matrix of
c---  second derivatives), invert it with GJ.
        if (inv .eq. 1) then
            do 20 k=n, 1, -1
                p = S(1,1)
                if (p .le. 0.) then
                    err = 2
                    return
                endif
                do 10 i=2, n
                    q = S(i,1)
                    h(i) = q / p
                    if (i .le. k) h(i) = - h(i)
                    do 5 j=2, i
                        S(i-1, j-1) = S(i,j) + q*h(j)
5                   continue
10              continue
                S(n,n) = 1. / p
                do 15 i=2, n
                    S(n,i-1) = h(i)
15              continue
20          continue
            do 30 i=1, n
                do 25 j=i+1,n
                    S(i,j) = S(j,i)
25              continue
30          continue
        endif

c---  Do the Cholesky decomposition recursively.
        d1 = 1.
c        d2 = 0
        do 100 i=1, n
            do 80 j=i, n
                x = S(i,j)
                if (i .gt. 1) then
                    do 35 k=i-1, 1, -1
                        x = x - L(j,k) * L(i,k)
35                  continue
                endif
                if (i .eq. j) then
                    d1 = d1 * x
                    if (x .eq. 0) then
                        err = 1
                        return
                    endif
c40                  if (abs(d1) .ge. 1.) then
c                        d1 = d1 * 0.0625
c                        d2 = d2 + 4
c                        go to 40
c                    endif
c60                  if (abs(d1) .lt. 0.0625) then
c                        d1 = d1 * 16.
c                        d2 = d2 - 4
c                        go to 60
c                    endif
                    if (x .lt. 0.) then
                        err = 1
                        return
                    endif
                    L(i,i) = sqrt(x)
                else
                    L(j,i) = x / L(i,i)
                    L(i,j) = 0.
                endif
80          continue
100     continue

c---  For large matrices, use the commented code above and return
c---  the log of the determinant.
c        det = log(d1) + d2*log(2.)
        det = d1

        norm = n / 2.
        norm = 1. / ((2.*pi)**norm * sqrt(det))

        return
        end

c-----------------------------------------------------------------------
        subroutine mvnsamp(n, mu, L, snsamps, samp)

c---    Return a sample from a multivariate normal distribution.
c
c       Input:
c           n           dimension of MVN dist'n
c           mu          the mean vector (length n)
c           L           n x n lower triangular matrix from mvgini
c           snsamps     vector of n standard normal samples
c
c       Output:
c           samp        an n-vector with the sample
c
c       Algorithm:
c           See Rubinstein's "Simulation and the Monte Carlo Method,"
c           p. 65.

        integer n, i, j
        real*8 mu(n), L(n,n), samp(n), snsamps(n)

cf2py  intent(in) n, snsamps, mu, L
cf2py  intent(out) samp

c---  Initialize to the independent normal samples.
        do 20 i=1, n
            samp(i) = snsamps(i)
20      continue

c---  Build the shifted, correlated deviates using L and mu.
c---  Work down from sample n to avoid needing working space.
        do 60 i=n, 1, -1
            samp(i) = mu(i) + samp(i) * L(i,i)
            do 40 j=1, i-1
                samp(i) = samp(i) + L(i,j) * samp(j)
40          continue
60      continue

        return
        end

c-----------------------------------------------------------------------
        subroutine mvnsampq(n, mu, L, snsamps, samp, Q)

c---    Return a sample from a multivariate normal distribution,
c       and the value of its associated quadratic form.
c
c       Input:
c           n           dimension of MVN dist'n
c           mu          the mean vector (length n)
c           L           n x n lower triangular matrix from mvgini
c           snsamps     vector of n standard normal samples
c
c       Output:
c           samp        an n-vector with the sample
c           Q           the quadratic form corresponding to the sample,
c                       Q = (samp-mu).S^-1.(samp-mu)
c
c       Algorithm:
c           See Rubinstein's "Simulation and the Monte Carlo Method,"
c           p. 65.

        integer n, i, j
        real*8 mu(n), L(n,n), samp(n), Q, snsamps(n)

cf2py  intent(in) n, snsamps, mu, L
cf2py  intent(out) samp, Q

c---  Initialize to the independent normal samples.
        do 20 i=1, n
            samp(i) = snsamps(i)
20      continue

c---  Now build the correlated deviates out of them using L and mu.
c---  Do it starting with n to avoid needing working space.
c---  Calculate Q along the way.
        Q = 0.
        do 60 i=n, 1, -1
            Q = Q + samp(i)**2
            samp(i) = mu(i) + samp(i) * L(i,i)
            do 40 j=1, i-1
                samp(i) = samp(i) + L(i,j) * samp(j)
40          continue
60      continue

        return
        end

c-----------------------------------------------------------------------
        function mvnden(n, samp, mu, norm, L, work)

c---    This returns the density associated with a sample from a
c       multivariate normal.  If the scalar value, Q, of the quadratic
c       form is already known (e.g. from mvnsamplq), it is much quicker
c       to just use norm*exp(-Q/2).
c
c       Arguments:
c           n       dimension of multivariate Gaussian
c           samp    an n-vector sample from the mvg
c           mu      the mean vector for the mvg
c           norm    the normalization constant; see mvgini
c           L       the Cholesky decomp. of the cov. matrix; see mvgini
c           work    a vector of working space
c
c       Algorithm:
c           The argument of the exponential is (samp-mu)S^-1(samp-mu).
c           Here S^-1(samp-mu) is calculated from L using the LINPAK 
c           cholsol algrorithm in Num. Mat. 7, p. 362, 1965, here
c           translated to FORTRAN and simplified for use with a single
c           "right hand side".

        real*8 mvnden
        integer n
        real*8 samp(n), mu(n), L(n,n), work(n), pi, z
        real*8 arg, norm
        integer i, k
        parameter (pi = 3.141592653589793d0)

c---  Don't declare work as cache; for repeated calls, the calling
c---  Python function should just create the space once.
cf2py  intent(in)  n, samp, mu, norm, L, work
cf2py  intent(out) mvnden

c---  Load the work array with (samp-mu).
        do 20 i=1, n
            work(i) = samp(i) - mu(i)
20      continue

c---  Calculate E^-1 (samp-mu) using the cholsol algorithm.
        do 40 i=1, n
            z = work(i)
            do 30 k=i-1, 1, -1
                z = z - L(i,k) * work(k)
30          continue
            work(i) = z / L(i,i)
40      continue
        do 80 i=n, 1, -1
            z = work(i)
            do 60 k=i+1, n
                z = z - L(k,i) * work(k)
60          continue
            work(i) = z / L(i,i)
80      continue

c---  Now dot it into (samp-mu).
        arg = 0.
        do 100 i=1, n
            arg = arg + (samp(i) - mu(i)) * work(i)
100     continue

c---  Now we've got it...
        mvnden = norm * exp(-arg/2.)

        return
        end

c-----------------------------------------------------------------------
        function mvnqform(n, samp, mu, L, work)

c---    "MultiVariate Normal Quadratic FORM"
c
c       Arguments:
c           n       dimension of multivariate Gaussian
c           samp    an n-vector sample from the mvg
c           mu      the mean vector for the mvg
c           L       the Cholesky decomp. of the cov. matrix; see mvgini
c           work    a vector of working space
c
c       Algorithm:
c           The quadratic form is is (samp-mu)S^-1(samp-mu).
c           Here S^-1(samp-mu) is calculated from L using the LINPAK 
c           cholsol algrorithm in Num. Mat. 7, p. 362, 1965, here
c           translated to FORTRAN and simplified for use with a single
c           "right hand side".

        real*8 mvnqform
        integer n
        real*8 samp(n), mu(n), L(n,n), work(n), pi, z
        real*8 arg
        integer i, k
        parameter (pi = 3.141592653589793d0)

c---  Don't declare work as cache; for repeated calls, the calling
c---  Python function should just create the space once.
cf2py  intent(in)  n, samp, mu, L, work
cf2py  intent(out) mvnqform

c---  Load the work array with (samp-mu).
        do 20 i=1, n
            work(i) = samp(i) - mu(i)
20      continue

c---  Calculate E^-1 (samp-mu) using the cholsol algorithm.
        do 40 i=1, n
            z = work(i)
            do 30 k=i-1, 1, -1
                z = z - L(i,k) * work(k)
30          continue
            work(i) = z / L(i,i)
40      continue
        do 80 i=n, 1, -1
            z = work(i)
            do 60 k=i+1, n
                z = z - L(k,i) * work(k)
60          continue
            work(i) = z / L(i,i)
80      continue

c---  Now dot it into (samp-mu).
        arg = 0.
        do 100 i=1, n
            arg = arg + (samp(i) - mu(i)) * work(i)
100     continue

c---  Now we've got it...
        mvnqform = arg

        return
        end

c-----------------------------------------------------------------------
        function mvnkde(ndim, ncomp, pt, nodes, norm, L, work, scale)

c---    "MultiVariate Normal Kernel Density Estimate"

        real*8 mvnkde
        integer ndim, ncomp
        real*8 pt(ndim), nodes(ndim, ncomp), L(ndim,ndim), work(ndim)
        real*8 norm, scale
        integer i
        real*8 Q, scale2
        real*8 mvnqform

c---  Don't declare work as cache; for repeated calls, the calling
c---  Python function should just create the space once.
cf2py  intent(in)  ndim, ncomp, pt, nodes, norm, L, work, scale
cf2py  intent(out) mvnkde

        mvnkde = 0.
        scale2 = scale**2
        do 100 i=1, ncomp
c            print *, i, ndim, nodes(1,i), nodes(2,i)
            Q = mvnqform(ndim, pt, nodes(1,i), L, work) / scale2
            mvnkde = mvnkde + exp(-Q/2.)
100     continue
        mvnkde = mvnkde / (norm*ncomp*scale)

        return
        end

c-----------------------------------------------------------------------
        subroutine mvtini(n, nu, H, L, norm, work, err)

c---    This calculates the matrix L needed to sample a multivariate 
c       Student dist. of n parameters with negative normalized Hessian 
c       matrix H.  It just performs a Cholesky decomposition of
c       H to find the matrix L such that
c                   inv H = L L^
c       where L^ is the transpose of L.
c
c       Note that the MVT may not have finite 2nd moments; this is
c       why the correlation structure is best specified by the
c       normalized Hessian, rather than a covariance matrix.
c
c       Input:
c               n       Length of vector of deviates needed.
c               nu      # of dof for mvt.
c               H       The n by n NORMALIZED -Hessian matrix at mode.
c               ndim    H and L are assumed dimensioned (ndim,ndim).
c               work    Vector workspace: real*8 work(n)
c
c       Output:
c               L       The n by n lower triangular Cholesky decomp.
c                       of H.
c               norm    The normalization constant for the mvt.
c               H       Replaced with its inverse.
c               err     Error flag:  0 = OK, 1 = Cholesky failed,
c                       2 = S inversion failed.
c
c       For very large nu, it may be necessary to compute norm using
c       logarithms.

        integer n
        real*8 H(n,n), L(n,n), work(n), norm
        integer nu, err, i, j
        real*8 fac, hgamma, rtpi
        parameter (rtpi = 1.7724538509d0)

cf2py  intent(in) n, nu
cf2py  intent(in,out) H
cf2py  intent(out) L, norm, err
cf2py  intent(cache,hide) work

c---  Do the decomposition of H.
        call mvnini(n, H, L, norm, 1, work, err)
        if (err .ne. 0) return

c---  Rescale L by sqrt((n+nu)/nu), since the quadratic form for mvt
c---  is -nu/(n+nu) times its normalized Hessian, and L is the
c---  "square root" of its inverse.
        fac = n + nu
        fac = sqrt(fac/nu)
        do 20 i=1, n
            do 10 j=1, i
                L(i,j) = fac * L(i,j)
10          continue
20      continue

c---  Calculate the normalization.  It is 
c---     nu**(nu/2) gam((nu+n)/2) rt(det(Q)) / ( pi**(n/2) gam(nu/2) )
c---  Note the pi factor and the det factor are returned by mvgini,
c---  with extra powers of 2, and a missing "fac" factor, since
c---  det(Q) = (nu/(n+nu))**n det(H).
        norm = norm * 2.**(n/2.) / fac**n
        norm = norm * nu**(nu/2.) * hgamma(nu+n) / hgamma(nu)

        return
        end

c-----------------------------------------------------------------------
        subroutine mvtsamp(n, mu, L, nu, snsamps, gamsamp, samp)

c---    Return a sample from a multivariate Student t distribution.
c
c       Input:
c           n           length of vector of samples
c           mu          the mean vector (length n)
c           L           n x n lower triangular matrix from mvtini
c           nu          # of dof for t distribution
c           snsamps     vector of n standard normal samples
c           gamsamp     sample from a gamma dist'n with unit scale
c                       and shape parameter alpha = nu/2
c
c       Output:
c           samp        an n-vector of the MVT samples
c
c       Algorithm:
c           This is based on a cute, well-known relationship between
c           the multivariate t distribution and the multivariate 
c           Gaussian and univariate inverse gamma distributions.
c           Crudely, mvt is mvg with its covariance matrix scaled
c           by a scalar, and marginalized with respect to this scalar,
c           using the inverse gamma distribution.  See Zellner's
c           book, "An Intro. to Bayesian Inference in Econometrics,"
c           for a derivation and analysis; the algorithm is outlined
c           in L. Bauwens, "Bayesian Full Information Analysis of 
c           Simultaneous Equation Models Using Integration by
c           Monte Carlo."
c
c           An alternative, less efficient, method is described in
c           H. K. van Dijk and T. Kloek, J. Econometrics, 14, p. 307,
c           1980.  I think it does something similar, using a 
c           relationship between the normal and half-integral gamma 
c           distributions to scale by inverse gamma.

        integer n
        real*8 mu(n), L(n,n), samp(n), sig, snsamps(n), gamsamp
        integer nu, i

cf2py  intent(in) n, mu, L, nu, snsamps, gamsamp
cf2py  intent(out) samp

c---  First sample sig from an inverse gamma distribution with unit 
c---  scale.  If x is distributed as gamma with scale = 2/(nu*s*s) and
c---  alpha = nu/2, then sig = 1/rt(x) is distributed as IG
c---  with nu dof and scale s.
        sig = 2. * gamsamp / nu
        sig = 1. / sqrt(sig)

c---  Now sample mvg, and scale the results by sig.
        call mvnsamp(n, mu, L, snsamps, samp)
        do 20 i=1, n
            samp(i) = mu(i) + sig * (samp(i)-mu(i))
20      continue

        return
        end

c-----------------------------------------------------------------------
        subroutine mvtsampq(n, mu, L, nu, snsamps, gamsamp, samp, Q)

c---    Return a sample from a multivariate Student t distribution,
c       and the value of its associated quadratic form.
c
c       Input:
c           n           length of vector of samples
c           mu          the mean vector (length n)
c           L           n x n lower triangular matrix from mvtini
c           nu          # of dof for t distribution
c           snsamps     vector of n standard normal samples
c           gamsamp     sample from a gamma dist'n with unit scale
c                       and shape parameter alpha = nu/2
c
c       Output:
c           samp        an n-vector of the MVT samples
c           Q           the quadratic form corresponding to the sample,
c                       Q = (samp-mu).V.(samp-mu).  The matrix V is
c                       related to the normed Hessian, H, by
c                           V = -2 nu /(n+nu) H 
c
c       Algorithm:
c           This is based on a cute, well-known relationship between
c           the multivariate t distribution and the multivariate 
c           Gaussian and univariate inverse gamma distributions.
c           Crudely, mvt is mvg with its covariance matrix scaled
c           by a scalar, and marginalized with respect to this scalar,
c           using the inverse gamma distribution.  See Zellner's
c           book, "An Intro. to Bayesian Inference in Econometrics,"
c           for a derivation and analysis; the algorithm is outlined
c           in L. Bauwens, "Bayesian Full Information Analysis of 
c           Simultaneous Equation Models Using Integration by
c           Monte Carlo."
c
c           An alternative, less efficient, method is described in
c           H. K. van Dijk and T. Kloek, J. Econometrics, 14, p. 307,
c           1980.  I think it does something similar, using a 
c           relationship between the normal and half-integral gamma 
c           distributions to scale by inverse gamma.

        integer n
        real*8 mu(n), L(n,n), samp(n), snsamps(n), gamsamp, Q, sig
        integer nu, i

cf2py  intent(in) n, mu, L, nu, snsamps, gamsamp
cf2py  intent(out) samp, Q

c---  First sample sig from an inverse gamma distribution with unit 
c---  scale.  If x is distributed as gamma with scale = 2/(nu*s*s) and
c---  alpha = nu/2, then sig = 1/rt(x) is distributed as IG
c---  with nu dof and scale s.
        sig = 2. * gamsamp / nu
        sig = 1. / sqrt(sig)

c---  Now sample mvg, and scale the results by sig.
        call mvnsamp(n, mu, L, snsamps, samp)
        do 20 i=1, n
            samp(i) = mu(i) + sig * (samp(i)-mu(i))
20      continue
        Q = sig**2 * Q

        return
        end

c-----------------------------------------------------------------------
        function mvtden (n, samp, mu, nu, norm, L, work)

c---    This returns the density of the mvt distribution with n 
c       parameters and nu degrees of freedom, given a sample and
c       the parameters.  It is patterned after mvnden, above.
c
c       Arguments:
c           n       # of params.
c           samp    The n sampled parameters.
c           mu      The means.
c           nu      # of dof.
c           norm    normalization constant, from mvtini.
c           L       From mvtini, above.
c           work    A work array of length n.
c
c       The formula is:
c
c               p(x) = norm * (nu + Q)**((n+nu)/2)
c
c       where  
c           norm =  nu**(nu/2) gam((nu+n)/2) rt(det(V)) / 
c                    ( pi**(n/2) gam(nu/2) )
c           Q = (x-mode).V.(x-mode)
c           x is a vector of parameter values
c           V is a PDS matrix related to the normalized Hessian by
c               V = -2 nu/(n+nu) H   
c                 = -2 nu/(n+nu) [(d/dxi)(d/dxj) p(x)]/p(x) at mode
c
c       For large nu, it may be necessary to compute logarithmically;
c       there are ratios of large numbers involved.


        real*8 mvtden
        integer n, nu
        real*8 samp(n), mu(n), L(n,n), work(n), z
        real*8 arg, norm
        integer i, k

cf2py  intent(in) n, samp, mu, nu, norm, L, work
cf2py  intent(out) mvtden

c---  Load the work array with (samp-mu).
        do 20 i=1, n
            work(i) = samp(i) - mu(i)
20      continue

c---  Calculate V (samp-mu) using the cholsol algorithm.
        do 40 i=1, n
            z = work(i)
            do 30 k=i-1, 1, -1
                z = z - L(i,k) * work(k)
30          continue
            work(i) = z / L(i,i)
40      continue
        do 80 i=n, 1, -1
            z = work(i)
            do 60 k=i+1, n
                z = z - L(k,i) * work(k)
60          continue
            work(i) = z / L(i,i)
80      continue

c---  Now dot it into (samp-mu).
        arg = 0.
        do 100 i=1, n
            arg = arg + (samp(i) - mu(i)) * work(i)
100     continue

c---  Now we've got it...
        mvtden = norm / sqrt( (nu + arg)**(n+nu) )

        return
        end

c-----------------------------------------------------------------------
        function mvtdenq (n, nu, norm, Q)

c---    This returns the density of the mvt distribution with n 
c       parameters and nu degrees of freedom, given the value,
c       Q, of the quadratic form for the parameters.
c
c       Arguments:
c           n       # of params.
c           nu      $ of dof.
c           norm    normalization constant, from mvtini.
c           Q       scalar value of the quadratic form.
c
c       The formula is:
c
c               p(x) = norm * (nu + Q)**((n+nu)/2)
c
c       where  
c           norm =  nu**(nu/2) gam((nu+n)/2) rt(det(V)) / 
c                    ( pi**(n/2) gam(nu/2) )
c           Q = (x-mode).V.(x-mode)
c           x is a vector of parameter values
c           V is a PDS matrix related to the normalized Hessian by
c               V = - nu/(n+nu) H   
c                 = - nu/(n+nu) [(d/dxi)(d/dxj) p(x)]/p(x) at mode


        real*8 mvtdenq
        real*8 norm, Q
        integer n, nu

cf2py  intent(in) n, nu, norm, Q
cf2py  intent(out) mvtdenq

        mvtdenq = norm / sqrt( (nu + Q)**(n+nu) )

        return
        end

c-----------------------------------------------------------------------
        function hgamma(n)

c---    Return the gamma function of n/2, with n an integer >= 1.
c
c       For n even,  gamma(n/2) = (n/2 - 1)!
c       For n odd,   gamma(n/2) = (n/2-1)*(n/2-2)*...*(1/2)*rt(pi)

        real*8 hgamma
        integer n, ntop, j
        real*8 a(65), gammln, arg
        data ntop,a(1),a(2)/2, 1.7724538509, 1./
        save a, ntop

        if (n .le. ntop) then
            hgamma = a(n)
        else if (n .le. 64) then
            do 20 j=ntop+1, n
                a(j) = (j/2. - 1.) * a(j-2)
20          continue
            ntop = n
            hgamma = a(n)
        else if (n .gt. 64) then
            arg = n / 2.
            hgamma = exp(gammln(arg))
        endif

        return
        end

c-----------------------------------------------------------------------
        function gammln(xx)

c---    Taken from 'Numerical Recipes,' this function calculates
c       the log gamma function.

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
