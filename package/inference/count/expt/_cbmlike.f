c-----------------------------------------------------------------------
c
c       CBMLike:  Constant Background Marginal Likelihood
c
c       Functions for calculating the Poisson likelihood for
c       a signal rate when data are collected with a constant
c       unknown background rate, estimated by observations
c       in "off-source" intervals.
c
c       NOTE:  These docs (at the start of this file) are not
c       completely updated wrt the following code; trust each
c       subprogram's docs over this header!
c
c       The CBM likelihood is proportional to equation (12.18) of
c       "The Promise of Bayesian Inference for Astrophysics"
c       (in *Statistical Challenges in Modern Astronomy*, ed.
c       E. Feigelson and G. Babu, Springer-Verlag, 1992,
c       pages 275-305).  That equation gives the normalized marginal
c       posterior density for the signal.  These subroutines calculate 
c       the marginal likelihood for the signal.  This is not
c       normalized over signal values.  It is also not
c       normalized over data (on/off count) values; it is
c       missing data-dependent factors that would make it a
c       rigorously normalized sampling distribution.
c
c       The data that the user must provide are:
c           n_on        Integer # of counts seen on-source
c           t_on        On-source interval size (eg time for
c                       rate measurements, time*area for flux
c                       measurements)
c           n_off        Integer # of counts seen off-source
c           t_off        Off-source interval size
c
c       The user must also set two parameters:
c           s                The signal strength for which the likelihood
c                       is to be calculated
c           cut                A parameter specifying the accuracy of
c                       the calculation, in [0., 1.). 
c
c       The calculation is done by summing terms in the sum
c       given in equation (12.18) mentioned above.  If CUT = 0.,
c       the sum is done with all (n_on + 1) terms, giving an
c       exact result (to machine precision).  Otherwise, the
c       sum starts with the largest term (call it LTERM), and
c       adds decreasing terms until the new terms fall below
c       CUT*LTERM.  Terms in the sum are calculated recursively.
c
c       Note that the C_i factor in the sum is the probability
c       that i of the n_on counts are actually from the signal,
c       not the background.  For large values of n_on or n_off,
c       many C_i coefficients will be negligibly small, possibly
c       causing an underflow.  There is no use keeping such terms,
c       so I strongly recommend using CUT > 0.  This obviously
c       also saves computational time, though in none of my own
c       applications has this likelihood calculation demanded
c       excessive time (with, say, CUT = 1.e-5).
c
c       When the number of counts on- or off-source is large,
c       the value of the likelihood can be very small, possibly
c       smaller than allowed by Fortran77 double precision.
c       Thus these subroutines calculate a LOGARITHM of the 
c       SCALED likelihood.  The user must first calculate a scale
c       factor using SL_SCL, and then calculate the scaled log
c       likelihood using SLLIKE.  The scale factor depends only
c       on the data, not on the values of S or CUT, hence the
c       separate functions.  
c
c       Thus, once the data are known, the user should execute
c       the following:
c
c           double precision sl_scl
c           double precision lscl, t_on, t_off
c           integer n_on, n_off
c
c           lscl = sl_scl(n_on, t_on, n_off, t_off)
c
c       Then for each value of S for which the likelihood is needed,
c       execute (for example):
c
c           double precision sllike
c           double precision s, cut, sll
c           integer nt
c
c           cut = 1.e-5
c           sll = sllike(s, n_on, ta_on, n_off, ta_off, lscl, cut, nt)
c
c       SLL will be the LOG of the SCALED likelihood.  If you need
c       the actual value of the log likelihood for some reason,
c       just add LSCL to SLL.  But since the value of LSCL depends
c       only on the data and not on the value of S, no common
c       calculations (parameter estimates, credible regions, Bayes
c       factors, odds ratios) are afffected by using only the
c       scaled value.
c
c       NT is set by the function and indicates the number of terms
c       used in the calculation.  NT <= N_ON + 1.  It may be of
c       interest for monitoring how many terms are necessary at
c       a particular level of accuracy (value of CUT).
c
c       If the user knows beforehand that scaling is not necessary
c       (e.g., the n_on and n_off values are relatively small),
c       the small amount of overhead associated with sl_scl can be
c       eliminated by eliminating the call to sl_scl and simply
c       setting LSCL=0.
c
c       If the derivative of the scaled log likelihood wrt S is also
c       needed, use the following in place of the last set of
c       Fortran77 lines:
c
c           double precision s, cut, sll, dsll
c           integer nt
c
c           cut = 1.e-5
c           call cbmlike (s, n_on, ta_on, n_off, ta_off, lscl, cut, 
c                         sll, dsll, nt)
c
c       SLL is as above; DSLL is its derivative wrt S.
c
c       Note that these functions are DOUBLE PRECISION.
c
c       Also, note that these functions call a function, GAMMLN,
c       that returns the log of the gamma function.  This must be
c       a real*8 function of a real*8 argument.  An appropriate
c       implementation is the GAMMLN function in the "Numerical
c       Recipes" books, changed to be real*8.  It is included
c       here.
c
c       Algorithm and code by Tom Loredo, 1992.
c       Modifications:
c           09 Apr 97:  Added documentation section
c           21 Apr 97:  Added derivative calculation
c           14 May 97:  Added remark about requiring GAMMLN to docs
c           13 Dec 04:  Modified for Python
c           03 Jan 06:  Further changes for Python; n_off=0 bugfix
c                       in sigctp
c
c-----------------------------------------------------------------------
        subroutine sigctp(n_on, t_on, n_off, t_off, c)

c       SIGnal CounT Probabilities
c
c       Calculate the c coefficients in the series expression for the
c       for the marginal likelihood of the signal.  The dimension of c
c       is n_on + 1, and the first element is c(0).
c
c       c(i) is the probability that i of the n_on counts are from
c       the signal.
c
c       The coefficients are calculated recursively, going up and
c       down from the largest one so roundoff error accumulates only
c       on the negligible terms.
c
c       This assumes a flat prior for the background rate.

cf2py   intent(in) n_on, t_on, n_off, t_off
cf2py   intent(out) c

        integer n_on, n_off, i, mode
        real*8 t_on, t_off, c(0:n_on), tfac, sum, n

        tfac = 1. + t_off/t_on
        n = n_on + n_off
        mode = n_on - t_on*n_off/t_off + 1.
        mode = max(mode, 0)
        mode = min(mode, n_on)
        sum = 1.
        c(mode) = 1.
        do 20 i=mode-1, 0, -1
            c(i) = c(i+1) * ( ((n - i)/(n_on - i)) / tfac )
            sum = sum + c(i)
20      continue
        do 30 i=mode+1, n_on
            c(i) = c(i-1) * (tfac * ((n_on - i + 1.)/(n - i + 1.)) )
            sum = sum + c(i)
30      continue
        do 40 i=0, n_on
            c(i) = c(i) / sum
40      continue

        return
        end

c-----------------------------------------------------------------------
        function psigc(s, n_on, t_on, c)

c---        Probability for SIGnal from Coefficients
c
c       Calculate the probability density for the signal, s, using the
c       c coefficients from coeff, above.

cf2py   intent(in) s, n_on, t_on, c

        integer n_on, i
        real*8 psigc, s, t_on, c(0:n_on), ppois, st

        psigc = 0.
        st = s*t_on
        do 20 i=0, n_on
            psigc = psigc + c(i) * t_on * ppois(i, st)
20      continue

        return
        end

c-----------------------------------------------------------------------
        function slml_offset (n_on, t_on, n_off, t_off)

c---        "Constant Background Marginal Log Likelihood OFFset"
c
c       Return the log of a scale factor for use in SLIKE, to
c       prevent overflow/underflow.

cf2py   intent(in) n_on, t_on, n_off, t_off

c+++  Arguments:
        real*8 slml_offset, t_on, t_off
        integer n_on, n_off

c+++  Locals:
        real*8 stt, nn, nns, disc, arg1, arg2, arg3
        integer km
        real*8 gammln

c---  The factor is just the largest term in the SLIKE series 
c---  for a value of s near the peak.  The logarithm is returned.

c---  stt is the value of the signal, s, at the estimated peak,
c---  times the sum of on and off times.
        stt = (n_on/t_on - n_off/t_off)*(t_on + t_off)

c---  Identify the index, km, for the largest term in the likelihood
c---  by solving a quadratic.
        nn = n_on + n_off
        nns = nn + stt - 1.
        disc = nns*nns + 4.*(nn - n_on*stt)
        if (disc .ge. 0.) then
            km = 0.5* (nns - sqrt(disc)) + 1.
            km = min(km, n_on)
        else
            km = n_on
        endif

c---  If km < 0, use the 0 term
        if (km .le. 0) then
            arg2 = n_on + 1.
            arg1 = arg2 + n_off
            slml_offset = gammln(arg1) - gammln(arg2)
        else
            arg2 = n_on - km + 1.
            arg1 = arg2 + n_off
            arg3 = km + 1.
            slml_offset = km*log(stt) + gammln(arg1) - gammln(arg2) - 
     +               gammln(arg3) - (n_on - n_off*t_on/t_off)
        endif

        return
        end

c-----------------------------------------------------------------------
        subroutine vslml_offset (nc, n_on, t_on, n_off, t_off, offset)

c---        "Vector Constant Background Marginal Log Likelihood OFFSET"
c
c       Return the log of a scale factor for use in SLIKE, to
c       prevent overflow/underflow, for each channel in a
c       vector of on/off measurements

c+++  Arguments:
        real*8 t_on, t_off
        integer nc, n_on(nc), n_off(nc)
        real*8 offset(nc)

c+++  Locals:
        integer i
        real*8 slml_offset

cf2py   intent(in) nc, n_on, t_on, n_off, t_off
cf2py   intent(out) offset

        do 20 i=1, nc
            offset(i) = slml_offset(n_on(i), t_on, n_off(i), t_off)
20        continue

        return
        end

c-----------------------------------------------------------------------
        function slmlike(s, n_on, t_on, n_off, t_off, offset, cut, nt,
     *                         ierr)

c---        "Signal Marginal Log LIKElihood"
c
c       Return the log likelihood for the signal, calculating with all
c       terms up to a factor of CUT from the leading term.
c
c       offset is the log of a scale factor (see slml_offset, above),
c       used to prevent underflow/overflow.
c
c       On return, nt is the number of terms used; nt <= n_on+1.
c
c       exp(slmlike) is NOT normalized w.r.t. s!  It's a log LIKELIHOOD, 
c       not a posterior density.
c
c       It is also NOT normalized w.r.t. n_on, n_off.  It is a
c       likelihood function, missing several data-dependent
c       factors that would make it a normalized sampling dist'n.

cf2py   intent(in) s, n_on, t_on, n_off, t_off, offset, cut
cf2py   intent(out) nt, ierr

c+++  Arguments:
        real*8 slmlike, s, t_on, t_off, offset, cut
        integer n_on, n_off, nt, ierr

c+++  Locals:
        real*8 slike, stt, nn, nns, disc
        real*8 term, lterm, low, arg1, arg2, arg3
        integer km, k
        real*8 gammln

c---  Error flag:
        ierr = 0

c---  Treat s=0 specially, to avoid log(0) calls.
        nt = 1
        if (s .eq. 0.) then
            arg1 = n_on + n_off + 1.
            arg2 = n_on + 1.
            slmlike = gammln(arg1) - gammln(arg2) - offset
            return
        endif

c---  Begin by finding the largest term.
        stt = s * (t_on + t_off)
        nn = n_on + n_off
        nns = nn + stt - 1.
        disc = nns*nns + 4.*(nn - n_on*stt)
        if (disc .ge. 0.) then
            km = 0.5* (nns - sqrt(disc)) + 1.
            km = min(km, n_on)
        else
            km = n_on
        endif
        arg2 = n_on - km + 1.
        arg1 = arg2 + n_off
        arg3 = km + 1.
        lterm = km*log(stt) + gammln(arg1) - gammln(arg2) - 
     +          gammln(arg3) - s*t_on - offset

c---  It is conceivable it could underflow, if we are evaluating
c---  the likelihood extremely far from the peak.  It should
c---  never overflow, but we check anyway.
        if (lterm .lt. -700. .or. lterm .gt. 700.) then
            slike = 0.
            print *, 'Under/overflow lterm in slike:'
            write(*,'(g12.4,2i6,2g12.4)') s, n_on, n_off, t_on, t_off
            ierr = 1
            return
c           pause 'Under/overflow in slmlike!'
        else
            lterm = exp(lterm)
        endif
        low = cut * lterm
        if (cut .gt. 0. .and. low .eq. 0.) then
            print *, 'Vanishing cutoff in slmlike: ',
     *         cut, lterm, cut*lterm
        endif

c----  First work down from the largest term.
        slike = lterm
        term = lterm
        k = km
20          k = k - 1
            if (k .ge. 0.) then
                term = term * (nn - k) * (k + 1.) / (stt*(n_on - k))
                slike = slike + term
                nt = nt + 1
                if (term .gt. low) goto 20
            endif

c---  Now work up from the largest term.
        term = lterm
        k = km
40          k = k + 1
            if (k .le. n_on) then
                term = term * stt * (n_on - k + 1.) / (k*(nn - k + 1.))
                slike = slike + term
                nt = nt + 1
                if (term .gt. low) goto 40
            endif

c---  Take the log, and that's it!
        slmlike = log(slike)

        return
        end

c-----------------------------------------------------------------------
        subroutine vslmlike (nc, s, n_on, t_on, n_off, t_off, offset,  
     *                       cut, slml)

c---        "Vector Signal Marginal Log LIKElihood"
c
c       Return the log likelihood for the signal values in nc channels
c       of on/off data, calculating with all
c       terms up to a factor of CUT from the leading term.
c
c       off is the log of a scale factor (see slml_off, above),
c       used to prevent underflow/overflow.
c
c       exp(slml) is NOT normalized w.r.t. s!  It's a log LIKELIHOOD, 
c       not a posterior density.
c
c       It is also NOT normalized w.r.t. n_on, n_off.  It is a
c       likelihood function, missing several data-dependent
c       factors that would make it a normalized sampling dist'n.

c***  Add ierr return argument.

c+++  Arguments:
        integer nc, n_on(nc), n_off(nc)
        real*8 s(nc), t_on, t_off, offset(nc), cut, slml(nc)

cf2py   intent(in) nc, s, n_on, t_on, n_off, t_off, offset, cut
cf2py   intent(out) slml

c+++  Locals:
        integer i, nt, ierr
        real*8 slmlike

        do 20 i=1, nc
            slml(i) = slmlike(s(i), n_on(i), t_on, n_off(i), t_off, 
     *                        offset(i), cut, nt, ierr)
20        continue

        return
        end

c-----------------------------------------------------------------------
        subroutine slmliked (s, n_on, t_on, n_off, t_off, offset, cut,
     *                      sll, dsll, nt, ierr)

c---        "Signal Marginal Log LIKElihood and Derivative"
c
c       Return the scaled log likelihood for the signal (SLL), 
c       calculating with all terms up to a factor of CUT from the 
c       leading term.  Also calculate the derivative of SLL
c       with respect to s (returned as DSLL).
c
c       offset is the log of a scale factor (see slml_off, above),
c       used to prevent underflow/overflow.
c
c       On return, nt is the number of terms used; nt <= n_on+1.
c
c       SLL is NOT normalized w.r.t. s!  It's a log LIKELIHOOD, not
c       a posterior density.
c
c       noff=0 gives nan

cf2py   intent(in) s, n_on, t_on, n_off, t_off, offset, cut
cf2py   intent(out) sll, dsll, nt, ierr

c+++  Arguments:
        real*8 s, t_on, t_off, offset, cut, sll, dsll
        integer n_on, n_off, nt, ierr

c+++  Locals:
        real*8 slike, dsl, stt, nn, nns, disc
        real*8 term, lterm, low, arg1, arg2, arg3
        integer km, k
        real*8 gammln

c---  Error flag:
        ierr = 0

c---  Treat s=0 specially, to avoid log(0) calls.
        nt = 1
        if (s .eq. 0.) then
            arg1 = n_on + n_off + 1.
            arg2 = n_on + 1.
            sll = gammln(arg1) - gammln(arg2) - offset
            arg1 = arg1 - 1.
            arg2 = arg2 - 1.
            term = log(1. + t_off/t_on) + gammln(arg1) - gammln(arg2)
     *          - offset
            dsll = t_on * (exp(term-sll) - 1.)
            return
        endif

c---  Begin by finding the largest term.
        stt = s * (t_on + t_off)
        nn = n_on + n_off
        nns = nn + stt - 1.
        disc = nns*nns + 4.*(nn - n_on*stt)
        if (disc .ge. 0.) then
            km = 0.5* (nns - sqrt(disc)) + 1.
            km = min(km, n_on)
        else
            km = n_on
        endif
        arg2 = n_on - km + 1.
        arg1 = arg2 + n_off
        arg3 = km + 1.
        lterm = km*log(stt) + gammln(arg1) - gammln(arg2) - 
     +          gammln(arg3) - s*t_on - offset
        if (lterm .lt. -700. .or. lterm .gt. 700.) then
            slike = 0.
            print *, 'Under/overflow lterm in cbmllike:'
            write(*,'(g12.4,2i6,2g12.4)') s, n_on, n_off, t_on, t_off
            ierr = 1
            return
c           pause 'Under/overflow in CBMLLIKE!'
        else
            lterm = exp(lterm)
        endif
        low = cut * lterm
        if (cut .gt. 0. .and. low .eq. 0.) then
            print *, 'Vanishing lterm in CBMLLIKE: ',
     *         cut, lterm, cut*lterm
        endif
        dsl = lterm * (km/s - t_on)

c----  First work down from the largest term.
        slike = lterm
        term = lterm
        k = km
20            k = k - 1
            if (k .ge. 0.) then
                term = term * (nn - k) * (k + 1.) / (stt*(n_on - k))
                slike = slike + term
                dsl = dsl + term*(k/s - t_on)
                nt = nt + 1
                if (term .gt. low) goto 20
            endif

c---  Now work up from the largest term.
        term = lterm
        k = km
40            k = k + 1
            if (k .le. n_on) then
                term = term * stt * (n_on - k + 1.) / (k*(nn - k + 1.))
                slike = slike + term
                dsl = dsl + term*(k/s - t_on)
                nt = nt + 1
                if (term .gt. low) goto 40
            endif

c---  Take the log, and that's it!
        dsll = dsl / slike
        sll = log(slike)

        return
        end

c-----------------------------------------------------------------------
        function ppois(n,dm)

c---    This returns the probability of sampling n events from a
c       Poisson distribution with mean dm.
c
c       P(n,m) = exp(-dm) dm**n / n!
c
c       The factorial part is from Numerical Recipes.

cf2py   intent(in) n, dm

        real*8 ppois, dm, a, dnfac, arg, gammln
        integer n, ntop, j
        dimension a(33)
        data a(1),ntop/1.,0/
        save a, ntop

        if (dm .eq. 0.) then
            if (n .eq. 0.) then
                ppois = 1.
            else
                ppois = 0.
            endif
        else if (n .le. ntop) then
            dnfac = a(n+1)
            ppois = dm**n * exp(-dm) / dnfac
        else if (n .le. 32) then
            do 20 j=ntop+1, n
                a(j+1) = j * a(j)
20          continue
            ntop = n
            dnfac = a(n+1)
            ppois = dm**n * exp(-dm) / dnfac
        else if (n .gt. 32) then
            dnfac = gammln(dble(n)+1.)
            arg = n * log(dm) - dm - dnfac
            ppois = exp(arg)
        endif

        return
        end

c-----------------------------------------------------------------------
        function lppois(n,dm)

c---    This returns the log probability of sampling n events from a
c       Poisson distribution with mean dm.
c
c       P(n,m) = n*log(dm) - dm -log(n!)
c
c       The factorial part is from Numerical Recipes.

cf2py   intent(in) n, dm

        real*8 lppois, dm, a, b, dnfac, gammln
        integer n, ntop, j
        dimension a(33), b(33)
        data a(1),ntop/1.,0/
        save a, ntop

        if (dm .eq. 0.) then
            if (n .eq. 0.) then
                lppois = 0.
            else
                lppois = -1.d100
            endif
        else if (n .le. ntop) then
            dnfac = b(n+1)
            lppois = n*log(dm) - dm - dnfac
        else if (n .le. 32) then
            do 20 j=ntop+1, n
                a(j+1) = j * a(j)
                b(j+1) = log(a(j+1))
20          continue
            ntop = n
            dnfac = b(n+1)
            lppois = n*log(dm) - dm - dnfac
        else if (n .gt. 32) then
            dnfac = gammln(dble(n)+1.)
            lppois = n*log(dm) - dm - dnfac
        endif

        return
        end

c-----------------------------------------------------------------------
        function gammln(xx)

c---    Adapted from 'Numerical Recipes,' this function calculates
c       the log gamma function.

cf2py   intent(in) xx

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
