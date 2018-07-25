	program moments

c---	Calculate moments of 1 or 2 variables from MCMC output.
c	Use only the equilibrated part of the time series, as
c	determined by DIAG.  Write out a set of approximately
c	independent samples based on the correlation length of
c	the variables.
c
c	7 Jul 93  TJL

	integer maxt, maxp, maxo, maxr, maxm, maxs
	parameter (maxt = 16384, maxp = 15, maxr = 50, maxs = 500)
	parameter (maxo = 5, maxm = (maxo*(maxo+3))/2)
	integer npar, nruns, lrun, intvl, nt, d, nmom, m, ntot
	integer time(maxt), ivals(maxr), del(maxr)
	real*8 pars(maxp), vals(maxt), xvals(maxt), yvals(maxt)
	real*8 wrk1(maxt), wrk2(maxt), xsamp(maxs), ysamp(maxs)
	real*8 rcl(maxr,maxm), rmom(maxr,maxm), rsig(maxr,maxm)
	real*8 mom(maxm), msig(maxm)
	real*8 clmax(maxr)
	integer p1, p2
	common /parcb/ p1, p2

	integer r, p, ord, n, t, ivals(maxr), nreq, sper, skip, ns
	real*8 clen, var
	character*60 dfile, dgfile, ofile
	logical eq(maxr)
	real*8 cfuns

c---  Read params.
	p2 = 0
	call opar(1,'moments.pf')
	call rdpw('data.file',dfile)
	call rdpw('diag.file',dgfile)
	call rdpw('output.file',ofile)
	call rdpi('dimension',d)
	if (d .lt. 1 .or. d .gt. 2) pause '1- or 2-d only!'
	call rdpia('params',p1,d)
	call rdpi('order',ord)
	nmom = ord*(ord+3)/2
	if (nmom .gt. maxm) pause 'Too many moments!'
	call cpar

c---  Read the diag header, and determine how many steps to delete.
c---  Delete enough so that *all* parameters are separately in EQ.
	call opar(1,dgfile)
	call rdpi('params',npar)
	if (npar .gt. maxp) pause 'Param # too high!'
	call rdpi('runs',nruns)
	call rdpi('lrun',lrun)
	call rdpi('intvl',intvl)
	ntot = lrun / intvl + 1
	call rdpi('nreq',nreq)
	if (nreq .le. 0) then
	    print *, 'No runs in EQ.'
	    stop
	endif
	call rdplvec('eq',eq,nruns)
	do 40 p=1, npar
	    call rdpivec('r.del', ivals, nruns)
	    do 20 r=1, nruns
	        ivals(r) = ivals(r)/intvl + 1
	        if (p .eq. 1) then
	            del(r) = ivals(r)
	        else
	            del(r) = max(del(r),ivals(r))
	        endif
c	        del(r) = max(del(r),ntot/2)
20	    continue
40	continue
	call cpar
	    
c---  Read the data header.  Check for consistency with diag info.
	call opar(1,dfile)
	call rdpi('npar',n)
	if (n .ne. npar) pause 'Bad npar in data file!'
	call rdpi('nruns',n)
	if (n .ne. nruns) pause 'Bad nruns in data file!'
	call rdpi('lrun',n)
	if (n .ne. lrun) pause 'Bad lrun in data file!'
	call rdpi('intvl',n)
	if (n .ne. intvl) pause 'Bad intvl in data file!'

c===  Go thru the data 1 run at a time, using only EQ runs.
	ns = 0
	sper = maxs / nreq
	do 200 r=1, nruns
	    clmax(r) = 0.
	    if (.not. eq(r)) goto 200
	    print *, 'Run # ',r,' of ',nruns,', deleting ',del(r)

c===  Read in the evenly spaced samples.  Skip the noneq part.
c===  Keep only the 1 or 2 params requested.
	    call rdpi('time',time(1))
	    if (time(1) .ne. 0) pause 'Run alignment problem!'
	    nt = ntot - del(r)
	    do 60 t=1, del(r)
	        call rdpvec('samp',pars,npar)
60	    continue
	    do 80 t=1, nt
	        if (mod(t,1000) .eq. 0) print *, '  Reading sample ',t
	        call rdpi('time',time(t))
	        call rdpvec('samp',pars,npar)
	        xvals(t) = pars(p1)
	        if (d .eq. 2) yvals(t) = pars(p2)
80	    continue

c===  Now work thru the moments one at a time.
	    do 160 m=1, nmom

c>>>  First, make a time series of values of the moment.
	        do 100 t=1, nt
	            vals(t) = cfuns(m,xvals(t),yvals(t))
c	            write(9,'(2i5,3(1pg11.3))')m,t,xvals(t),yvals(t),
c     *  vals(t)
100	        continue

c>>>  Now find the mean, sigma, and corr length with an AR fit.
	        rmom(r,m) = 0.
	        do 120 t=1, nt
	            rmom(r,m) = rmom(r,m) + vals(t)
120	        continue
	        rmom(r,m) = rmom(r,m) / nt
	        do 140 t=1, nt
	            vals(t) = vals(t) - rmom(r,m)
140	        continue
	        call bestar (nt, vals, ord, var, clen, wrk1, wrk2)
	        if (ord .ge. 30) print *, 'Order >= 30!'
	        rsig(r,m) = sqrt(var)
	        rcl(r,m) = clen*intvl
	        clmax(r) = max(clmax(r),rcl(r,m))

c>>>  Go do the next moment.
	        write(*,225) m, rmom(r,m), rsig(r,m)
160	    continue

c===  Collect some nearly independent samples, then go for next run.
	    skip = clmax(r)
	    skip = max(skip, nt/sper)
	    do 180 t=1, nt
	        if (mod(t,skip) .eq. 0) then
	            ns = ns + 1
	            xsamp(ns) = xvals(t)
	            if (d .eq. 2) ysamp(ns) = yvals(t)
	        endif
180	    continue
200	continue
	call cpar
225	format('  Moment ',i3,', = ',1pg12.4,' +- ',1pg12.4)

c---  Find the moments & sigs averaged over all runs.
	do 240 m=1, nmom
	    mom(m) = 0.
	    msig(m) = 0.
	    do 230 r=1, nruns
	        if (eq(r)) then
	            mom(m) = mom(m) + rmom(r,m)
	            msig(m) = msig(m) + rsig(r,m)**2
	        endif
230	    continue
	    mom(m) = mom(m) / nreq
	    msig(m) = msig(m) / nreq
	    msig(m) = sqrt(msig(m))
240	continue

c---  Write out the results.
	call openf(1,ofile)
	write(1,300) npar, nruns, lrun, intvl, nreq, 
     *          d, p1, p2, ord, nmom
300	format('params ',i4,/,
     *         'runs   ',i4,/,
     *         'lrun   ',i6,/,
     *         'intvl  ',i4,/,
     *         'nreq   ',i4,/,/,
     *         'dimen  ',i4,/,
     *         'p1.p2  ',2i4,/,
     *         'order  ',i4,/,
     *         'nmom   ',i4,/)

	call fwrows ('moments ', mom, nmom, 5, 1, '(1pg12.4)')
	write(1,*)
	call fwrows ('sigs    ', msig, nmom, 5, 1, '(1pg12.4)')
	write(1,*)
	call fwrows ('clens   ', clmax, nruns, 5, 1, '(1pg12.4)')
	write(1,*)

c---  Now write the set of almost-independent samples.
	write(1,'(a,i6)') 'nsamp  ',ns
	write(1,*)
	call fwrows ('xsamp   ', xsamp, ns, 5, 1, '(1pg12.4)')
	if (d .eq. 2) then
	    write(1,*)
	    call fwrows ('ysamp   ', ysamp, ns, 5, 1, '(1pg12.4)')
	endif
	close(1)

	end

c-----------------------------------------------------------------------
        function cfuns(i, x, y)

c---    Constraint FUNctionS.
c
c	This returns the value of the function for the i'th moment, 
c	numbered thru various orders consecutively thru the following 
c	table:
c
c	order   terms
c	  0       1
c	  1       x        y
c	  2       x**2     xy      y**2
c         3       x**3     x**2y   xy**2     y**3
c         4       x**4     x**3y   x**2y**2  xy**3  y**4
c	etc.
c
c	The "1" is term 0, "x" is term 1, etc., so, for example,
c	(x**2 * y) is returned for i=7.
c	
c	There are O*(O+3)/2 terms (not counting term 0)
c	thru order O.

c---  Arguments:
        real*8 cfuns, x, y
        integer i

c---  Locals:
        real*8 s
        integer ord, m, n

c---  i=0 is easy!
	if (i .eq. 0) then
	    cfuns = 1.
	    return
	endif

c---  First find the order, given by the integer part of the solution
c---  of a quadratic for the number of terms before ord.
	s = 1. + 8.*i
	ord = 0.5 * (sqrt(s)-1.)

c---  Find the powers for x**m * y**n.
	m = (ord*(ord+3))/2 - i
	n = ord - m

c---  Here's the value.
	if (m .eq. 0) then
	    cfuns = y**n
	else if (n .eq. 0) then
	    cfuns = x**m
	else
	    cfuns = x**m * y**n
	endif
c        write(9,'(4i3,3(1pg12.4))') i,ord,m,n,x,y,cfuns

        return
        end

