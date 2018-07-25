	program diag

c---	Perform diagnostics on Markov chain Monte Carlo output,
c	determining the transient length and the mean of
c	the stationary part of the time series.
c
c	1 Jul 93  TJL

	integer maxt, maxp, maxo, csiz, maxr
	parameter (maxt = 16384, maxp = 15, maxr = 50)
	parameter (maxo = 25, csiz = (maxo*(maxo+1))/2)
	integer np, npar, nruns, lrun, intvl, nt
	integer time(maxt), md
	real*8 pars(maxp), pvals(maxt,maxp), mean, vals(maxt)
	real*8 wrk1(maxt), wrk2(maxt), wrk3(maxt)
	real*8 rmean(maxr,maxp), rar(maxr,maxp), rcl(maxr,maxp)
	real*8 rsp(maxr,maxp)
	real*8 pmean(maxp), par(maxp), psp(maxp)
	integer rdel(maxr,maxp), rord(maxr,maxp), rt0(maxr,maxp)
	real*8 ebar, esig, mass(maxp)
	common /parcb/ ebar, esig

	integer r, p, ord, ndat, t0, t, ivals(maxr), nreq
	real*8 clen, var
	character*60 ifile, ofile
	logical ok, eq(maxr)

c---  Read params.
	call opar(1,'diag.pf')
	call rdpw('input.file',ifile)
	call rdpw('output.file',ofile)
	call cpar

c---  Read the data header.
	call opar(1,ifile)
	call rdpi('npar',npar)
	if (npar .gt. maxp) pause 'Param # too high!'
	call rdpi('nruns',nruns)
	call rdpi('lrun',lrun)
	call rdpi('intvl',intvl)
	nt = lrun / intvl + 1
	if (nt .gt. maxt) pause 'Too many times!'
	if ((nt-1)*intvl .ne. lrun) then
	    print *, 'Skipping final sample.'
	endif
	write(*,'(a,i5,a)') 'Using ',nt,' samples.'
	call rdpi('mdstps',md)
	call rdpda('ebar.esig',ebar,2)
	call rdpvec('masses',mass,np)

c===  Go thru the data 1 run at a time.
	nreq = 0
	do 200 r=1, nruns
	    print *, 'Run # ',r,' of ',nruns,'...'
	    eq(r) = .true.

c===  Read in the evenly spaced samples.
	    call rdpi('time',time(1))
	    if (time(1) .ne. 0) pause 'Run alignment problem!'
	    call rdpvec('samp',pars,npar)
	    do 20 p=1, npar
	        pvals(1,p) = pars(p)
20	    continue
	    do 60 t=2, nt
	        if (mod(t,1000) .eq. 0) print *, '  Reading sample ',t
	        call rdpi('time',time(t))
	        call rdpvec('samp',pars,npar)
	        do 40 p=1, npar
	            pvals(t,p) = pars(p)
40	        continue
60	    continue

c===  Now work thru the parameters one at a time.
	    do 140 p=1, npar

c>>>  First find the correlation length; as a rule, we'll throw
c>>>  out the first 2 correlation lengths.  Estimate it by subtracting
c>>>  the mean and doing an AR fit.
	        mean = 0.
	        do 80 t=1, nt
	            mean = mean + pvals(t,p)
80	        continue
	        mean = mean / nt
	        do 100 t=1, nt
	            vals(t) = pvals(t,p) - mean
100	        continue
	        call bestar (nt, vals, ord, var, clen, wrk1, wrk2)
	        t0 = 2 * clen + 1
	        ndat = nt - t0 + 1
	        if (ndat .lt. 1) then
	            rt0(r,p) = intvl*(t0-1) + time(1)
	            rdel(r,p) = time(nt)
	            rmean(r,p) = 0.
	            rar(r,p) = 0.
	            rsp(r,p) = 0.
	            rord(r,p) = ord
	            rcl(r,p) = clen
	            write(*,210) p, rt0(r,p)
	            eq(r) = .false.
	            goto 140
	        endif
	        rt0(r,p) = time(t0)

c>>>  Now run diagnostics on the rest.
	        call odiag(ndat, pvals(t0,p), rdel(r,p), rmean(r,p), 
     *            rar(r,p), rsp(r,p), rord(r,p), rcl(r,p), 
     *            wrk1, wrk2, wrk3, ok)
	        if (.not. ok) then
	            rdel(r,p) = time(nt)
	            rmean(r,p) = 0.
	            rar(r,p) = 0.
	            rsp(r,p) = 0.
	            rord(r,p) = 0
	            rcl(r,p) = 0.
	            write(*,220) p, rt0(r,p)
	            eq(r) = .false.
	            goto 140
	        endif
	        rdel(r,p) = rdel(r,p) + t0 - 1
	        rdel(r,p) = time(rdel(r,p))

c>>>  Go do the next param or next run, then close the data file.
	        write(*,225) p, rt0(r,p), rdel(r,p)
140	    continue
	    if (eq(r)) nreq = nreq + 1
200	continue
	call cpar
210	format('  Param ',i3,', t0 = ',i6,', t0 too big')
220	format('  Param ',i3,', t0 = ',i6,', no EQ')
225	format('  Param ',i3,', t0 = ',i6,', del = ',i6)
	print *, 'Runs with EQ for all params: ',nreq

c---  Find the mean & sigs averaged over all runs.
	do 240 p=1, npar
	    pmean(p) = 0.
	    par(p) = 0.
	    psp(p) = 0.
	    if (nreq .gt. 0) then
	        do 230 r=1, nruns
	            if (.not. eq(r)) goto 230
	            pmean(p) = pmean(p) + rmean(r,p)
	            par(p) = par(p) + rar(r,p)**2
	            psp(p) = psp(p) + rsp(r,p)**2
230	        continue
	        pmean(p) = pmean(p) / nreq
	        par(p) = par(p) / nreq
	        par(p) = sqrt(par(p))
	        psp(p) = psp(p) / nreq
	        psp(p) = sqrt(psp(p))
	    endif
240	continue

c---  Write out the results.
	call openf(1,ofile)
	write(1,300) npar, nruns, lrun, intvl, md, ebar, esig, nreq
300	format('params ',i4,/,
     *         'runs   ',i4,/,
     *         'lrun   ',i6,/,
     *         'intvl  ',i4,/,
     *         'mdstps ',i4,/,
     *         'ebar   ',1pg12.4,/,
     *         'esig   ',1pg12.4,/,/,
     *         'nreq   ',i4)

	call fwlrows('eq  ', eq, nruns, 20, 1, 'l3')
	write(1,*)

	if (nreq .gt. 0) then
	    call fwrows ('p.mean   ', pmean, npar, 5, 1, '(1pg12.4)')
	    write(1,*)
	    call fwrows ('p.arsig  ', par, npar, 5, 1, '(1pg12.4)')
	    write(1,*)
	    call fwrows ('p.spsig  ', psp, npar, 5, 1, '(1pg12.4)')
	    write(1,*)
	    write(1,*)
	endif

	do 400 p=1, npar
	    write(1,'(a,i5)') 'param  ',p
	    do 310 r=1, nruns
	        ivals(r) = rt0(r,p)
310	    continue
	    call fwirows('r.t0  ', ivals, nruns, 5, 1, '(i11)')
	    do 315 r=1, nruns
	        ivals(r) = rdel(r,p)
315	    continue
	    call fwirows('r.del ', ivals, nruns, 5, 1, '(i11)')
	    do 320 r=1, nruns
	        vals(r) = rmean(r,p)
320	    continue
	    call fwrows('r.mean    ', vals, nruns, 5, 1, '(1pg11.3)')
	    do 330 r=1, nruns
	        vals(r) = rar(r,p)
330	    continue
	    call fwrows('r.arsig   ', vals, nruns, 5, 1, '(1pg11.3)')
	    do 340 r=1, nruns
	        vals(r) = rsp(r,p)
340	    continue
	    call fwrows('r.spsig   ', vals, nruns, 5, 1, '(1pg11.3)')
	    do 350 r=1, nruns
	        ivals(r) = rord(r,p)
350	    continue
	    call fwirows('r.ord ', ivals, nruns, 5, 1, '(i11)')
	    do 360 r=1, nruns
	        vals(r) = rcl(r,p) * intvl
360	    continue
	    call fwrows('r.clen    ', vals, nruns, 5, 1, '(1pg11.3)')
	    write(1,*)
400	continue
	close(1)

	end
