	program acorr

c---	Calculate the avg. autocorrelation function from several run 
c	of Metropolis samples.
c
c	04 May 93  TJL

	integer maxt, maxp
	parameter (maxt = 8192, maxp = 15)
	integer np, npar, nruns, lrun, intvl, nt, i, lag, ind
	integer t(maxt), md, run, tmin, nuse
	real*8 pars(maxp), vals(maxt), mean, msqr, c, sig, c0, var0, var
	real*8 dc, rn, ebar, esig, mass(maxp), rnr
	real*8 corr(0:maxt), sc(0:maxt)
	common /parcb/ ebar, esig

	print *, 'Enter param #:'
	read(*,*) np
	print *, 'Enter minimum t:'
	read(*,*) tmin
	call opar(1,'hmc.dat')
	call rdpi('npar',npar)
	if (np .gt. npar .or. npar .gt. maxp) pause 'Param # too high!'
	call rdpi('nruns',nruns)
	if (nruns .eq. 1) then
	    print *, 'Only one run!  Quitting...'
	    stop
	endif
	print *, 'Runs: ',nruns
	call rdpi('lrun',lrun)
	call rdpi('intvl',intvl)
	nt = lrun / intvl + 1
	if (nt .gt. maxt) pause 'Too many times!'
	if ((nt-1)*intvl .ne. lrun) then
	    print *, 'Skipping final sample.'
	endif
	write(*,'(a,i5,a)') 'Reading ',nt,' samples.'
	call rdpi('mdstps',md)
	call rdpda('ebar.esig',ebar,2)
	call rdpvec('masses',mass,np)

c---  Initialize.
	do 5 i=0, nt-1
	    corr(i) = 0.
	    sc(i) = 0.
5	continue

c---  Loop over the runs.
	do 400 run=1, nruns

c---  Read in the evenly spaced samples.
	    nuse = 0
15	    call rdpi('time',t(1))
	    if (t(1) .ne. 0) goto 15
	    call rdpvec('samp',pars,npar)
	    if (t(1) .ge. tmin) then
	        nuse = nuse + 1
	        vals(nuse) = pars(np)
	    endif
	    do 20 i=2, nt
	        call rdpi('time',t(i))
	        call rdpvec('samp',pars,npar)
	        if (t(i) .ge. tmin) then
	            nuse = nuse + 1
	            vals(nuse) = pars(np)
	        endif
20	    continue

c---  Find the mean & mean square for the run.
	    mean = 0.
	    msqr = 0.
	    do 40 i=1, nuse
	        mean = mean + vals(i)
	        msqr = msqr + vals(i)**2
40	    continue
	    mean = mean / nuse
	    msqr = msqr / nuse
	    sig = sqrt(msqr - mean**2)
	    if (mod(run,10) .eq. 0) then
	        write(*,'(a,i4,2(1pg12.4))') 'run, mean, sig: ',
     *            run, mean, sig
	    endif

c---  Do the calculation.  For each lag, find the correlation & 
c---  its uncertainty.
	    do 60 lag=0, nuse-1
	        c = 0.
	        var = 0.
	        do 50 i=1, nuse-lag
	            ind = i + lag
	            dc = (vals(i)-mean) * (vals(ind)-mean)
	            c = c + dc
	            var = var + dc*dc
50	        continue
	        rn = nuse - lag
	        c = c / rn
	        var = var/rn - c**2
	        if (lag .eq. 0) then
	            c0 = c
	            var0 = var
	        endif
	        c = c / c0
	        var = ( var/c0**2 - var0*(c/c0)**2 ) / rn
	        var = max(var,0.)
	        corr(lag) = corr(lag) + c
	        sc(lag) = sc(lag) + c*c
60	    continue
400	continue

c---  Write out the results.
	rnr = nruns
	call openf(1,'acorr.dat')
	write(1,100) np, nruns, md, ebar, esig, mass(np), tmin
100	format('param  ',i4,/,
     *         'runs   ',i4,/,
     *         'mdstps ',i4,/,
     *         'ebar   ',1pg12.4,/,
     *         'esig   ',1pg12.4,/,
     *         'mass   ',1pg12.4,/,
     *         'tmin   ',i5)
	do 500 lag=0, nuse-1
	    corr(lag) = corr(lag) / rnr
	    sig = sqrt( (sc(lag)/rnr - corr(lag)**2)/rnr)
	    write(1,'(3(1pg13.4))') 1.*t(lag+1), corr(lag), sig
500	continue
	close(1)
	end

