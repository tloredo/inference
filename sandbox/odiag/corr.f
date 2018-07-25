	program corr

c---	Calculate the autocorrelation function of a single run 
c	of Metropolis samples.
c
c	29 Apr 92  TJL

	integer maxt, maxp
	parameter (maxt = 16384, maxp = 15)
	integer np, npar, nruns, lrun, intvl, nt, i, lag, ind
	integer t(maxt), md, run
	real*8 pars(maxp), vals(maxt), mean, msqr, c, sig, c0, var0, var
	real*8 dc, rn, ebar, esig, mass(maxp)
	common /parcb/ ebar, esig

	print *, 'Enter param #:'
	read(*,*) np
	call opar(1,'hmc.dat')
	call rdpi('npar',npar)
	if (np .gt. npar .or. npar .gt. maxp) pause 'Param # too high!'
	call rdpi('nruns',nruns)
	if (nruns .eq. 1) then
	    run = 1
	else
	    print *, 'Enter run # (<= ',nruns,'):'
	    read(*,*) run
	endif
c	if (nruns .gt. 1) print *, 'Using 1 run of',nruns
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

c---  Skip runs, if necessary.
	do 10 i=1, run-1
5	    call rdpi('time',t(1))
	    if (t(1) .ne. 0) goto 5
10	continue

c---  Read in the evenly spaced samples.
15	call rdpi('time',t(1))
	if (t(1) .ne. 0) goto 15
	call rdpvec('samp',pars,npar)
	vals(1) = pars(np)
	do 20 i=2, nt
	    if (mod(i,1000) .eq. 0) print *, 'Reading sample ',i
	    call rdpi('time',t(i))
	    call rdpvec('samp',pars,npar)
	    vals(i) = pars(np)
c	    vals(i) = sin(6.283*t(i)/100.)
20	continue
	call cpar

c---  Find the mean & mean square.
	mean = 0.
	msqr = 0.
	do 40 i=1, nt
	    mean = mean + vals(i)
	    msqr = msqr + vals(i)**2
40	continue
	mean = mean / nt
	msqr = msqr / nt
	sig = sqrt(msqr - mean**2)
	write(*,'(a,2(1pg12.4))') 'Mean, sig: ', mean, sig

c---  Do the calculation.  For each lag, find the correlation & 
c---  its uncertainty.
	call openf(1,'corr.dat')
	write(1,45) np, run, md, ebar, esig, mass(np)
45	format('param  ',i4,/,
     *         'run    ',i4,/,
     *         'mdstps ',i4,/,
     *         'ebar   ',1pg12.4,/,
     *         'esig   ',1pg12.4,/,
     *         'mass   ',1pg12.4,/)

	do 60 lag=0, nt-1
	    if (mod(lag,50) .eq. 0) print *, 'Doing lag ',lag
	    c = 0.
	    var = 0.
	    do 50 i=1, nt-lag
	        ind = i + lag
	        if (ind .gt. nt) pause 'Lag too large!'
	        dc = (vals(i)-mean) * (vals(ind)-mean)
	        c = c + dc
	        var = var + dc**2
50	    continue
	    rn = nt-lag
	    c = c / rn
	    var = var/rn - c**2
	    if (lag .eq. 0) then
	        c0 = c
	        var0 = var
	    endif
	    c = c / c0
	    var = ( var/c0**2 + var0*(c/c0)**2 ) / rn
c	    var = max(var,0.)
	    write(1,'(3(1pg13.4))') 1.*t(lag+1), c, sqrt(var)
60	continue
	close(1)
	end

