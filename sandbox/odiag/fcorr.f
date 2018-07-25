	program fcorr

c---	Calculate the autocorrelation function of a single run 
c	of Metropolis samples, using FFTs.
c
c	6 Jul 92  TJL

	integer maxt, maxp, work
	parameter (maxt = 2048, maxp = 15, work = 2*maxt)
	integer np, npar, nruns, lrun, intvl, nt, i, lag, len
	integer t(maxt)
	real*8 pars(maxp), vals(maxt), mean, msqr, c, sig
	real*8 ans(work)

	print *, 'Enter param #:'
	read(*,*) np
	call opar(1,'hmc.dat')
	call rdpi('npar',npar)
	if (np .gt. npar .or. npar .gt. maxp) pause 'Param # too high!'
	call rdpi('nruns',nruns)
	if (nruns .gt. 1) print *, 'Using 1 run of',nruns
	call rdpi('lrun',lrun)
	call rdpi('intvl',intvl)
	nt = lrun / intvl + 1
	if (nt .gt. maxt/2) pause 'Too many times!'
	if ((nt-1)*intvl .ne. lrun) then
	    print *, 'Skipping final sample.'
	endif
	write(*,'(a,i5,a)') 'Using ',nt,' samples.'
	if (nt .le. 512) then
	    len = 512
	else if (nt .le. 1024) then
	    len = 1024
	endif

c---  Read in the evenly spaced samples.
	do 20 i=1, nt
	    call rdpi('time',t(i))
	    call rdpvec('samp',pars,npar)
	    vals(i) = pars(np)
c	    vals(i) = sin(6.283*t(i)/100.)
20	continue
	call cpar

c---  Zero-pad the vals array.
	do 30 i=nt+1, len
	    vals(i) = 0.
30	continue

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

c---  Subtract off mean, divide by sig.
	do 50 i=1, nt
	    vals(i) = (vals(i) - mean) / sig
50	continue

c---  Cross correlate vals with itself.
	call correl(vals, vals, len, ans)
	
c---  Write it out.
	call openf(1,'fcorr.dat')
	do 60 lag=0, nt-1
	    c = ans(lag+1) / ans(1)
	    write(1,'(3(1pg13.4))') 1.*t(lag+1), c
60	continue
	close(1)
	end

