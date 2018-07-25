	program pickt

c---	Pick samples at a particular simulation time from a
c	Metropolis time series of samples.
c
c	01 May 92  TJL

	integer maxp
	parameter (maxp = 10)
	integer it, npar, ns, t, nruns
	real*8 pars(maxp), p
	logical rdpfnd
	character*60 ifile
	common /pickcb/ pars

	print *, 'Enter file to read:'
	read(*,'(a)') ifile
	print *, 'Enter time to select:'
	read(*,*) it
	call opar(1,ifile)
	call rdpi('npar',npar)
	if (npar .gt. maxp) pause 'Param # too high!'
	call rdpi('nruns',nruns)

	call openf(2,'pickt.dat')
	write(2,10) it, npar, nruns
10	format('time  ',i5,/,'npar  ',i5,/,'nmc  ',i5,/)
	ns = 0
20	    call rdpi('time',t)
	    if (rdpfnd()) then
	        if (t .eq. it) then
	            call rdpvec('samp',pars,npar)
	            call rdpd('log.p',p)
	            p = exp(p)
                    call wrrows ('samp ', pars, npar, 5, 2)
                    write(2,'(a,2(1pg14.6))') 'p.pj  ',p,1.
	            ns = ns + 1
	        endif
	        goto 20
	    endif

	write(*,*) 'Samples found: ',ns
	call cpar
	close(2)
	end

