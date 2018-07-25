	program bint

c---	Bin samples at a particular simulation time from a
c	Metropolis time series of samples.
c
c	09 Jul 92  TJL

	integer maxp, maxs, maxb
	parameter (maxp = 10, maxs = 1000, maxb = 50)
	integer it, npar, ns, t, np, bins(maxb), nb, i, n
	real*8 pars(maxp),  par(maxs), lo, hi, bw
	logical rdpfnd
	character*60 ifile
	common /bintcb/ pars, lo, hi

	call opar(1,'bint.pf')
	call rdpw('input.file',ifile)
	call rdpi('time',it)
	call rdpi('param',np)
	call rdpda('bin.range',lo,2)
	call rdpi('bins',nb)
	call cpar

	call opar(1,ifile)
	call rdpi('npar',npar)
	if (np .gt. npar .or. npar .gt. maxp) pause 'Param # too high!'

c---  Read samples.
	ns = 0
20	    call rdpi('time',t)
	    if (rdpfnd()) then
	        if (t .eq. it) then
	            call rdpvec('samp',pars,npar)
	            ns = ns + 1
	            par(ns) = pars(np)
	        endif
	        goto 20
	    endif
	write(*,*) 'Samples found: ',ns
	call cpar

c---  Bin the data.
	do 100 i=1, nb
	    bins(i) = 0
100	continue
        bw = (hi - lo) / nb
        do 120 i=1, ns
            n = (par(i) - lo) / bw + 1
            bins(n) = bins(n) + 1
120     continue
	write(*,200) (bins(i),i=1,nb)
200	format('bins  ',10i6)

	end

