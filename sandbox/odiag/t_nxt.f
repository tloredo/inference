	program t_nxt

c---	"Test NXTarc"
c
c	Compare nxtarc with memcof, to make sure it gives the
c	same coefficients.

	integer ndat, t, cur, p, csiz, off, s, ord, seed, r, tt
	parameter (ndat = 1024, p = 8, csiz = (p*(p+1))/2)
	real*8 data(ndat), coefs(csiz), wrk(ndat), wrks(ndat)
c	real*8 d(p)
	real*8 wrk3(ndat), vars(p), cor(5), sig, xms, AIC, AIClo, bsig
	integer besto, nft, tlo, nwt, nu, del, ord
	real*8 wt(ndat), sum, mean, pow, smean, ln
	real*8 dmean,  arsig, spsig, clen
	logical ok
	real*8 gasdev

c---  Simulate correlated data.  Start far from the mean to put
c---  an initial transient in the series, and further add a linear
c---  transient.
	cor(1) = .4
	cor(2) = .2
	cor(3) = .1
	cor(4) = .05
	cor(5) = .05
	sig = 1.
	mean = 2.
	sum = ndat
	sum = (1.-(cor(1)+cor(2)+cor(3)+cor(4)+cor(5)))*sqrt(sum)
	print *, 'True sig = ', sig / sum
	call drinit(1,seed)
	do 10 t=1, 5
	    data(ndat-t+1) = 2.*t*sig*mean
10	continue
	do 40 tt=1, ndat
	    t = tt
	    if (t .gt. ndat) t = tt - ndat
	    data(t) = mean + sig*gasdev(seed)
	    do 20 s=1, 5
	        r = t - s
	        if (r .lt. 1) r = ndat - (s-t)
	        data(t) = data(t) + cor(s)*(data(r)-mean)
20	    continue
	    data(t) = data(t) + 2.*sig*max(100.-t,0.)/100.
40	continue

c---  Run diagnostics.
	call odiag(ndat, data, del, dmean, arsig, spsig, ord, clen,
     *                  wrk, wrks, wrk3, ok)
	if (.not. ok) pause 'odiag failed!'
	write(*,43) del, dmean, ord, clen, arsig, spsig
43	format('Stationary after time ',i6,' with mean ',1pg12.4,/,
     *         '  AR order, clen, sig: ',i4,2(1pg12.4),/,
     *         '  Spectral sig:    ',1pg12.4,/)

c---  Find the sample mean and subtract it.
	smean = 0.
	do 45 t=1, ndat
	    smean = smean + data(t)
45	continue
	smean = smean / ndat
	do 50 t=1, ndat
	    data(t) = data(t) - smean
50	continue

c---  Calculate AR coefs for a bunch of orders.
	do 60 ord=1, p
	    cur = ord - 1
	    call nxtarc (ndat, data, cur, vars, coefs, wrk, wrks)
60	continue

c---  Compare with MEMCOF.
	ln = ndat
	ln = log(ln)
	do 100 ord=1, p
c	    print *, 'Order ',ord,':'
	    cur = ord - 1
	    call nxtarc (ndat, data, cur, vars, coefs, wrk, wrks)
	    off = (ord*(ord-1))/2
c	    call memcof(data,ndat,ord,xms,d)
c	    write(*,120) xms, vars(ord), xms-vars(ord)
	    sum = 0.
	    do 80 s=1, ord
	        sum = sum + coefs(off+s)
c	        if (abs(coefs(off+s)-d(s)) .gt. 1.e-13)
c     *            write(*,140) s, d(s), coefs(off+s), d(s)-coefs(off+s)
80	    continue
	    sig = vars(ord) / (ndat*(1.-sum)**2)
	    sig = sqrt(sig)
c	    AIC = ndat*log(vars(ord)) + 2.*ord
	    AIC = ndat*log(vars(ord)) + ln*ord
	    write(*,150) ord, xms, sig, AIC
	    if (ord .eq. 1) then
	        besto = ord
	        AIClo = AIC
	        bsig = sig
	    else
	        if (AIC .lt. AIClo) then
	            AIClo = AIC
	            besto = ord
	            bsig = sig
	        endif
	    endif
100	continue
120	format('  Variance:  ',3(1pg12.4))
140	format('  Coef:  ',i4,3(1pg12.4))
150	format('  ord, xms, sig, AIC:  ',i4,3(1pg12.4))

c---  Print the best order.
	print *, 'Best order:  ',besto
	write(*,240) smean, bsig
240	format('Mean, sig: ',2(1pg12.4))
c	write(*,'(a,1pg12.4)') 'Sample mean: ',smean
	off = (besto*(besto-1))/2
	do 250 s=1, besto
	    write(*,140) s, coefs(off+s)
250	continue

c---  Use the H&W spectral method.
	do 255 t=1, ndat
	    data(t) = data(t) + smean
255	continue
	call specvar (ndat, data, sig, nu, wrk)
	sig = sqrt(sig)
	write(*,'(a,1pg12.4,i5)') 'Spectral sigma:  ',sig,nu
	if (sig .gt. 0.) then
	    call drend(1,seed)
	    stop
	endif

c---  Use the best order to "prewhiten" the time series.
	tlo = besto + 1
	nwt = ndat - tlo + 1
	nwt = ndat
	do 260 t=1, nwt
	    wt(t) = sig*gasdev(seed)
c	    wt(t) = 0.
c	    do 250 s=1, besto
c	        wt(t) = wt(t) + coefs(off+s)*data(tlo+t-1-s)
c250	    continue
c	    wt(t) = data(tlo+t-1) - wt(t)
260	continue

c---  Pad it with 0s.
	do 280 t=nwt+1, ndat
	    wt(t) = 0.
280	continue

c---  Do an FFT of the prewhitened residuals.
	nft = ndat
	call realft(wt, nft, 1)

c---  Write out power at 1st few frequencies.
	do 300 t=1, 10
	    if (t .eq. 1) then
	        pow = wt(1)**2
	    else
	        pow = wt(2*t-1)**2 + wt(2*t)**2
	    endif
	    write(*,320) t, pow
300	continue
320	format(' Freq: ',i4,'  power = ',1pg12.4)
	print *

c---  Use wt(1) in place of var(besto) in the sig calculation.
	sum = 0.
	do 340 s=1, p
	    sum = sum + coefs(off+s)
340	continue
	sig = abs(wt(1) / (1. - sum))
	print *, 'Prewhitening results: '
	write(*,240) mean, sig

	call drend(1,seed)

	end

