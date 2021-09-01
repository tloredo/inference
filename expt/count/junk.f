	real*8 s, t, tb, sf, term, c(0:100), d, nn
	integer k, n, nu, km1, km2

	write(*,'(a,:)') 'Enter s: '
	read(*,*) s
	n = 100
	t = 1.
	nu =50
	tb = 1.
	call setc(n,t,nu,tb,c,100)
	sf = s*(t+tb)
	nn = n + nu + sf - 1.
	d = nn*nn - 4.*n*sf + 4.*(n+nu)
	if (d .lt. 0.) then
	    km1 = n
	    km2 = 0
	else
	    km1 = 0.5 * (nn - sqrt(d)) + 1.
	    km2 = 0.5 * (nn + sqrt(d)) + 1.
	endif
c	if (4.*sf*(n+nu) .gt. (n+sf)**2) then
c	    km1 = n
c	    km2 = 0
c	else
c	    km1 = 0.5 * (n + sf - sqrt((n+sf)**2-4*sf*(n+nu)))
c	    km2 = 0.5 * (n + sf + sqrt((n+sf)**2-4*sf*(n+nu)))
c	endif
	print *, 'km = ',km1, km2
	term = 1.
	write(*,'(i5,2(1pg12.4))') 0, term, c(0)
	do 20  k=1, n
c	    term = term * sf *(n+nu-k) / ((n-k)*k)
	    term = term * sf *(n-k+1.) / ((n+nu-k+1.)*k)
	    write(*,'(i5,2(1pg12.4))') k, term, c(k)
20	continue
c	term = term * sf * nu / n
c	write(*,'(i5,2(1pg12.4))') n, term, c(n)
	end

c-----------------------------------------------------------------------
        subroutine setc(n_on, t_on, n_off, t_off, c, cdim)

c---    Initialize the c array of coefficients in the Poisson series
c       for the density of the signal.  The dimension of c(0:cdim)
c	must be >= n_on + 1, and the first element should be c(0).
c
c	c(i) is the probability that i of the n_on counts are from
c	the background.
c
c	The coefficients are calculated recursively, going up and
c	down from the largest one so roundoff error accumulates only
c	on the negligible terms.

        integer n_on, n_off, cdim, i, mode
        real*8 t_on, t_off, c(0:cdim), tfac, sum, n

        tfac = 1. + t_off/t_on
        n = n_on + n_off
        mode = n_on - t_on*n_off/t_off + 1.
        mode = max(mode,0)
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

