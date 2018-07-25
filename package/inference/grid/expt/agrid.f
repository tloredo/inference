	program agrid

c	"Analyze GRID"
c
c	Analyze a log posterior grid produced by grbd,
c	finding the global likelihood and values
c	defining the boundaries of credible regions.
c
c	We presume the prior is flat over the grid; ie, if
c	the grid is linear, the prior is constant; if logarithmic,
c	it is log-constant.
c
c	Tested with test_agrid.f for 1 and 2 grids, linear, log,
c	and mixed axes.
c
c	08 Sep 95  TJL

c+++  Parameters:
	integer nlev
	real*8 tol
	parameter (nlev = 3, tol = 0.001)
	real*8 level(nlev), rgauss(nlev,2)
	data level/0.6827, 0.9545, 0.9973/
	data rgauss/0.61, 0.14, 0.01, 0.32, 0.05, 0.003/

c+++  Globals:
	integer nxymax, ngmax
	parameter (nxymax = 100, ngmax = 2)
	integer ng, anx(ngmax), any(ngmax)
	real*8 pprob(nxymax,nxymax,ngmax), dx(ngmax), dy(ngmax)
	real*8 norm, cfrac
	common /gridcb/ pprob, dx, dy, norm, cfrac, ng, anx, any

c+++  Locals:
	character ifile(ngmax)*60, args*80, flag*60
	integer nx, ny, ndim, cloff
	real*8 xl, xh, yl, yh, xlo, xhi, ylo, yhi
	logical success, ok, xlog, ylog
	real*8 lpmax, arg, prior
	real*8 gnorm, normy, lgl
	real*8 lpc1, lpc2, lr1, lr2, lpcrit(nlev)
	integer i, j, n, nl, nd, nxv(nlev)
	real*8 xarr(nxymax), yarr(nxymax), xvals(nxymax,nlev), del
	logical use_fp

c+++  Functions:
	integer nclarg, nwlen
	logical rdpfnd
	real*8 lratio, zbrent2
	external lratio

c---  Check the command line for info about the prior.
	ng = nclarg()
	print *, ng
        if (ng .lt. 1) pause 'Not enough arguments!'
	if (ng .ge. 1) then
	    call clarg(1,flag)
	    if (flag(1:3) .eq. '-pl') then
	        ng = ng - 1
	        if (ng .lt. 1) pause 'Not enough arguments!'
	        cloff = 1
	        use_fp = .false.
	        print *, 'Treating grid as prior*likelihood.'
	    else
	        cloff = 0
	        use_fp = .true.
	        print *, 'Multiplying by flat, proper prior.'
	    endif
	endif
	write(*,*)

c---  Get the # of input file names.
        if (ng .gt. ngmax) pause 'Too many grids specified!'

c---  Read in the grids.
        do 5 n=1, ng
            call clarg(n+cloff,ifile(n))
            call opar(1,ifile(n))

c>>>  See if the grids are 1d or 2d.
	    call rdpi('num.stepped.syms',nd)
	    if (.not. rdpfnd()) pause 'Bad input file!'
	    if (n .eq. 1) then
	        ndim = nd
	    else
	        if (nd .ne. ndim) pause 'Incompatible ndim!'
	    endif

c>>>  Read in the grid data.
	    if (ndim .eq. 1) then
	        ok = .true.
	        call fndcmd(1,'step',args,ok)
	        call rddarg(args,xl,1,ok)
	        call rddarg(args,xh,1,ok)
	        call rdiarg(args,nx,1,ok)
	        if (.not. ok) pause 'Bad 1d step data!'
	        if (n .eq. 1) then
	            xlog = nx .lt. 0
	            xlo = xl
	            xhi = xh
	        else
	            if (xlog .and. (nx .gt. 0)) pause 'Incompat. nx!'
	            if (.not. xlog .and. (nx .lt. 0)) 
     *                  pause 'Incompat. nx!'
	            xlo = min(xl,xlo)
	            xhi = max(xh,xhi)
	        endif
	        call rdpcol ('#', '<', '>', 2, pprob(1,1,n), 
     *                anx(n), nxymax)
	        print *, 'Read ',anx(n),' 1-d posterior values.'
	        if (anx(n) .ne. abs(nx)) pause 'nx mismatch!'
	        ny = 1
	    else
	        call rdadat (pprob, nxymax, ngmax, n, nx, xl, xh,
     *                     ny, yl, yh)
	        if (n .eq. 1) then
	            xlog = nx .lt. 0
	            xlo = xl
	            xhi = xh
	            ylog = ny .lt. 0
	            ylo = yl
	            yhi = yh
	        else
	            if (xlog .and. (nx .gt. 0)) pause 'Incompat. nx!'
	            if (.not. xlog .and. (nx .lt. 0)) 
     *                  pause 'Incompat. nx!'
	            if (ylog .and. (ny .gt. 0)) pause 'Incompat. ny!'
	            if (.not. ylog .and. (ny .lt. 0)) 
     *                  pause 'Incompat. ny!'
	            xlo = min(xl,xlo)
	            xhi = max(xh,xhi)
	            ylo = min(yl,ylo)
	            yhi = max(yh,yhi)
	        endif
	    endif
	    call cpar
	    anx(n) = abs(nx)
	    any(n) = abs(ny)
	    if (xlog) then
	        dx(n) = log(xh/xl) / (anx(n) - 1.)
	    else
	        dx(n) = (xh - xl) / (anx(n) - 1.)
	    endif
	    if (ndim .eq. 2) then
	        if (ylog) then
	            dy(n) = log(yh/yl) / (any(n) - 1.)
	        else
	            dy(n) = (yh - yl) / (any(n) - 1.)
	        endif
	    else
	        dy(n) = 1.
	    endif
5	continue

c---  Find the prior factor.
	if (xlog) then
	    prior = 1. / log(xhi/xlo)
	else
	    prior = 1. / abs(xhi - xlo)
	endif
	if (ndim .eq. 2) then
	    if (ylog) then
	        prior = prior / log(yhi/ylo)
	    else
	        prior = prior / abs(yhi - ylo)
	    endif
	endif
	if (.not. use_fp) prior = 1.

c---  Find the maximum value of the log posterior.
	lpmax = pprob(1,1,1)
	do 25, n=1, ng
	    do 20 i=1, anx(n)
	        do 10 j=1, any(n)
	            if (pprob(i,j,n) .gt. lpmax) lpmax = pprob(i,j,n)
10	        continue
20	    continue
25	continue

c---  Now change the array to the posterior, scaled by the mode.
c---  Find the global likelihood along the way.
	norm = 0.
	if (ndim .eq. 1) then
	    do 35 n=1, ng
	        gnorm = 0.
	        do 30 i=1, anx(n)
	            arg = pprob(i,1,n) - lpmax
	            if (arg .gt. -200.) then
	               pprob(i,1,n) = exp(arg)
	            else
	                pprob(i,1,n) = 0.
	            endif
	            if (i .eq. 1 .or. i .eq. anx(n)) then
	                gnorm = gnorm + 0.5 * pprob(i,1,n)
	            else
	                gnorm = gnorm + pprob(i,1,n)
	            endif
30	        continue
	        norm = norm + gnorm*dx(n)
35	    continue
	else
	    do 55 n=1, ng
	        gnorm = 0.
	        do 50 i=1, anx(n)
	            normy = 0.
	            do 40 j=1, any(n)
	                arg = pprob(i,j,n) - lpmax
	                if (arg .gt. -200.) then
	                    pprob(i,j,n) = exp(arg)
	                else
	                    pprob(i,j,n) = 0.
	                endif
	                if (j .eq. 1 .or. j .eq. any(n)) then
	                    normy = normy + 0.5 * pprob(i,j,n)
	                else
	                    normy = normy + pprob(i,j,n)
	                endif
40	            continue
	            if (i .eq. 1 .or. i .eq. anx(n)) then
	                gnorm = gnorm + 0.5 * normy
	            else
	                gnorm = gnorm + normy
	            endif
50	        continue
	        norm = norm + gnorm*dx(n)*dy(n)
55	    continue
	endif
	lgl = log(prior*norm) + lpmax

c---  Loop thru levels for which we want credible regions.
	do 200 nl=1, nlev

c>>>  Bracket the root, starting with guesses for the bounding density c>>>  (ratio) based on the Gaussian distribution.  We work in terms
c>>>  of the log(ratio).
	    cfrac = level(nl)
	    lpc2 = log(rgauss(nl,ndim))
	    lpc1 = lpc2 - 0.05
	    call zbrac2 (lratio, lpc1, lpc2, lr1, lr2, success)
	    if (.not. success) then
	        print *, 'NLEV = ',nlev
	        pause 'Bracketing failed...'
	        stop
	    endif

c>>>  Find the critical level with Brent's method.
	    lpcrit(nl) = zbrent2 (lratio, lpc1, lpc2, lr1, lr2, 
     *                           tol, success)
	    if (.not. success) pause 'ZBRENT2 failed!'

c>>>  For 1d arrays, find the boundaries of the region.
	    if (ndim .eq. 1 .and. ng .eq. 1) then
	        if (xlog) then
	            del = log(xhi/xlo) / (anx(1)-1)
	        else
	            del = (xhi - xlo) / (anx(1)-1)
	        endif
	        do 100 i=1, anx(1)
	            if (xlog) then
	                xarr(i) = log(xlo) + (i-1)*del
	            else
	                xarr(i) = xlo + (i-1)*del
	            endif
	            yarr(i) = log(pprob(i,1,1))
100	        continue
	        call xvalues (xarr, yarr, anx(1), 
     *                      lpcrit(nl), xvals(1,nl), nxv(nl))
	        if (xlog) then
	            do 120 i=1, nxv(nl)
	                xvals(i,nl) = exp(xvals(i,nl))
120	            continue
	        endif
	    endif

200	continue

c---  Write it all to stdout.
	write(*,'(/,a)') 'AGRID results for following grid files:'
	do 300 n=1, ng
	    write(*,305) ifile(n)(1:nwlen(ifile)),anx(n),any(n)
300	continue
305	format('gfile  ',a,2i5)
	write(*,*)
	write(*,310) ndim, lpmax
310	format('ndim  ',i5,/,
     *         'lpmax   ',1pg14.6,/)
	write(*,320)'x.range  ', xlo, xhi, 'x.log  ', xlog
320	format(a, 2(1pg12.4),/,
     *         a, l)
	if (ndim .gt. 1) write(*,320) 'y.range  ', ylo, yhi, 
     *          'y.log  ', ylog
	write(*,*)
	write(*, '(a,1pg14.6,/)') 'log.gl  ', lgl
	do 330 nl=1, nlev
	    write(*,'(a,2(1pg12.4))') 'lev.lr  ', level(nl), lpcrit(nl)
	    if (ndim .eq. 1 .and. ng .eq. 1) then
	        write(*, '(6hxvals ,4(1pg12.4))') 
     *             (xvals(i,nl),i=1,nxv(nl))
	    endif
330	continue

	end

c-----------------------------------------------------------------------
        subroutine rdadat(array, ndim, ngmax, ng,
     *        nx, xlo, xhi, ny, ylo, yhi)

c---	"ReaD Array DATa"

c+++  Arguments:
	integer ndim, ngmax, ng, nx, ny
	real*8 array(ndim,ndim,ngmax), xlo, xhi, ylo, yhi

c+++  Locals:
	real*8 a, b
	common /rangcb/ a, b
        integer i
        logical rdpfnd

c---  Read in the ranges and steps associated with the array.
        call rdpi('x.points',nx)
        call rdpda('x.range',a,2)
        xlo = a
        xhi = b
        call rdpi('y.points',ny)
        call rdpda('y.range',a,2)
        ylo = a
        yhi = b

c---  Read in the array itself.
        if (abs(nx) .gt. ndim .or. abs(ny) .gt. ndim)
     +    pause 'Array too large!'
        do 20 i=1, abs(ny)
            call rdpvec('row', array(1,i,ng), abs(nx))
            if (.not. rdpfnd()) pause 'RDADAT:  Could not find row!'
20      continue
        write(*,'(a,i3,a,i3,a)') 
     +        ' Read data for ',abs(nx),' x ',abs(ny),' array.'

        return
        end

c-----------------------------------------------------------------------
	function qgt (array, nxymax, nx, ny, crit)

c---	"Quadrature Greater Than"
c
c	Return the quadrature (by the trapezoid rule) of the
c	contents of array above the value crit.
c
c	If ny=1, a 1-d quadrature rule is used (of course!).
c
c	Step size factors (e.g., dx*dy) are *not* included in the
c	result, so this subroutine needn't know these.

c+++  Arguments:
	integer nxymax, nx, ny
	real*8 qgt, array(nxymax,nxymax), crit

c+++  Locals:
	integer i, j
	real*8 lo, hi, qy

	qgt = 0.
	if (ny .eq. 1) then
	    do 70 i=1, nx-1
	        if (array(i,1) .gt. crit .and. 
     *              array(i+1,1) .gt. crit) then
		    qgt = qgt + 0.5 * (array(i,1) + array(i+1,1))
		else if (array(i,1) .gt. crit .or.
     *                   array(i+1,1) .gt. crit) then
		    hi = max(array(i,1), array(i+1,1))
		    lo = min(array(i,1), array(i+1,1))
		    qgt = qgt + 0.5 * (hi + crit) *
     *                     (hi - crit) / (hi - lo)
		endif
70	    continue
	else
	    do 90 i=1, nx
	        qy = 0.
	        do 80 j=1, ny-1
	            if (array(i,j) .gt. crit .and. 
     *                  array(i,j+1) .gt. crit) then
		        qy = qy + 0.5*(array(i,j) + array(i,j+1))
		    else if (array(i,j) .gt. crit .or.
     *                       array(i,j+1) .gt. crit) then
		        hi = max(array(i,j), array(i,j+1))
		        lo = min(array(i,j), array(i,j+1))
		        qy = qy + 0.5 * (hi +crit) *
     *                         (hi -crit) / (hi - lo)
		    endif
80	        continue
	        if (i .eq. 1 .or. i .eq. nx) then
	            qgt = qgt + 0.5 * qy
	        else
	            qgt = qgt + qy
	        endif
90	    continue
	endif

	return
	end
	
c-----------------------------------------------------------------------
	function lratio (lpcrit)

c---	"Log RATIO"
c
c	Logarithm of the ratio of the fraction of the posterior
c	with log(density) > lcrit to the desired fraction.

c+++  Arguments:
	real*8 lratio, lpcrit

c+++  Globals:
	integer nxymax, ngmax
	parameter (nxymax = 100, ngmax = 2)
	integer ng, anx(ngmax), any(ngmax)
	real*8 pprob(nxymax,nxymax,ngmax), dx(ngmax), dy(ngmax)
	real*8 norm, cfrac
	common /gridcb/ pprob, dx, dy, norm, cfrac, ng, anx, any

c+++  Locals:
	real*8 pcrit, tot
	integer n

c+++  Functions:
	real*8 qgt

	pcrit = exp(lpcrit)
	tot = 0.
	do 10 n=1, ng
	    tot = tot + dx(n)*dy(n)*
     *            qgt(pprob(1,1,n), nxymax, anx(n), any(n), pcrit)
10	continue
	if (tot .le. 0.) then
	    lratio = -500.
	else
	    lratio = tot / norm
	    lratio = log(lratio/cfrac)
	endif

	return
	end

c-----------------------------------------------------------------------
        subroutine xvalues(xarr, yarr, n, ycrit, xvals, nx)

c---  	"X VALUES"
c
c	Return the values of x for which y=ycrit, linearly
c	interpolating the x and y arrays.

        real*8 xarr(*), yarr(*), xvals(*), ycrit
        integer n, nx, i

	nx = 0
	do 20 i=1, n-1
	    if (yarr(i) .lt. ycrit .and. yarr(i+1) .gt. ycrit) then
	        nx = nx + 1
	        xvals(nx) = xarr(i) + (xarr(i+1)-xarr(i)) *
     *            (ycrit - yarr(i)) / (yarr(i+1) - yarr(i))
	    else if (yarr(i+1) .lt. ycrit .and. yarr(i) .gt. ycrit) then
	        nx = nx + 1
	        xvals(nx) = xarr(i) + (xarr(i+1)-xarr(i)) *
     *            (ycrit - yarr(i)) / (yarr(i+1) - yarr(i))
	    endif
20	continue

        return
        end
