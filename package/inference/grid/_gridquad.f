c-----------------------------------------------------------------------
	function qgt1d (array, nx, crit)

c---	"Quadrature Greater Than - 1 Dimensional"
c
c	Return the quadrature (by the trapezoid rule) of the
c	contents of array above the value crit.
c
c	Step size factors (i.e., dx) are *not* included in the
c	result, so this subroutine needn't know these.

c+++  Arguments:
	integer nx
	real*8 qgt1d, array(nx), crit

cf2py intent(in) nx, array, crit

c+++  Locals:
	integer i
	real*8 lo, hi

	qgt1d = 0.
	do 70 i=1, nx-1
	    if (array(i) .gt. crit .and. array(i+1) .gt. crit) then
		qgt1d = qgt1d + 0.5 * (array(i) + array(i+1))
	    else if (array(i) .gt. crit .or.
     *               array(i+1) .gt. crit) then
		hi = max(array(i), array(i+1))
		lo = min(array(i), array(i+1))
		qgt1d = qgt1d + 0.5 * (hi + crit) *
     *                    (hi - crit) / (hi - lo)
	    endif
70	continue

	return
	end
	
c-----------------------------------------------------------------------
	function qgt2d (array, nx, ny, crit)

c---	"Quadrature Greater Than - 2 Dimensional"
c
c	Return the quadrature (by the trapezoid rule) of the
c	contents of array above the value crit.
c
c	Step size factors (i.e., dx*dy) are *not* included in the
c	result, so this subroutine needn't know these.

c+++  Arguments:
	integer nx, ny
	real*8 qgt2d, array(nx, ny), crit

cf2py intent(in) nx, ny, array, crit

c+++  Locals:
	integer i, j
	real*8 lo, hi, qy

	qgt2d = 0.
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
     *                   (hi -crit) / (hi - lo)
		endif
80	    continue
	    if (i .eq. 1 .or. i .eq. nx) then
	        qgt2d = qgt2d + 0.5 * qy
	    else
	        qgt2d = qgt2d + qy
	    endif
90	continue

	return
	end
	
c-----------------------------------------------------------------------
        subroutine xvalues(xarr, yarr, n, ycrit, xvals, nxvmax, nx, ok)

c---  	"X VALUES"
c
c	Return the values of x for which y=ycrit, linearly
c	interpolating the x and y arrays.

        integer n, nxvmax, nx, i
        real*8 xarr(n), yarr(n), xvals(nxvmax), ycrit
        logical ok

cf2py intent(in) n, nxv
cf2py intent(in) xarr, yarr, ycrit
cf2py intent(out) xvals, nx, ok

	ok = .true.
	nx = 0
	do 20 i=1, n-1
	    if (yarr(i) .lt. ycrit .and. yarr(i+1) .gt. ycrit) then
	        nx = nx + 1
	        xvals(nx) = xarr(i) + (xarr(i+1)-xarr(i)) *
     *            (ycrit - yarr(i)) / (yarr(i+1) - yarr(i))
	    else if (yarr(i+1) .lt. ycrit .and. yarr(i) .gt. ycrit) then
	        nx = nx + 1
	        if (nx .gt. nxvmax) then
	            ok = .false.
	            return
	        endif
	        xvals(nx) = xarr(i) + (xarr(i+1)-xarr(i)) *
     *            (ycrit - yarr(i)) / (yarr(i+1) - yarr(i))
	    endif
20	continue

        return
        end
