	subroutine mk2dc(n, array)
	integer n, i, j
	real*8 array(n,n)

cf2py  intent(in) n
cf2py  intent(c,out) array

	do 20 i=1, n
	    do 10 j=1, n
	        array(i,j) = 10*i + j
10	    continue
20	continue

	return
	end

	subroutine mk2df(n, array)
	integer n, i, j
	real*8 array(n,n)

cf2py  intent(in) n
cf2py  intent(out) array

	do 20 i=1, n
	    do 10 j=1, n
	        array(i,j) = 10*i + j
10	    continue
20	continue

	return
	end

c-----------------------------------------------------------------------
        subroutine ninfol(n, mu, lfmax, dp, iter, info, ierr, 
     *                    h, d2f)


c	Skeleton to search for memory problem.
c
c	PROBLEM:  Dimensions of h, d2f should be ITER!!!!

c+++  Arguments:
	integer n, iter
	real*8 mu(n), lfmax, dp(n)
	real*8 info(n,n), ierr(n,n)
	real*8 h(n), d2f(n)

cf2py  intent(in) lf, n, mu, lfmax, dp, iter
cf2py  intent(out) info, ierr
cf2py  intent(in,out) h, d2f

c+++  Locals:
        real*8 ratio, one
        parameter (ratio = 0.707, one = 1.d0)

        real*8 mu1, mu2, di, extrap
        integer i, j, k

c---  First find the derivatives along the diagonal (unmixed).
        do 100 i=1, n
            mu1 = mu(i)

c>>>  Begin by finding iter second differences.
            do 20 j=1, iter
                h(j) = dp(i) * ratio**(iter-j)
                mu(i) = mu1 + h(j)
                d2f(j) = j
                mu(i) = mu1 - h(j)
                d2f(j) = (d2f(j) + j+.5 - 2.*lfmax)/h(j)**2
20          continue

c>>>  Now extrapolate them to h = 0.
	    if (iter .gt. 1) then
                extrap = 5.
                di = .01
                info(i,i) = -extrap
                ierr(i,i) = abs(di)
            else
                info(i,i) = -d2f(1)
                ierr(i,i) = 0.
            endif

c>>>  End loop over the diagonal.
            mu(i) = mu1
100     continue

c---  Now find the mixed derivatives
        do 200 i=1, n-1
            mu1 = mu(i)
            do 180 j=i+1, n
                mu2 = mu(j)

c>>>  Find iter second differences of f in terms of a
c>>>  new "diagonal" variable.
                do 120 k=1, iter
                    h(k) = ratio**(iter-k)
                    mu(i) = mu1 + h(k)/sqrt(info(i,i))
                    mu(j) = mu2 + h(k)/sqrt(info(j,j))
                    d2f(k) = k
                    mu(i) = mu1 - h(k)/sqrt(info(i,i))
                    mu(j) = mu2 - h(k)/sqrt(info(j,j))
                    d2f(k) = (d2f(k) + k+.5 - 2.*lfmax)/h(k)**2
120             continue

c>>>  Now extrapolate them to h = 0, and eliminate the
c>>>  effects of the transformation.
	        if (iter .gt. 1) then
                    extrap = 10.
                    di = .01
                    info(i,j) = -(extrap + 2.) * 
     *                           sqrt(info(i,i)*info(j,j)) / 2.
	            ierr(i,j) = abs(di)*sqrt(info(i,i)*info(j,j)) / 2.
                    info(j,i) = info(i,j)
                    ierr(j,i) = ierr(i,j)
                else
                    info(i,j) = -(d2f(1) + 2.) * 
     *                           sqrt(info(i,i)*info(j,j)) / 2.
	            ierr(i,j) = 0.
                    info(j,i) = info(i,j)
                    ierr(j,i) = 0.
                endif

c>>>  End the loops.
                mu(j) = mu2
180         continue
            mu(i) = mu1
200     continue

        return
        end

