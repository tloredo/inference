	program rep4

c	Calculate odds for repetition for a selection of up to 4
c	bursts, exactly (such as the 4 consecutive repeaters
c	of October 1996).
c
c	11 Dec 96

c+++  Globals:
	integer nbmax
	parameter (nbmax = 4)
	real*8 drxn(3,nbmax), kappa(nbmax)
	integer nb
	common /burstcb/ drxn, kappa, nb

c+++  Locals:
	real*8 pi, onedeg, pi4
	parameter (pi = 3.141592653589793d0, onedeg = pi/180.)
	parameter (pi4 = 4*pi)
	real*8 long, lat, sig, cd, sys
	common /parcb/ long, lat, sig
	integer i
	character*60 ofile

c---  Read in burst directions and uncertainties.
	call opar(1,'rep4.pf')
	call rdpw('output.file',ofile)
	call rdpd('sys.err', sys)
	call rdpi('#.bursts',nb)
	if (nb .gt. nbmax) pause 'Too many bursts!'
	sys = sys*onedeg
	do 20 i=1, nb
	    call rdpda('burst',long,3)
	    long = long*onedeg
	    lat = lat*onedeg
	    sig = sqrt((sig*onedeg)**2 + sys**2)
	    cd = cos(lat)
	    drxn(1,i) = cos(long) * cd
	    drxn(2,i) = sin(long) * cd
	    drxn(3,i) = sin(lat)
	    call sig_to_kappa (sig, kappa(i))
20	continue
	call cpar

c---  Write out results on the fly...
	call openf(1, ofile)
	write(1,'(a,i4)') '#.bursts  ',nb
	write(1,*)
	if (nb .eq. 3) call rep_3
	if (nb .eq. 4) call rep_4

	close(1)
	end

c-----------------------------------------------------------------------
	subroutine rep_3

c	Calculate repeating odds for partitions of 3 bursts.

c+++  Globals:
	integer nbmax
	parameter (nbmax = 4)
	real*8 drxn(3,nbmax), kappa(nbmax)
	integer nb
	common /burstcb/ drxn, kappa, nb

c+++  Locals:
	real*8 singlet, doublet(nbmax,nbmax), triplet
	real*8 lk111, lk12, lk3
	integer i, j, n1, n2, n3

c+++  Functions:
	real*8 doublet_l, triplet_l

c---  Compile the doublet factors.
	n2 = 0
	do 220 i=1, nb-1
	    do 210 j=i+1, nb
	        n2 = n2 + 1
	        doublet(i,j) = doublet_l(i, j)
	        print *, i,j, doublet(i,j)
210	    continue
220	continue

c---  The rest are simple...
	n1 = 1
	singlet = 1.
	n3 = 1
	triplet = triplet_l(1,2,3)

c---  The nonrepeating likelihood is pretty simple!  (We renormalized
c---  to avoid 1/4*pi factors.)
	lk111 = singlet**3
	write(1,'(a, 1pg12.4)') '1,2,3     ', lk111
	write(1,*)

c---  Average the 1+2 partitions.  The singlet factors are unity.
	lk12 = 0
	do 410 i=1, nb-1
	    do 400 j=i, nb
	    lk12 = lk12 + doublet(i,j)
400	continue
410	continue
	lk12 = lk12 / 3
	write(1,'(a,1pg12.4)') '1,(2,3)   ', doublet(2,3)
	write(1,'(a,1pg12.4)') '2,(1,3)   ', doublet(1,3)
	write(1,'(a,1pg12.4)') '3,(1,2)   ', doublet(1,2)
	write(1,'(a,1pg12.4)') 'O(1+2) = ',lk12
	write(1,*)

c---  The 3 term is hard!  NOT!
	lk3 = triplet
	write(1,'(a,1pg12.4)') '(1,2,3)   ', triplet
	write(1,'(a,1pg12.4)') 'O(3) =   ',lk3
	write(1,*)

	return
	end

c-----------------------------------------------------------------------
	subroutine rep_4

c	Calculate repeating odds for partitions of 4 bursts.

c+++  Globals:
	integer nbmax
	parameter (nbmax = 4)
	real*8 drxn(3,nbmax), kappa(nbmax)
	integer nb
	common /burstcb/ drxn, kappa, nb

c+++  Locals:
	real*8 singlet, doublet(nbmax,nbmax)
	real*8 triplet(nbmax,nbmax,nbmax), quadruplet
	real*8 lk1111, lk112, lk22, lk13, lk4
	integer i, j, k, n1, n2, n3, n4

c+++  Functions:
	real*8 doublet_l, triplet_l, quadruplet_l

c---  Compile the doublet factors.
	n2 = 0
	do 220 i=1, nb-1
	    do 210 j=i+1, nb
	        n2 = n2 + 1
	        doublet(i,j) = doublet_l(i, j)
	        print *, i,j, doublet(i,j)
210	    continue
220	continue
	print *, n2, ' doublets.'

c---  Compile the triplet factors.
	n3 = 0
	do 330 i=1, nb-2
	    do 320 j=i+1, nb-1
	        do 310 k=j+1, nb
	            n3 = n3 + 1
	            triplet(i,j,k) = triplet_l(i, j, k)
	        print *, i,j,k, triplet(i,j,k)
310	        continue
320	    continue
330	continue
	print *, n3, ' triplets.'

c---  The rest are simple...
	n1 = 1
	singlet = 1.
	n4 = 1
	quadruplet = quadruplet_l(1,2,3,4)

c---  The nonrepeating likelihood is pretty simple!  (We renormalized
c---  to avoid 1/4*pi factors.)
	lk1111 = singlet**4
	write(1,'(a, 1pg12.4)') '1,2,3,4     ', singlet**4
	write(1,*)

c---  Average the 1+1+2 partitions.  The singlet factors are unity.
	lk112 = 0
	do 410 i=1, nb-1
	    do 400 j=i, nb
	    lk112 = lk112 + doublet(i,j)
400	continue
410	continue
	lk112 = lk112 / 6
	write(1,'(a,1pg12.4)') '1,2,(3,4)   ', doublet(3,4)
	write(1,'(a,1pg12.4)') '1,3,(2,4)   ', doublet(2,4)
	write(1,'(a,1pg12.4)') '1,4,(2,3)   ', doublet(2,3)
	write(1,'(a,1pg12.4)') '2,3,(1,4)   ', doublet(1,4)
	write(1,'(a,1pg12.4)') '2,4,(1,3)   ', doublet(1,3)
	write(1,'(a,1pg12.4)') '3,4,(1,2)   ', doublet(1,2)
	write(1,'(a,1pg12.4)') 'O(1+1+2) = ',lk112
	write(1,*)

c---  Average the 2+2 partitions.
	lk22 = ( doublet(1,2)*doublet(3,4) + doublet(1,3)*doublet(2,4) +
     +           doublet(1,4)*doublet(2,3) ) / 3.
	write(1,'(a,1pg12.4)') '(1,2),(3,4) ',doublet(1,2)*doublet(3,4)
	write(1,'(a,1pg12.4)') '(1,3),(2,4) ',doublet(1,3)*doublet(2,4)
	write(1,'(a,1pg12.4)') '(1,4),(2,3) ',doublet(1,4)*doublet(2,3)
	write(1,'(a,1pg12.4)') 'O(2+2) =   ',lk22
	write(1,*)

c---  Average the 1+3 partitions.  The singlet factors are unity.
	lk13 = ( triplet(2,3,4) + triplet(1,3,4) + triplet(1,2,4) +
     +           triplet(1,2,3) ) / 4.
	write(1,'(a,1pg12.4)') '1,(2,3,4)   ', triplet(2,3,4)
	write(1,'(a,1pg12.4)') '2,(1,3,4)   ', triplet(1,3,4)
	write(1,'(a,1pg12.4)') '3,(1,2,4)   ', triplet(1,2,4)
	write(1,'(a,1pg12.4)') '4,(1,2,3)   ', triplet(1,2,3)
	write(1,'(a,1pg12.4)') 'O(1+3) =   ',lk13
	write(1,*)

c---  The 4 term is hard!  NOT!
	lk4 = quadruplet
	write(1,'(a,1pg12.4)') '(1,2,3,4)   ', quadruplet
	write(1,'(a,1pg12.4)') 'O(4) =     ',lk4
	write(1,*)

	return
	end

c-----------------------------------------------------------------------
	subroutine sig_to_kappa (sig, kappa)

c---	"SIGMA TO KAPPA for fisher dist'n"
c
c	Sets the kappa parameter defining the Fisher dist'n for
c	the direction uncertainty for a burst.  sig is the 68.3%
c	confidence region angular radius in radians.

c+++  Arguments:
	real*8 sig, kappa

	real*8 zero, twopi
	parameter (zero = 0., twopi = 6.283185307179586d0)

c+++  Locals:
	real*8 one, minus, C, C1, eps
	integer nrmax
	parameter (one = 1., minus = -1., C = 0.683, C1 = 1.-C)
	parameter (eps = 1.e-5, nrmax = 25)
	integer iter
	real*8 kold, f, csig, arg, ekc, ekk, num, den

c---  Find kappa burst using Newton-Raphson.  Use
c---  the Gaussian (small angle) approximation as a starting guess.
	kold = 2.297707 / sig**2
	csig = cos(sig)
	iter = 0
10	continue
	if (iter .ge. nrmax) pause 'Too many iters for kappa!'
	arg = csig - 1.
	ekc = exp(kold*arg)
	if (kold .lt. 250.) then
	    ekk = exp(-2.*kold)
	else
	    ekk = 0.
	endif
	den = 1. - ekk
	num = 1. - ekc
	f = num / den - C
	kappa = kold -
     *          f * den / (2.*num*ekk/den - arg*ekc)
	iter = iter + 1
	if (abs(kappa-kold)/kappa .gt. eps) then
	    kold = kappa
	    goto 10
	endif

	return
	end

c-----------------------------------------------------------------------
	function doublet_l (i, j)

c---	"DOUBLET Likelihood"
c
c	Calculate the likelihood factor for a doublet.

c+++  Arguments:
	real*8 doublet_l
	integer i, j

c+++  Globals:
	integer nbmax
	parameter (nbmax = 4)
	real*8 drxn(3,nbmax), kappa(nbmax)
	integer nb
	common /burstcb/ drxn, kappa, nb

c+++  Locals:
	real*8 pi4
	parameter (pi4 = 4.*3.141592653589793d0)
	real*8 tot(3), R

c+++  Functions:
	real*8 lsinhc

	tot(1) = kappa(i)*drxn(1,i) + kappa(j)*drxn(1,j)
	tot(2) = kappa(i)*drxn(2,i) + kappa(j)*drxn(2,j)
	tot(3) = kappa(i)*drxn(3,i) + kappa(j)*drxn(3,j)
	R = sqrt(tot(1)**2 + tot(2)**2 + tot(3)**2)


	doublet_l = lsinhc(R) - lsinhc(kappa(i)) - lsinhc(kappa(j))
	doublet_l = exp(doublet_l)

	return
	end

c-----------------------------------------------------------------------
	function triplet_l (i, j, k)

c---	"TRIPLET Likelihood"
c
c	Calculate the likelihood factor for a triplet.

c+++  Arguments:
	real*8 triplet_l
	integer i, j, k

c+++  Globals:
	integer nbmax
	parameter (nbmax = 4)
	real*8 drxn(3,nbmax), kappa(nbmax)
	integer nb
	common /burstcb/ drxn, kappa, nb

c+++  Locals:
	real*8 pi4
	parameter (pi4 = 4.*3.141592653589793d0)
	real*8 tot(3), R
	integer n, indx(3)

c+++  Functions:
	real*8 lsinhc

	indx(1) = i
	indx(2) = j
	indx(3) = k

	tot(1) = 0.
	tot(2) = 0.
	tot(3) = 0.
	do 20 n=1, 3
	    tot(1) = tot(1) + kappa(indx(n))*drxn(1,indx(n))
	    tot(2) = tot(2) + kappa(indx(n))*drxn(2,indx(n))
	    tot(3) = tot(3) + kappa(indx(n))*drxn(3,indx(n))
20	continue
	R = sqrt(tot(1)**2 + tot(2)**2 + tot(3)**2)

	triplet_l = lsinhc(R)
	do 40 n=1, 3
	    triplet_l = triplet_l - lsinhc(kappa(indx(n)))
40	continue
	triplet_l = exp(triplet_l)

	return
	end

c-----------------------------------------------------------------------
	function quadruplet_l (i, j, k, l)

c---	"QUADRUPLET Likelihood"
c
c	Calculate the likelihood factor for a quadruplet.

c+++  Parameters:
	real*8 quadruplet_l
	integer i, j, k, l

c+++  Globals:
	integer nbmax
	parameter (nbmax = 4)
	real*8 drxn(3,nbmax), kappa(nbmax)
	integer nb
	common /burstcb/ drxn, kappa, nb

c+++  Locals:
	real*8 pi4
	parameter (pi4 = 4.*3.141592653589793d0)
	real*8 tot(3), R
	integer n, indx(4)

c+++  Functions:
	real*8 lsinhc

	indx(1) = i
	indx(2) = j
	indx(3) = k
	indx(4) = l

	tot(1) = 0.
	tot(2) = 0.
	tot(3) = 0.
	do 20 n=1, 4
	    tot(1) = tot(1) + kappa(indx(n))*drxn(1,indx(n))
	    tot(2) = tot(2) + kappa(indx(n))*drxn(2,indx(n))
	    tot(3) = tot(3) + kappa(indx(n))*drxn(3,indx(n))
20	continue
	R = sqrt(tot(1)**2 + tot(2)**2 + tot(3)**2)

	quadruplet_l = lsinhc(R)
	do 40 n=1, 4
	    quadruplet_l = quadruplet_l - lsinhc(kappa(indx(n)))
40	continue
	quadruplet_l = exp(quadruplet_l)

	return
	end

c-----------------------------------------------------------------------
	function lsinhc (x)

c	Compute log(sinh(x)/x).

c+++  Arguments:
	real*8 lsinhc, x

c+++  Locals:
	real*8 log2
	parameter (log2 = 0.6931471805599453)

	lsinhc = x - log2 - log(x)
	if (x .lt. 20.) lsinhc = lsinhc - exp(-2.*x)

	return
	end
