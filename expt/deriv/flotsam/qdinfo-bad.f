c-----------------------------------------------------------------------
	subroutine qdinfol0(functn, np, p, func, g, h, ndim, simp, info)

c---	"Quick and Dirty INFOrmation matrix from Log likelihood"
c
c	Use second differences to estimate the info matrix.
c
c	The simplex expansion idea (to overcome roundoff error) is
c	adapted from the "quadratic surface fit" termination step
c	of a Nelder-Mead simplex implementation archived
c	at StatLib, by D. E. Shaw, R. W. M. Wedderburn, and A. Miller.
c	Their surface fit does little more than take second differences,
c	as done here.
c
c	functn		real*8 function functn (point)
c	nop		# of parameters
c	p(i)		The current min loc'n (presum. in simplex)
c	func		Function value at maximum
c	g(i,j)		Coordinate j of simplex point i (nop+1 pts)
c	h(i)		values on the simplex
c	ndim		Physical dimension for g, h
c	simp		Criterion for expanding simplex to overcome
c			rounding errors before fitting; set to at
c			least 1e3*rounding errors.
c	info(i,j)	returns the information matrix
c
c	Created:  May 97  Tom Loredo

c+++  Arguments:
	integer np, neval, ndim
	real*8 p(ndim), g(ndim+1,ndim), h(ndim+1), simp, info(ndim,ndim)
	real*8 functn
	external functn

c+++  Locals:
	integer npmax, n2max
	parameter (npmax = 20, n2max = npmax*(npmax+1)/2)
	real*8 zero, half, one, two, three
	parameter (zero = 0.d0, half = 0.5d0)
	parameter (one = 1.d0, two = 2.d0, three = 3.d0)
	integer i, j, np1, i1, i2, j1, k
	real*8 pstst(npmax), pstar(npmax), aval(npmax)
	real*8 func
	real*8 test, a0, hstst
	real*8 d1, d2

c---  Check array size constraints.
	if (np .gt. npmax) pause 'QSFIT:  Too many params!'

c---  Assume all parameters vary.
	np1 = np + 1

C
C     EXPAND THE SIMPLEX, IF NECESSARY, TO OVERCOME ROUNDING
C     ERRORS.
C
      NEVAL=0
      DO 490 I=1,NP1
  470   TEST=ABS(H(I)-FUNC)
        IF(TEST.GE.SIMP) GO TO 490
        DO 480 J=1,NP
          G(I,J)=(G(I,J)-P(J))+G(I,J)
          PSTST(J)=G(I,J)
  480   CONTINUE
        H(I) = FUNCTN(PSTST)
        NEVAL=NEVAL+1
        GO TO 470
  490 CONTINUE
C
C     FUNCTION VALUES ARE CALCULATED AT AN ADDITIONAL NAP POINTS.
C
      DO 510 I=1,NP
        I1=I+1
        DO 500 J=1,NP
  500   PSTAR(J)=(G(1,J)+G(I1,J))*HALF
        AVAL(I) = FUNCTN(PSTAR)
        NEVAL=NEVAL+1
  510 CONTINUE
C
C     THE MATRIX OF ESTIMATED SECOND DERIVATIVES IS CALCULATED AND ITS
C     LOWER TRIANGLE STORED IN BMAT.
C
      A0=H(1)
      DO 540 I=1,NP
        I1=I-1
        I2=I+1
        IF(I1.LT.1) GO TO 540
            d1 = g(i+1,i) - g(1,i)
        DO 530 J=1,I1
          J1=J+1
                d2 = g(j+1,j) - g(1,j)
          DO 520 K=1,NP
  520     PSTST(K)=(G(I2,K)+G(J1,K))*HALF
          HSTST = FUNCTN(PSTST)
          NEVAL=NEVAL+1
          info(i,j)= -4.d0*(HSTST+A0-AVAL(I)-AVAL(J)) / (d1*d2)
                if (i .ne. j) info(j,i) = info(i,j)
  530   CONTINUE
  540 CONTINUE
      DO 550 I=1,NP
        I1=I+1
            d1 = g(i+1,i) - g(1,i)
        info(i,i)= - 4.d0*(H(I1)+A0-TWO*AVAL(I)) / d1**2
  550 CONTINUE

	return
	end

