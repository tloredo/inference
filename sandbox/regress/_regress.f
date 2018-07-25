c -----------------------------------------------------------------------
      subroutine bess ( npts,y1,y1err,y2,y2err,cerr,
     &     a,b,avar,bvar,bvar_ifab,xi,zeta )
c -----------------------------------------------------------------------
c     Do the entire regression calculation for 4 slopes:
c     OLS(Y|X), OLS(X|Y), bisector, orthogonal
c

cf2py intent(in) npts, y1, y1err, y2, y2err, cerr
cf2py intent(out) a, b, avar, bvar, bvar_ifab
cf2py intent(cache,hide) xi, zeta

      implicit double precision (c,x,y,z)
      integer i,npts,nmod
      parameter (nmod = 4)
      dimension y1(npts),y1err(npts),y2(npts),y2err(npts)
      dimension cerr(npts),xi(npts,nmod),zeta(npts,nmod)
      double precision a(4),b(4),avar(4),bvar(4),bvar_ifab(4)
      double precision y1av,y2av,y1var,y2var,y1stddev,y2stddev
      double precision sig11var,sig22var,sig12var,sign
      double precision zetaav(4),zetadisp(4),xidisp(4),xiav(4)
      double precision covar_y1y2,covb1b2
c     
c     calculate sigma's for datapoints using length of conf. intervals
c
      sig11var = 0.
      sig22var = 0.
      sig12var = 0.
      do 2 i=1,npts
         sig11var = sig11var + y1err(i)**2
         sig22var = sig22var + y2err(i)**2
c
c     old version; mistake in treating cerr like a standard deviation
c     instead of a covariance
c        sig12var = sig12var + cerr(i)**2 
c                
c     new verions: Apr 28 1998
c
         sig12var = sig12var + cerr(i)
 2    continue
      sig11var = sig11var/real(npts)
      sig22var = sig22var/real(npts)
      sig12var = sig12var/real(npts)
c     
c     calculate means and variances
c     
      call class(y1,npts,y1av,y1stddev)
      call class(y2,npts,y2av,y2stddev)
      y1var = y1stddev**2
      y2var = y2stddev**2
      covar_y1y2 = 0.
      do 5 i=1,npts
         covar_y1y2 = (y1(i)-y1av)*(y2(i)-y2av) + covar_y1y2
 5    continue
      covar_y1y2 = covar_y1y2/real(npts)
c      
c     compute the regression slopes for OLS(Y2|Y1), OLS(Y1|Y2),
c     bisector, and orthogonal.
c      
      b(1) = (covar_y1y2 - sig12var)/(y1var - sig11var)
      b(2) = (y2var - sig22var)/(covar_y1y2 - sig12var)
      b(3) = ( b(1)*b(2) - 1.0 
     &     + sqrt((1.0 + b(1)**2)*(1.0 + b(2)**2)) ) /
     &     (b(1)+b(2))
      if ( covar_y1y2.lt.0. ) then
         sign = -1.
      else
         sign = 1.
      endif
      b(4) = 0.5*((b(2)-(1./b(1))) 
     &     + sign * sqrt(4.+(b(2)-(1./b(1)))**2))
c      
c     compute intercepts for above 4 cases:
c      
      do 10 i=1,4
         a(i) = y2av - b(i)*y1av
 10   continue
c      
c     set up variables to calculate standard deviations of slope
c     and intercept (MAB renamed: chi -> xi, xi -> zeta to be consistent 
c     with text.)
c
      do 15 i=1,npts
         xi(i,1) = ( (y1(i)-y1av) * (y2(i)-b(1)*y1(i)-a(1)) + 
     &        b(1)*y1err(i)**2 ) / (y1var-sig11var)
         zeta(i,1) = y2(i) - b(1)*y1(i) - y1av*xi(i,1)
         xi(i,2) = ( (y2(i)-y2av) * (y2(i)-b(2)*y1(i)-a(2)) - 
     &        y2err(i)**2 ) / covar_y1y2
         zeta(i,2) = y2(i) - b(2)*y1(i) - y1av*xi(i,2)
         xi(i,3) = xi(i,1) * 
     &        (1.+b(2)**2)*b(3) / 
     &        ((b(1)+b(2))*sqrt((1.+b(1)**2)*(1.+b(2)**2))) + 
     &        xi(i,2) * 
     &        (1.+b(1)**2)*b(3) /
     &        ((b(1)+b(2))*sqrt((1.+b(1)**2)*(1.+b(2)**2)))
         zeta(i,3) = y2(i) - b(3)*y1(i) - y1av*xi(i,3)
         xi(i,4) = xi(i,1) * 
     &        b(4)/(b(1)**2*sqrt(4.+(b(2)-1./b(1))**2)) +
     &        xi(i,2)*b(4)/sqrt(4.+(b(2)-1./b(1))**2)
         zeta(i,4) = y2(i) - b(4)*y1(i) - y1av*xi(i,4)
 15   continue
c      
c     calculate variance for all a and b
c
      do 20 i=1,4
         call nclass(npts,i,xi,xiav,xidisp)
         call nclass(npts,i,zeta,zetaav,zetadisp)
         bvar(i) = xidisp(i)**2/real(npts)
         avar(i) = zetadisp(i)**2/real(npts)
 20   continue
c      
c     alternate slope variances for b3 and b4 via IFAB formulae;
c     requires calculating first covariance for b1,b2
c
      covb1b2 = xidisp(1)*xidisp(2)/real(npts)
      bvar_ifab(3) = ( b(3)**2 / 
     &     (((b(1)+b(2))**2)*(1.+b(1)**2)*(1.+b(2)**2)) ) *
     &     ( (1.+b(2)**2)**2*bvar(1) + (1.+b(1)**2)**2*bvar(2) + 
     &     2.*(1.+b(1)**2)*(1.+b(2)**2)*covb1b2 )
      bvar_ifab(4) = ( b(4)**2 / 
     &     (4.*b(1)**2 + (b(1)*b(2)-1.)**2 ) ) *
     &     ( bvar(1)/b(1)**2 + 2.*covb1b2 + b(1)**2*bvar(2) )
c
      return
      end
c -----------------------------------------------------------------------
      subroutine class (x,n,xav,xstddev)
c -----------------------------------------------------------------------
c     Calculate mean and standard deviation of the array X of N numbers.
c     
      implicit double precision (x)
      dimension x(n)
      integer i
      double precision xav,xstddev
c
      xav = 0.0
      xstddev = 0.0
      
      do 10 i=1,n
         xav = xav + x(i)
 10   continue
      xav = xav/dfloat(n)
c      
      do 20 i=1,n
         xstddev = (x(i) - xav)**2 + xstddev
 20   continue
      xstddev = sqrt(xstddev/dfloat(n))
c      
      return
      end
c -----------------------------------------------------------------------
      subroutine nclass (npts,i,x,xav,xstddev)
c -----------------------------------------------------------------------
c     Calculate mean and standard deviation of the array X of N numbers.
c     
      implicit double precision (x)
      integer npts, nmod
      parameter (nmod = 4)
      dimension x(npts,nmod),xav(nmod),xstddev(nmod)
      integer i,j
c
      xav(i) = 0.0
      xstddev(i) = 0.0
c
      do 10 j=1,npts
         xav(i) = xav(i) + x(j,i)
 10   continue
      xav(i) = xav(i)/dfloat(npts)
c      
      do 20 j=1,npts
         xstddev(i) = (x(j,i) - xav(i))**2 + xstddev(i)
 20   continue
      xstddev(i) = sqrt(xstddev(i)/dfloat(npts))
c
      return
      end
c -----------------------------------------------------------------------------
      subroutine wlss ( npts,x1,y2,y2err, a,avar,b,bvar,evar,res,ierr )
c -----------------------------------------------------------------------------
c     Core subroutine: Calculate weighted least-squares for model
c     with intrinsic scatter.
c

cf2py intent(in) npts, x1, y2, y2err
cf2py intent(out) a, avar, b, bvar, evar, ierr
cf2py intent(cache,hide) res

      implicit real*8 (x,y,r)
      integer npts
      dimension x1(npts),y2(npts),y2err(npts),res(npts)
      double precision a,avar,b,bvar
      double precision avec(6),bvec(6),sigavec(6),sigbvec(6)
      double precision resav,resvar,sig22var,evar,estarvar
      double precision wsum,wsumx1,wsumx1y2,wsumy2,wsumx1x1,denom
c      
c     Step 1: call sixlin to get linear regression coefficients from OLS.
c     (NB: Only need j=1 case: OLS(Y|X))
c
      call sixlin(npts,x1,y2,avec,sigavec,bvec,sigbvec,ierr)
      if (ierr .eq. 0) return
c
c     Step 2: Calculate residuals and average residual
c     
      resav = 0.
      do 20 i=1,npts
         res(i) = y2(i) - avec(1) - bvec(1)*x1(i)
         resav = resav + res(i)
 20   continue
      resav = resav/dfloat(npts)
c
c     Step 3: Obtain estimate of variance of the intrinsic scatter
c
      resvar = 0.
      sig22var = 0.
      do 30 i=1,npts
         resvar = resvar + (res(i)-resav)**2
         sig22var = sig22var + y2err(i)**2
 30   continue
      sig22var = sig22var/dfloat(npts)
      resvar = resvar/dfloat(npts)
      evar = resvar - sig22var
c     
c     Prepare for b and a calculation
c     
      wsum = 0.
      wsumx1 = 0.
      wsumx1y2 = 0.
      wsumy2 = 0.
      wsumx1x1 = 0.
      do 40 i=1,npts
         estarvar = evar + y2err(i)**2
         wsum = wsum + 1./estarvar
         wsumx1 = wsumx1 + x1(i)/estarvar
         wsumx1y2 = wsumx1y2 + x1(i)*y2(i)/estarvar
         wsumy2 = wsumy2 + y2(i)/estarvar
         wsumx1x1 = wsumx1x1 + x1(i)*x1(i)/estarvar
 40   continue
c
c     Calcalculate b and a estimates and estimated variances
c
      denom =  wsum * wsumx1x1 - wsumx1**2
      b = ( wsum * wsumx1y2 - wsumx1 * wsumy2 ) / denom
      a = ( wsumx1x1 * wsumy2 - wsumx1 * wsumx1y2 ) / denom
      bvar = wsum / denom
      avar = wsumx1x1 / denom
c
      return
      end
C*******************************************************************
C************************* Subroutine datstt ***********************
C*******************************************************************
C
       Subroutine datstt(N,X,Y,Xavg,Varx,Yavg,Vary,Varxy,Rho,ierr)
C
C     *  Subroutine To Compute Simple Statistical Properties Of A  *
C     *  Bivariate Data Set.  It Gives The Variance In X, Variance *
C     *  In Y, And Pearson'S Linear Correlation Coefficient.       *
C
C     If Sxy=0 (invalid), stop calc'n and return ierr=0; else ierr=1.

cf2py intent(in) N, X, Y
cf2py intent(out) Xavg, Varx, Yavg, Vary, Varxy, Rho, ierr

       Implicit Real*8(A-H,O-Z)
       Dimension X(N),Y(N)
C
C
C     *         Initializations
C
       S1 = 0.0
       S2 = 0.0
       Sxx = 0.0
       Syy = 0.0
       Sxy = 0.0
       Rn = Dfloat(N)
       ierr = 1
C
C     *         Compute Averages And Sums
C
       Do 100 I=1,N
         S1 = S1 + X(I)
         S2 = S2 + Y(I)
  100  Continue
       Xavg = S1/Rn
       Yavg = S2/Rn
       Do 200 I = 1, N
       Sxx  = Sxx  + (X(I) - Xavg)**2
       Syy  = Syy  + (Y(I) - Yavg)**2
       Sxy  = Sxy  + (X(I) - Xavg)*(Y(I) - Yavg)
 200   Continue
       If(Sxy.eq.0.0) Then
         ierr = 0
         return
       Endif
C
C     *         Compute And Return Results
C
       Varx = Sxx/Rn
       Vary = Syy/Rn
       Varxy = Sxy/Rn
       Rho = Sxy/(Dsqrt(Sxx*Syy))
C
       Return
       End
C************************************************************************
C************************* Subroutine sixlin ****************************
C************************************************************************
C
      Subroutine sixlin(N,X,Y,A,Siga,B,Sigb,ierr)
C
C     *                     Six Linear Regressions
C     *     Written By T. Isobe, G. J. Babu And E. D. Feigelson
C     *               Center For Space Research, M.i.t.
C     *                             And
C     *              The Pennsylvania State University
C     *
C     *                   Rev. 1.0,   September 1990
C     *
C     *       This Subroutine Provides Linear Regression Coefficients
C     *    Computed By Six Different Methods Described In Isobe,
C     *    Feigelson, Akritas, And Babu 1990, Astrophysical Journal
C     *    And Babu And Feigelson 1990, Subm. To Technometrics.
C     *    The Methods Are Ols(Y/X), Ols(X/Y), Ols Bisector, Orthogonal,
C     *    Reduced Major Axis, And Mean-Ols Regressions.
C     *
C     *    Input
C     *         X(I) : Independent Variable
C     *         Y(I) : Dependent Variable
C     *            N : Number Of Data Points
C     *
C     *    Output
C     *         A(J) : Intercept Coefficients
C     *         B(J) : Slope Coefficients
C     *      Siga(J) : Standard Deviations Of Intercepts
C     *      Sigb(J) : Standard Deviations Of Slopes
C     *     Where J = 1, 6.
C     *
C     *    Error Returns
C     *         Calculation Is Stopped When Division By Zero Will
C     *         Occur (I.e. When Sxx, Sxy, Or Either Of The Ols
C     *         Slopes Is Zero).  Returns ierr=0; otherwise ierr=1.
C
C

cf2py intent(in) N, X, Y
cf2py intent(out) A, Siga, B, Sigb, ierr

      Implicit Real*8(A-H,O-Z)
      Dimension X(N),Y(N)
      Dimension A(6),B(6),Siga(6),Sigb(6)
C     
C     *         Initializations
C
      S1 = 0.0
      S2 = 0.0
      Sxx = 0.0
      Syy = 0.0
      Sxy = 0.0
      Sum1 = 0.0
      Sum2 = 0.0
      Sum3 = 0.0
      Rn = N
      ierr = 1
      Do 50 J = 1, 6
         A(J)    = 0.0
         Siga(J) = 0.0
         B(J)    = 0.0
         Sigb(J) = 0.0
 50   Continue
C     
C     *         Compute Averages And Sums
C
      Do 100 I=1,N
         S1 = S1 + X(I)
         S2 = S2 + Y(I)
 100  Continue
      Xavg = S1/Rn
      Yavg = S2/Rn
      Do 200 I = 1, N
         X(I) = X(I) - Xavg
         Y(I) = Y(I) - Yavg
         Sxx  = Sxx  + X(I)**2
         Syy  = Syy  + Y(I)**2
         Sxy  = Sxy  + X(I)*Y(I)
 200  Continue
      If(Sxy.eq.0.0) Then
         ierr = 0
         return
      Endif
      Sign = 1.0
      If(Sxy .lt. 0.0) Sign = -1.0
C
C     *               Compute The Slope Coefficients
C
      B(1) = Sxy/Sxx
      B(2) = Syy/Sxy
      B(3) = (B(1)*B(2) - 1.0
     +     + Dsqrt((1.0 + B(1)**2)*(1.0 + B(2)**2)))/(B(1) + B(2))
      B(4) = 0.5*(B(2) - 1.0/B(1)
     +     + Sign*Dsqrt(4.0 + (B(2) - 1.0/B(1))**2))
      B(5) = Sign*Dsqrt(B(1)*B(2))
      B(6) = 0.5*(B(1) + B(2))
C     
C     *            Compute Intercept Coefficients
C     
      Do 300 J = 1, 6
         A(J) = Yavg - B(J)*Xavg
 300  Continue
C     
C     *     Prepare For Computation Of Variances
C     
      Gam1 = B(3)/((B(1) + B(2))
     +     *Dsqrt((1.0 + B(1)**2)*(1.0 + B(2)**2)))
      Gam2 = B(4)/(Dsqrt(4.0*B(1)**2 + (B(1)*B(2) - 1.0)**2))
      Do 400 I = 1, N
         Sum1 = Sum1 + (X(I)*(Y(I) - B(1)*X(I)))**2
         Sum2 = Sum2 + (Y(I)*(Y(I) - B(2)*X(I)))**2
         Sum3 = Sum3 + X(I)*Y(I)*(Y(I) - B(1)*X(I))*(Y(I) - B(2)*X(I))
 400  Continue
      Cov = Sum3/(B(1)*Sxx**2)
C     
C     *    Compute Variances Of The Slope Coefficients
C     
      Sigb(1) = Sum1/(Sxx**2)
      Sigb(2) = Sum2/(Sxy**2)
      Sigb(3) = (Gam1**2)*(((1.0 + B(2)**2)**2)*Sigb(1)
     +     + 2.0*(1.0 + B(1)**2)*(1.0 + B(2)**2)*Cov
     +     + ((1.0 +B(1)**2)**2)*Sigb(2))
      Sigb(4) = (Gam2**2)*(Sigb(1)/B(1)**2 + 2.0*Cov
     +     + B(1)**2*Sigb(2))
      Sigb(5) = 0.25*(B(2)*Sigb(1)/B(1)
     +     + 2.0*Cov + B(1)*Sigb(2)/B(2))
      Sigb(6) = 0.25*(Sigb(1) + 2.0*Cov + Sigb(2))
C     
C     *   Compute Variances Of The Intercept Coefficients
C     
      Do 500 I = 1, N
         Siga(1) = Siga(1) + ((Y(I) - B(1)*X(I))
     +        *(1.0 - Rn*Xavg*X(I)/Sxx))**2
         Siga(2) = Siga(2) + ((Y(I) - B(2)*X(I))
     +        *(1.0 - Rn*Xavg*Y(I)/Sxy))**2
         Siga(3) = Siga(3) + ((X(I)*(Y(I)
     +        - B(1)*X(I))*(1.0 + B(2)**2)/Sxx
     +        + Y(I)*(Y(I) - B(2)*X(I))*(1.0 + B(1)**2)/Sxy)
     +        *Gam1*Xavg*Rn - Y(I) + B(3)*X(I))**2
         Siga(4) = Siga(4) + ((X(I)*(Y(I) - B(1)*X(I))/Sxx
     +        + Y(I)*(Y(I) - B(2)*X(I))*(B(1)**2)/Sxy)*Gam2
     +        *Xavg*Rn/Dsqrt(B(1)**2) - Y(I) + B(4)*X(I))**2
         Siga(5) = Siga(5) + ((X(I)*(Y(I)
     +        - B(1)*X(I))*Dsqrt(B(2)/B(1))/Sxx
     +        + Y(I)*(Y(I) - B(2)*X(I))*Dsqrt(B(1)/B(2))/Sxy)
     +        *0.5*Rn*Xavg - Y(I) + B(5)*X(I))**2
         Siga(6) = Siga(6) + ((X(I)*(Y(I) - B(1)*X(I))/Sxx
     +        + Y(I)*(Y(I) - B(2)*X(I))/Sxy)
     +        *0.5*Rn*Xavg - Y(I) + B(6)*X(I))**2
 500  Continue
C     
C     *  Convert Variances To Standard Deviations
C     
      Do 600 J = 1, 6
         Sigb(J) = Dsqrt(Sigb(J))
         Siga(J) = Dsqrt(Siga(J))/Rn
 600  Continue
C     
C     *  Return Data Arrays To Their Original Form
C
      Do 900 I = 1, N
         X(I) = X(I) + Xavg
         Y(I) = Y(I) + Yavg
 900  Continue
C     
      Return
      End
