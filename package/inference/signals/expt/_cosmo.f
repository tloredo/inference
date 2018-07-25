c-----------------------------------------------------------------------
c
c   Functions for basic FRW cosmology calculations in Python.
c
c        May 2006  TJL  (Adapted from code for DLW 2000, written in 1998.)
c        Mods:  
c
c-----------------------------------------------------------------------
        subroutine set_frw_cosmo(hh, OM, OL)

cf2py intent(in) hh, OM, OL

c+++  Arguments:
            real*8 hh, OM, OL

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 c_fac
        parameter (c_fac =  2.99792458d8)
c        parameter (c_fac =  2997.92458) in Mpc

        h = hh
        O_M = OM
        O_L = OL
        O_k = 1. - OM - OL
        rtO_k = sqrt(abs(O_k))
        c_h = c_fac / h

        return
        end
    
c-----------------------------------------------------------------------
        function ldist (z)

c---        "Luminosity DISTance"
c
c        The luminosity distance, in units of 10pc.
c
c        The dependence on the cosmology is via the common block.

cf2py intent(in) z
cf2py intent(out) ldist

c+++  Arguments:
        real*8 z, ldist

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 tol, zero, c_fac
        parameter (tol = 1.d-6, zero = 0.)
        parameter (c_fac =  2.99792458d8)
        real*8 dd_l

c+++  Functions:
        real*8 dd_l_ig
        external dd_l_ig

        call qromb(dd_l_ig, zero, z, dd_l, tol)
c        write(9,'(a,3(1pg12.4))') '***********', z, dd_l, 
c     *       rtO_k*dd_l/(2.*3.14159265358d0)
        if (O_k .lt. zero) then
            dd_l = (1.+z) * abs(sin(rtO_k*dd_l)) / rtO_k
        else if (O_k .gt. zero) then
            dd_l = (1.+z) * sinh(rtO_k*dd_l) / rtO_k
        else
            dd_l = (1.+z) * dd_l
        endif
        ldist = c_h * dd_l

        return
        end

c-----------------------------------------------------------------------
        function coord_int (z)

c---        "COORDinate INTegral"
c
c        Integral for the dimensionless coordinate distance;
c        multiply by l_Hub/R_0 to get the coordinate distance.
c
c        The dependence on the cosmology is via the common block.

cf2py intent(in) z
cf2py intent(out) coord_int

c+++  Arguments:
        real*8 z, coord_int

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 tol, zero, c_fac
        parameter (tol = 1.d-6, zero = 0.)
        parameter (c_fac =  2.99792458d8)

c+++  Functions:
        real*8 dd_l_ig
        external dd_l_ig

        call qromb(dd_l_ig, zero, z, coord_int, tol)

        return
        end

c-----------------------------------------------------------------------
        subroutine nldist (n, zvals, ldvals)

c---        "N Luminosity DISTances"
c
c        The luminosity distance, in units of 10pc.
c
c        The dependence on the cosmology is via the common block.

cf2py intent(in) n, zvals
cf2py intent(out) ldvals

c+++  Arguments:
        integer n
        real*8 zvals(n), ldvals(n)

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 tol, zero, c_fac
        parameter (tol = 1.d-6, zero = 0.)
        parameter (c_fac =  2.99792458d8)
        real*8 dd_l
        integer i

c+++  Functions:
        real*8 dd_l_ig
        external dd_l_ig

        do i=1, n
            call qromb(dd_l_ig, zero, zvals(i), dd_l, tol)
            if (O_k .lt. zero) then
                dd_l = (1.+zvals(i)) * abs(sin(rtO_k*dd_l)) / rtO_k
            else if (O_k .gt. zero) then
                dd_l = (1.+zvals(i)) * sinh(rtO_k*dd_l) / rtO_k
            else
                dd_l = (1.+zvals(i)) * dd_l
            endif
            ldvals(i) = c_h * dd_l
        end do

        return
        end

c-----------------------------------------------------------------------
        function mu_z (z)

c---        "distance modulus (mu) from redshift (z)"
c
c        The luminosity distance, in units of 10pc (as required
c        for the distance modulus), is calculated first.
c        Then mu = 5 log (d/10pc).
c
c        The dependence on the cosmology is via the common block.

cf2py intent(in) z
cf2py intent(out) mu_z

c+++  Arguments:
        real*8 mu_z, z

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 tol, zero, c_fac
        parameter (tol = 1.d-6, zero = 0.)
        parameter (c_fac =  2.99792458d8)
        real*8 dd_l

c+++  Functions:
        real*8 dd_l_ig
        external dd_l_ig

        call qromb(dd_l_ig, zero, z, dd_l, tol)
c        write(9,'(a,3(1pg12.4))') '***********', z, dd_l, 
c     *       rtO_k*dd_l/(2.*3.14159265358d0)
        if (O_k .lt. zero) then
            dd_l = (1.+z) * abs(sin(rtO_k*dd_l)) / rtO_k
        else if (O_k .gt. zero) then
            dd_l = (1.+z) * sinh(rtO_k*dd_l) / rtO_k
        else
            dd_l = (1.+z) * dd_l
        endif
        dd_l = c_h * dd_l
        if (dd_l .le. 0.) print *, '0 dd_l in mu_z: ', z, dd_l
        mu_z = 5.*log10(dd_l)

        return
        end

c-----------------------------------------------------------------------
        function mu_z_noh (z)

c---        "distance modulus (mu) from redshift (z) - NO H_0 factor"
c
c        The 5 log (d/10pc), but without the c/H_0 factor in d.
c
c        The dependence on the cosmology is via the common block.

cf2py intent(in) z
cf2py intent(out) mu_z_noh

c+++  Arguments:
        real*8 mu_z_noh, z

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 tol, zero, c_fac
        parameter (tol = 1.d-6, zero = 0.)
        parameter (c_fac =  2.99792458d8)
        real*8 dd_l

c+++  Functions:
        real*8 dd_l_ig
        external dd_l_ig

        call qromb(dd_l_ig, zero, z, dd_l, tol)
c        write(9,'(a,3(1pg12.4))') '***********', z, dd_l, 
c     *       rtO_k*dd_l/(2.*3.14159265358d0)
        if (O_k .lt. zero) then
            dd_l = (1.+z) * abs(sin(rtO_k*dd_l)) / rtO_k
        else if (O_k .gt. zero) then
            dd_l = (1.+z) * sinh(rtO_k*dd_l) / rtO_k
        else
            dd_l = (1.+z) * dd_l
        endif
        if (dd_l .le. 0.) print *, '0 dd_l in mu_z_noh: ', z, dd_l
        mu_z_noh = 5.*log10(dd_l)

        return
        end

c-----------------------------------------------------------------------
        function dd_l_ig (z)

c---        "Dimensionless Distance_Luminosity InteGrand"
c
c        The integrand for the luminosity distance integral,
c        not including the c/H_0 factor.

c+++  Arguments:
        real*8 dd_l_ig, z

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 z1

c        dd_l_ig = 1./ sqrt( (1. + z)**2 * (1. + O_M*z) - z*(2.+z)*O_L )
        z1 = 1. + z
        dd_l_ig = 1./ sqrt( (O_M*z1 + O_k)*z1**2 + O_L )
c        write(9,'(4(1pg12.4))') z, dd_l_ig

        return
        end

c-----------------------------------------------------------------------
        function E_Hub (z)

c---        "E HUBble"
c
c        The dimensionless Hubble parameter vs. z, E(z).

cf2py intent(in) z
cf2py intent(out) E_Hub

c+++  Arguments:
        real*8 E_Hub, z

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 z1

        z1 = 1. + z
        E_Hub = sqrt( (O_M*z1 + O_k)*z1**2 + O_L )

        return
        end

c-----------------------------------------------------------------------
        function vol_elem (z)

c---        "VOLume ELEMent"
c
c        The volume per unit redshift per unit steradian.
c
c        The dependence on the cosmology is via the common block.

cf2py intent(in) z
cf2py intent(out) vol_elem

c+++  Arguments:
        real*8 z, vol_elem

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 tol, zero, c_fac
        parameter (tol = 1.d-6, zero = 0.)
        parameter (c_fac =  2.99792458d8)
        real*8 dist, z1

c+++  Functions:
        real*8 dd_l_ig
        external dd_l_ig

c---  This call returns the l.o.s. comoving distance in dV.
        call qromb(dd_l_ig, zero, z, dist, tol)
c        write(9,'(a,3(1pg12.4))') '***********', z, dist, 
c     *       rtO_k*dV/(2.*3.14159265358d0)

c---  Find the transverse comoving distance.
        if (O_k .lt. zero) then
            dist = abs(sin(rtO_k*dist)) / rtO_k
        else if (O_k .gt. zero) then
            dist = sinh(rtO_k*dist) / rtO_k
        endif

        z1 = 1. + z
        vol_elem = c_h * dist**2 / sqrt( (O_M*z1 + O_k)*z1**2 + O_L )

        return
        end

c-----------------------------------------------------------------------
        function dlbt_z (z)

c---        "Dimensionless Look Back Time from redshift (z)"
c
c        Returns H_0*t, where t is the look-back time to redshift z.
c
c        The dependence on the cosmology is via the common block.

cf2py intent(in) z
cf2py intent(out) dlbt_z

c+++  Arguments:
        real*8 dlbt_z, z

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 tol, zero
        parameter (tol = 1.d-6, zero = 0.)

c+++  Functions:
        real*8 lbt_ig
        external lbt_ig

        call qromb(lbt_ig, zero, z, dlbt_z, tol)

        return
        end

c-----------------------------------------------------------------------
        function lbt_ig (z)

c---        "Look Back Time InteGrand"
c
c        The integrand for the dimensionless look-back time integral.

c+++  Arguments:
        real*8 lbt_ig, z

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 z1

        z1 = 1. + z
        lbt_ig = 1. / (z1 * sqrt( O_L + z1**2 * (z1*O_M + O_k) ) )

        return
        end

c-----------------------------------------------------------------------
        function dage ()

c---        "Dimensionless AGE"
c
c        Returns H_0*t0.
c
c        The dependence on the cosmology is via the common block.

cf2py intent(out) dage

c+++  Arguments:
        real*8 dage

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

c+++  Locals:
        real*8 tol, zero, one
        parameter (tol = 1.d-6, zero = 0., one = 1.)

c+++  Functions:
        real*8 age_ig
        external age_ig

        call qromb(age_ig, zero, one, dage, tol)
 
        return
        end

c-----------------------------------------------------------------------
        function age_ig (a)

c---        "AGE InteGrand"
c
c        The integrand for the age vs. scale factor.

c+++  Arguments:
        real*8 age_ig, a

c+++  Globals:
        real*8 h, O_M, O_L, O_k, rtO_k, c_h
        common /cosmo_cb/ h, O_M, O_L, O_k, rtO_k, c_h

        age_ig = sqrt( a / (O_M + a*(O_k + a**2*O_L)) )

        return
        end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
        subroutine set_ffid(Ffid)

c       Set the fiducial flux, setting the scale for luminosities.

cf2py intent(in) Ffid

c+++  Arguments:
            real*8 Ffid

c+++  Globals:
        real*8 F_fid
        common /ffid_cb/ F_fid

        F_fid = Ffid

        return
        end
    
c-----------------------------------------------------------------------
        subroutine set_lum(lum)

c       Set standard candle dimensionless luminosity.

cf2py intent(in) lum

c+++  Arguments:
            real*8 lum

c+++  Globals:
        real*8 lambda
        common /stdcdl_cb/ lambda

        lambda = lum

        return
        end
    
c-----------------------------------------------------------------------
        function flux2z(F)

c       Find redshift producing flux F from a standard candle source.
c       *** Note this assumes a flat cosmology! ***
c       (sin or sinh function needed for non-flat)

cf2py intent(in) F
cf2py intent(out) flux2z

c+++  Arguments:
            real*8 F, flux2z

c+++  Globals:
        real*8 F_fid
        common /ffid_cb/ F_fid
        real*8 lambda
        common /stdcdl_cb/ lambda

c+++  Parameters:
        real*8 c100, pi4, tol, ddel, zero
        integer itmax
        parameter (c100 = 2997.9246, pi4 = 4.*3.141592653589793d0)
        parameter (tol = 0.0001, ddel = 0.001, zero = 0., itmax = 25)

c+++  Locals:
        real*8 E_Hub, coord_int
        real*8 x, z, zp1, del, Icur, xcur
        integer iter

c---  Solve g(z) = x, for g(z) = (1+z)*I(z), by Newton-Raphson.

c---  This is the target x value:
        x = sqrt(lambda * F_fid / F)

c---  Start with a guess from low-z behavior, I(z) ~ z.
        z = -0.5 + 0.5*sqrt(1. + 4.*x)

c---  NR iterations:
        iter = 0
        del = min(0.003*z,ddel)
10      Icur = coord_int(z)
        zp1 = 1. + z
        xcur = zp1*Icur
        if (abs(x-xcur) .gt. tol*x) then
            iter = iter + 1
            if (iter .gt. itmax) pause 'flux2z:  Too many NRs for z!'
            dg = Icur + zp1/E_Hub(z)
            z = z + (x-xcur)/dg
            if (z .lt. 0.) z = 0.5*(zp1-1.)
            goto 10
        endif

        flux2z = z

        return
        end
    
c***********************************************************************
c
c       Following are 1-d Romberg integration routines, from 
c       Numerical Recipes.
c
c***********************************************************************
      SUBROUTINE QROMB(FUNC,A,B,SS,eps)

        external FUNC
        real*8 func, a, b, ss, dss, eps, zero
        integer j, l, jmax, jmaxp, k, km
      PARAMETER(JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1, zero=0.)
      real*8 h(jmaxp), s(jmaxp)

        if(a.eq.b) then
           ss=0.
           return
        endif

      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          L=J-KM
          CALL POLINT(H(L),S(L),K,zero,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
          if (SS .eq. 0.) return
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
11    CONTINUE
      PAUSE 'QROMB:  Too many steps.'
      END
 

c-----------------------------------------------------------------------
      SUBROUTINE TRAPZD(FUNC,A,B,S,N)

        save
        real*8 func, a, b, s, del, x, sum, tnm
        integer n, it, j
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
 
c-----------------------------------------------------------------------
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)

        real*8 X, Y, DY, DIF, DIFT, HO, HP, W, DEN
        integer N, NMAX, NS, I, M
      PARAMETER (NMAX=10) 
      real*8 XA(N),YA(N),C(NMAX),D(NMAX)

      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
c          IF(DEN.EQ.0.)PAUSE 'Polint error!'
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
