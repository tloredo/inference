C     -*- fortran -*-
C     This file is autogenerated with f2py (version:2)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrappsigc (psigcf2pywrap, s, n_on, t_on, c)
      external psigc
      real*8 s
      integer n_on
      real*8 t_on
      real*8 c(n_on + 1)
      real*8 psigcf2pywrap, psigc
      psigcf2pywrap = psigc(s, n_on, t_on, c)
      end


      subroutine f2pywrapslml_offset (slml_offsetf2pywrap, n_on, t
     &_on, n_off, t_off)
      external slml_offset
      integer n_on
      real*8 t_on
      integer n_off
      real*8 t_off
      real*8 slml_offsetf2pywrap, slml_offset
      slml_offsetf2pywrap = slml_offset(n_on, t_on, n_off, t_off)
      end


      subroutine f2pywrapslmlike (slmlikef2pywrap, s, n_on, t_on, 
     &n_off, t_off, offset, cut, nt, ierr)
      external slmlike
      real*8 s
      integer n_on
      real*8 t_on
      integer n_off
      real*8 t_off
      real*8 offset
      real*8 cut
      integer nt
      integer ierr
      real*8 slmlikef2pywrap, slmlike
      slmlikef2pywrap = slmlike(s, n_on, t_on, n_off, t_off, offse
     &t, cut, nt, ierr)
      end

