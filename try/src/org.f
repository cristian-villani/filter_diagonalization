      subroutine reai_(e,L,emin,emax,hm,dh)
      implicit none
      integer MAXL
      parameter(MAXL=100)
      integer L,i
      real*8 hevf,e(MAXL),emin,emax,hm,dh

      namelist /energies/ hevf,e

      hevf=1.0d0
      read(5,energies)
      if(hevf.eq.0.0d0) then
        print *,"ERROR: conversion factor hefv can not be zero!"
        stop
      end if

      print *, "energies given as input values:"
      print *, "==============================="
      print *, " "
      do i=1,L
        e(i)=e(i)/hevf
        if(e(i).lt.emin.or.e(i).gt.emax)
     .  print *,"Warning, the energy value is outside of the interval"
        write(*,'(i6,f12.6)') i,e(i)
        e(i)=(e(i)-hm)/dh
      end do

      return
      end

c     program test
c     implicit none
c     real*8 en(500)

c     call reai(en,3)

c     stop
c     end
