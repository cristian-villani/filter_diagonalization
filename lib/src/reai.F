      subroutine reai(energy,L,emin,emax,hm,dh)
      implicit none
#include "const.h"
c     integer MAXL
c     parameter(MAXL=100)
      integer L,i
      real*8 hevf,e(MAXL),emin,emax,hm,dh,energy(L)

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
        energy(i)=e(i)/hevf
        if(energy(i).lt.emin.or.energy(i).gt.emax)
     .  print *,"Warning, the energy value is outside of the interval"
        write(*,'(i6,f12.6)') i,energy(i)
        energy(i)=(energy(i)-hm)/dh
      end do

      return
      end
