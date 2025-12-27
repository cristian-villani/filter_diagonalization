      subroutine reai_(e,L,emin,emax,hm,dh)
      implicit none
      integer MAXL
      parameter(MAXL=100)
      integer L,i
      real*8 e(MAXL), emin, emax, hm, dh
      real*8 hevf

      hevf = 1.0d0

      print *, "energies given as input values:"
      print *, "==============================="

      do i = 1,L
         e(i) = e(i) / hevf
         if(e(i) .lt. emin .or. e(i) .gt. emax) then
             print *, "Warning, energy value is outside of the interval"
         endif
         write(*,'(i6,f12.6)') i,e(i)
         e(i) = (e(i) - hm)/dh
      end do

      return
      end
