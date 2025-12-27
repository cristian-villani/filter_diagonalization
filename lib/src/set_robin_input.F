      subroutine set_robin_input(dim1,L1,F_min1,F_max1,spec_min1,
     .                           spec_max1,prec1)
      implicit none
#include "filenames.h"
#include "spec.h"
#include "logicals.h"
#include "dimensions.h"
#include "f_minmax.h"
#include "gauss.h"
#include "cutprec.h"

      integer strlen
      integer L1,dim1
      real*8 F_min1,F_max1,spec_min1,spec_max1,prec1
c     namelist /blkinp/ F_min,F_max,L,Sfile,Hfile,psifile,
c    $     thresh,initialize,ana_only,restart,prec,scal,prtvec,vecfile,
c    $     hpsifile,filter_only,spec_min,spec_max,cut

      dim=dim1
      L=L1
      F_min= F_min1
      F_max=F_max1
      spec_min=spec_min1
      spec_max=spec_max1

      Sfile='blkfil_smat'
      Hfile='blkfil_hmat'
      psifile='blkfil_psimat'
      hpsifile='blkfil_hpsimat'
      vecfile='blkfil_vecs'

      prtvec=.true.
      thresh = 1.0D-4
c
c cut is set to "-1". If the input does not modify this value, then cut
c is optimized in the subroutine analyze.
c
      cut= -1.0d0
c     cut= 1.0D-9

      if(prec1.le.0.0d0) then
        prec= 1.0D-4
      else
        prec=prec1
      end if
      scal= 0.3D0
      
      packdim = L*(L+1)/2

      write(*,*) 'Input Parameters:'
      write(*,*) '-------------------------------------------'
      write(*,*) 'psifile    :', psifile(1:strlen(psifile))
      write(*,*) 'dim        :', dim
      write(*,*) 'L          :', L
      write(*,*) 'spec_min   :', spec_min
      write(*,*) 'spec_max   :', spec_max     
      write(*,*) 'F_min      :', F_min
      write(*,*) 'F_max      :', F_max
      write(*,*) 'scal       :', scal 
      write(*,*) 'thresh     :', thresh
      write(*,*) 'filter_only:', filter_only
      write(*,*) 'restart    :', restart
      write(*,*) 'initialize :', initialize
      if(.not.filter_only) then
        write(*,*) 'ana_only   :', ana_only
        write(*,*) 'prec       :', prec
        if(cut.ne.-1.0d0) then
           write(*,*) 'cut        :', cut
        else
           write(*,*) 'cut        : (optimized)'
        endif
        write(*,*) 'Sfile      :', Sfile(1:strlen(Sfile))
        write(*,*) 'Hfile      :', Hfile(1:strlen(Hfile))
        write(*,*) 'prtvec     :', prtvec
        if(prtvec) write(*,*) 'vecfile    :', vecfile(1:strlen(vecfile))
      endif
      write(*,*) '-------------------------------------------'

      return
      end
