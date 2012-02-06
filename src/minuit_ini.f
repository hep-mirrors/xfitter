      subroutine minuit_ini
C
C 2 Oct 2011: add extra parameters before starting minuit
C
      implicit none
C-------------------------------------
      character*72 minfile

      open ( 25, file='output/minuit.out.txt' )
      minfile='minuit.in.txt' 
      write(6,*) ' read minuit input params from file ',minfile
      call HF_errlog(12020504,
     +     'I: read minuit input params from file '//minfile) 
      open ( 24, file=minfile )
      open (  7, file='minuit.save.txt' )

      call mintio(24,25,7)

      return
      end

      subroutine ExtraParam
C--------------------------------------
C
C
C
C---------------------------------------
      implicit none
      include 'extrapars.inc'
      integer i, ierrf

C Add extra parameter:

      do i = 1,nExtraParam
         call mnparm(100+i,ExtraParamNames(i)
     $        ,ExtraParamValue(i)
     $        ,ExtraParamStep(i)
     $        ,ExtraParamMin(i)
     $        ,ExtraParamMax(i)
     $        ,ierrf)
         if (ierrf.ne.0) then
            print *,'Error adding extra parameter',i
            print *,'name, value, step, min, max are:',
     $           ExtraParamNames(i)
     $        ,ExtraParamValue(i)
     $        ,ExtraParamStep(i)
     $        ,ExtraParamMin(i)
     $        ,ExtraParamMax(i)
            print *,'Error code=',ierrf
            call HF_errlog(12020505,'F: Error in ExtraParam')
         else
            iExtraParamMinuit(i) = 100+i
         endif
      enddo
      end

      
