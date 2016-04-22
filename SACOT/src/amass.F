      FUNCTION AMASS(imassin, I)
C     Returns on-shell or MSbar mass 
c    of parton I for imassin=0 or 1, respectively. 
      implicit none
      integer imassin, i, ii, imass
      double precision amass

      Double precision POLEMASS, MSBMASS
      Logical MassesSet
      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), iMass, MassesSet
      
      AMASS = 0D0
      II = IABS(I)

c htj test      
!      print*,"==============="
!      print*,imassin,imass,MassesSet
!      print*,POLEMASS(1),POLEMASS(2),POLEMASS(3),POLEMASS(4)
!     > ,POLEMASS(5),POLEMASS(6)
!
!      print*,MSBMASS(1),MSBMASS(2),MSBMASS(3),MSBMASS(4)
!     > ,MSBMASS(5),MSBMASS(6)
!      stop

c htj test, comment out just for test
c      if (imassin.ne.imass.and..not. MassesSet) stop
c     >  'Quark masses are not synchronized; run MassSync to proceed'

      IF (II.GE.1 .and. II.LE.6)  THEN
        if (imassin .eq.0) then 
          AMASS = POLEMASS(II)
        else if (imassin .eq. 1) then 
          AMASS = MSBMASS(II) 
        else
          stop 'AMASS: no such mass type imassin'
        endif                   !imassin 
      ENDIF
      
      RETURN
      END
C
C***********************************************************
C
