      SUBROUTINE MASSSYNC

C =============================================================================
C GroupName: MASSSYNC
C Description: Routines created to check whether values of pole masses
C and MSbar masses are consistent "Synchronized". 
C It performs the conversion from pole masses and MSbar masses and vice-versa.
C ==============================================================================
C M.G. May 06 2011

      IMPLICIT NONE

      DOUBLE PRECISION qms,Factor, Small, alpi1, Del, Delta,X, alpi
      DOUBLE PRECISION POLEMASS, MSBMASS, AL, FSWITCH
      INTEGER IMass,Nfl,NFTOTAL,NORDER, N, N2, NFACTIVE
      integer lfit, lerr, lprt, loutput,lptp
      LOGICAL MassesSet, SET
      

      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), IMass, MassesSet
     >   /QCDPAR_CT / FSWITCH, AL, NFTOTAL, NORDER, SET
     >  /FitSwh/ Lfit, LErr, Lprt, Loutput
 

      double precision alam, amhat, amn
      integer nhq
      logical AmHatSet
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet


      data small /1.0d-5/ 

      DELTA(X)=1.2337d0 *X - 0.597d0*X**2 + 0.230d0*X**3

      
cpn11 Initialize alpha_s using the boundary condition from the previous run
      Call AlfSet(0d0,0d0) 

      lptp=Mod(Lprt,10) 
      if(lptp.eq. 1) write(*,*)
     >  'MASSSYNC: flavor   MS-bar Mass  Pole Mass (in GeV)'
 
      DO N = 1, NFTOTAL
cpn11   Decide if the input mass is an on-shell or MSbar mass
        IF(IMass.EQ.0) THEN
          QMS =POLEMASS(N)
        ELSE IF (IMass.eq.1) then
          QMS =MSBMASS(N)
        else
          stop 'MassSync: iMass is not defined'
        ENDIF !IMass 1

cpn11 Quarks with  masses lighter than AMN=0.5 GeV (in DatQCD) 
c     are light. Skip them.
        if (QMS.le.Amn) cycle

cpn11   Alpha_s is defined slightly above the threshold
        alpi1 = alpi( QMS +small )

        Factor = 1d0
        IF(IMass.EQ.0) THEN 
cmarco11     Convert the pole mass to the MSbar mass according to
c Chetyrkin et.al., Comput.Phys.Commun.133:43-65,2000; hep-ph/0004189,
c Eq.(14)

          if (norder.ge.2) then 
            Factor = Factor - 4d0/3d0* alpi1
c$$$          else if (norder.gt.2) then 
c$$$            print *, 
c$$$     >        'Conversion of the pole mass to MSbar mass ',
c$$$     >        'is not implemented for norder > 2'
c$$$            stop
          endif                 !norder 
            MSBMASS(N) = QMS*Factor
        ELSE
cmarco11     Convert the MS bar mass to the pole mass according to
c Chetyrkin et.al., Comput.Phys.Commun.133:43-65,2000; hep-ph/0004189,
c Eq.(17). Here the pole mass actually means the on-shell mass
          if (norder.ge.2)
     >      Factor = Factor + 4d0/3d0* alpi1

          if (norder .eq. 3) then
cpn11 Number of active flavors slightly above the current MSbar
cpn11 mass, according to the current factorization scheme
            NFACTIVE = NFL(QMS+small)
            
            Del=0d0
            DO N2=1,NFACTIVE
              Del = Del + Delta(MSBMASS(N2)/QMS)
            enddo               !N2
c 13.1454= 307d0/32d0 + 2*ZETA2 + 2d0/3d0*dlog(2d0)-1d0/6d0*ZETA3
c -1.04137 = -71d0/144d0 -1d0/3d0 *ZETA2
            Factor = Factor + alpi1**2 * 
     >        (13.1454d0 -1.04137d0* NFACTIVE + 4d0/3d0*Del)
          endif                 !norder.eq.3

          If (norder.gt.3) then
            print *, 
     >        'WARNING: Conversion of the pole mass to MSbar mass ',
     >        'is not implemented for norder > 2'
            stop
          endif 
           
          POLEMASS(N) = QMS*Factor

        ENDIF                   !IMass 2

cpn11   Print MS-bar and pole masses if requested from fitinp
      if(lptp.eq. 1) 
     >    write(*,'(14X,i1,2(3X,F10.5))') N,MSBMASS(N),POLEMASS(N)

      enddo                     !n=1,nftotal
      
      MassesSet=.TRUE.
      
      RETURN
      END !MassSync
