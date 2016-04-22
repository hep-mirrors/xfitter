      FUNCTION ALAMBD(AMU)
C                                                   -=-=- alambd
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Qcdpar
C Description: Many callable functions to return QCD parameters
C ListOfFiles: alambd nfl alamf amass amhatf nfltot amumin ch 
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Qcdpar.f,v 1.1 97/12/21 20:34:52 wkt Exp $ 
C #Log:	Qcdpar.f,v $
c Revision 1.1  97/12/21  20:34:52  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      ALAMBD = ALAMF(NFL(AMU))

      RETURN
      END
C
C***************************************************************
C
