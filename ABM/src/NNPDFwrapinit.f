***********************************************
*
*
*     NNPDFwrapinit.f
*
*     This routine computes, from a grid of NNPDF
*     evolved pdfs, the value of all xpdfs at x and Q
*     (LHAPDF convention) from replica KREP
*     
*     To be used from an external wrapper
*
************************************************

      subroutine InitNNPDFwrap(wrapfile,NMEM)

*     Loads the array of evolved PDF grids for the
*     LHAPDF-like pdf interpolation

      implicit none
      
      integer MXREP
      parameter(MXREP=1e3)
      integer NREP
      common/CNREP/NREP

      integer NPX,NPQ2,NPL,IX,IQ2
      parameter(NPX=60,NPQ2=50)
      parameter(NPL=3000)

      double precision Q2MIN,Q2MAX,XPDFMIN,XPDFMAX,XCH,Q2CH,XCH2
      parameter(Q2MAX=1d8,Q2CH=4d0)
      parameter(XPDFMIN=1d-9,XPDFMAX=1d0,XCH=1D-1,XCH2=0.75D0)
*
      double precision XG(NPX),Q2G(NPQ2),XPDFEV(NPX,NPQ2,-6:6,MXREP)
      common/CPDFGR/XPDFEV,XG,Q2G,Q2MIN,IX,IQ2

      integer IPDF,KREP,NLINESDES,ILINES,NXTMP,NQ2TMP,NMEM
      parameter(NLINESDES=11)
      double precision  XMINTMP,XMAXTMP,Q2MINTMP,Q2MAXTMP

      character*80 wrapfile, info

*     Read the NNPDF wrapper grid file

      write(*,*)"wrapfile=",wrapfile
      open(unit=10,status="old",
     1     file=wrapfile)


*     Read description of grid file
      do ILINES=1,NLINESDES
         read(10,*) info
         write(6,*) info
      enddo
      write(6,*) "  "
         
*     Read max and min extremes of grid
      read(10,*) XMINTMP,XMAXTMP,Q2MINTMP,Q2MAXTMP

*     Read and check grids
*     Read the grid in x
      read(10,*) NXTMP
      if(NXTMP.ne.NPX) then
         write(6,*) "Invalid value of NX!"
         call exit(-10)
      endif
      do IX=1,NPX
         read(10,*) XG(IX)
!         write(6,*) IX,XG(IX)
      enddo
*     Read the grid in Q2
      read(10,*) NQ2TMP
      if(NQ2TMP.ne.NPQ2) then
         write(6,*) "Invalid value of NQ2!"
         call exit(-10)
      endif
      do IQ2=1,NPQ2
         read(10,*) Q2G(IQ2)
!         write(6,*) IQ2,Q2G(IQ2)
*         write(6,*) Q2G(IQ2)
      enddo

*     Read the number of replicas
      read(10,*) NREP
      NMEM=NREP
*     Read the evolved xpdf grid for each replica
      do KREP=1,NREP
         do IX=1,NPX
            do IQ2=1,NPQ2
               read(10,*) ( XPDFEV(IX,IQ2,IPDF,KREP), IPDF=-6,6,1 )
            enddo
         enddo
      enddo
      close(10)
*
      return
      end


* -----------------------------------------------


      subroutine InitNNPDFwrap_mod(wrapfile,NMEM)

*     Loads the array of evolved PDF grids for the
*     LHAPDF-like pdf interpolation

      implicit none
      
      integer MXREP
      parameter(MXREP=1e3)
      integer NREP
      common/CNREP/NREP

      integer NPX,NPQ2,NPL,IX,IQ2,IREP
      parameter(NPX=60,NPQ2=50)
      parameter(NPL=3000)

      double precision Q2MIN,Q2MAX,XPDFMIN,XPDFMAX,XCH,Q2CH,XCH2
      parameter(Q2MAX=1d8,Q2CH=4d0)
      parameter(XPDFMIN=1d-9,XPDFMAX=1d0,XCH=1D-1,XCH2=0.75D0)
*
      double precision XG(NPX),Q2G(NPQ2),XPDFEV(NPX,NPQ2,-6:6,MXREP)
      common/CPDFGR/XPDFEV,XG,Q2G,Q2MIN,IX,IQ2

      integer IPDF,KREP,NLINESDES,ILINES,NXTMP,NQ2TMP,NMEM
      parameter(NLINESDES=11)
      double precision  XMINTMP,XMAXTMP,Q2MINTMP,Q2MAXTMP, ADUM

      character*53 wrapfile, info

*     Read the NNPDF wrapper grid file

      write(*,*)"wrapfile=",wrapfile
      open(unit=10,status="old",
     1     file=wrapfile)


*     Read description of grid file
      do ILINES=1,NLINESDES
         read(10,*) info
         write(6,*) info
      enddo
      write(6,*) "  "
         
*     Read max and min extremes of grid
      read(10,*) XMINTMP,XMAXTMP,Q2MINTMP,Q2MAXTMP

*     Read and check grids
*     Read the grid in x
      read(10,*) NXTMP
      if(NXTMP.ne.NPX) then
         write(6,*) "Invalid value of NX!"
         call exit(-10)
      endif
      do IX=1,NPX
         read(10,*) XG(IX)
!         write(6,*) IX,XG(IX)
      enddo
*     Read the grid in Q2
      read(10,*) NQ2TMP
      if(NQ2TMP.ne.NPQ2) then
         write(6,*) "Invalid value of NQ2!"
         call exit(-10)
      endif
      do IQ2=1,NPQ2
         read(10,*) Q2G(IQ2)
!         write(6,*) IQ2,Q2G(IQ2)
*         write(6,*) Q2G(IQ2)
      enddo

*     Read the number of replicas
      read(10,*) NREP
      NMEM=NREP
      IREP=0
*     Read the evolved xpdf grid for each replica
      do KREP=1,NREP
         if(KREP.lt.500) then
            do IX=1,NPX
               do IQ2=1,NPQ2
                  read(10,*) ( ADUM, IPDF=-6,6,1 )
               enddo
            enddo
         else
            IREP = IREP +1
            write(6,*) "IREP = ",IREP, KREP
            do IX=1,NPX
               do IQ2=1,NPQ2
                  read(10,*) ( XPDFEV(IX,IQ2,IPDF,IREP), IPDF=-6,6,1 )
               enddo
            enddo
         endif
      enddo
      close(10)
*
      return
      end


* -----------------------------------------------
