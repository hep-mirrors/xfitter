      subroutine DipolePrediction(IDataSet)
C=======================================================
C Calculate dipole model predictions
C 22 Nov 2011 P. Belov (pavel.belov@desy.de)
C=======================================================
      implicit none
C-------------------------------------------------------
#include "dipole.inc"
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "qcdnumhelper.inc"
#include "fcn.inc"
C-------------------------------------------------------
      integer IDataSet

      integer idxQ2, idxX, idxY, i,  idx, idxS

      integer NPmax
      parameter(NPmax=1000)

      double precision X(NPmax),Y(NPmax),Q2(NPmax)
      double precision yplus, yminus, S

      double precision F2,FL,XSec

      double precision
     & ft_l,
     & fl_l,
     & ft_c,
     & fl_c,
     $ temp_l,
     $ temp_c

C Functions:
      integer GetBinIndex
      integer GetInfoIndex

      write(*,*) 'Dipole prediction:: sig0, xlam, x0, xm, cBGK, eBGK'
      write(*,*)  sig0, xlam, x0, xm, cBGK, eBGK
      write(*,*) ' dipolemodel = ', dipolemodel

C
C define the type of the dipole model
C
      if (DipoleModel.eq.1.or.DipoleModel.eq.3.or.DipoleModel.eq.5) then
	 icharm = 1
      elseif (DipoleModel.eq.2.or.DipoleModel.eq.4) then
	 icharm = 0
      endif

C define SFUNCT_B inner parameter: charm mass
      cm = 1.30d0  ! please change if needed

C
C Get indexes for Q2, x and y bins:
C
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
      idxY = GetBinIndex(IDataSet,'y')
      idxS =  GetInfoIndex(IDataSet,'sqrt(S)')

      if (idxY.eq.0) then
         idxS =  GetInfoIndex(IDataSet,'sqrt(S)')
         if (idxS .gt. 0) then
            S = (DATASETInfo( GetInfoIndex(IDataSet,'sqrt(S)')
     $           , IDataSet))**2
         else
            print *,
     $ 'ERROR: DIS sample, neigher S nor y are defined !'
            print *,' !!! STOP STOP STOP STOP !!!'
            call HF_stop
         endif
      endif

      if (idxQ2.eq.0 .or. idxX.eq.0) then
         Return
      endif


C calculate predictions for each bin:
      do i=1,NDATAPOINTS(IDataSet)
C
C Reference from the dataset to a global data index:
C
         idx =  DATASETIDX(IDataSet,i)
C
C Local X,Y,Q2 arrays, used for caclulations:
C
         X(i)   = AbstractBins(idxX,idx)
         Q2(i)  = AbstractBins(idxQ2,idx)
         if (idxY.eq.0) then
            Y(i)   = Q2(i) / ( X(i) * S )
         else
            Y(i)   = AbstractBins(idxY,idx)
         endif



C
C Call the calculating subroutine
C
         CALL SFUNCT_B(X(i),Q2(i),ft_l,fl_l,ft_c,fl_c)

C         write(*,*) 'debug output: ',ft_l,fl_l,ft_c,fl_c

C Calculate temporary variables
	 temp_l = (ft_l + fl_l)
	 temp_c = (ft_c + fl_c)

C Get F2 and FL (we assume that R=0)
	 F2 = temp_l + temp_c
	 FL = (fl_l + fl_c) ! extra factor is absent

CC----------- output of the interim results -------
c      print *,'> x,Q2,F2,FL,Xsec: ',X(i),Q2(i),F2,FL,
c     $            F2 - Y(i)**2/(1.0d0+(1.0d0-Y(i))**2) * FL
CC-------------------------------------------------

C Get the double differential cross section
         XSec = F2 - Y(i)**2/(1.0d0+(1.0d0-Y(i))**2) * FL

C
C Store cross-section prediction in the global cross-sections table:
C
C        THEO(idx) =  XSec
C
        THEO(idx) = THEO(idx) + XSec


        ! COMMENT OUT
c        print *,'ss2',i,x(i),Q2(i),f2,fl


      enddo

      end



      subroutine LeaveOnlyValenceQuarks
C-------------------------------------------------------
C Leave only the contribution of the valence quarks
C for DGLAP+Dipole model fits.
C This subroutine ius called by subroutine fcn
C-------------------------------------------------------
      implicit none

#include "steering.inc"
#include "pdfparam.inc"

      parglue(1) = 0.0d0
      parubar(1) = 0.0d0
      pardbar(1) = 0.0d0
      parsea(1) = 0.0d0
      
      end


      subroutine RemoveOnlySeaQuarks
C-------------------------------------------------------
C Remove only sea quarks 09/05/2012
C 
C This subroutine ius called by subroutine fcn
C-------------------------------------------------------
      implicit none

#include "steering.inc"
#include "pdfparam.inc"

C      parglue(1) = 0.0d0
C      parubar(1) = 0.0d0
C      pardbar(1) = 0.0d0
      write(*,*) 'Remove sea quarks'
      parsea(1) = 0.0d0
      
      end



      subroutine DecodeDipolePar(pars)
C-------------------------------------------------------
C Read parameters for the dipole models from par array.
C This subroutine ius called by subroutine PDF_param_iteration.
C 22 Nov 2011 P.Belov
C upd: 01/12/2011
C-------------------------------------------------------
      implicit none 
#include "dipole.inc"
#include "extrapars.inc"
#include "steering.inc"

      integer idx
      integer GetParameterIndex
      double precision pars(*)

      write(*,*) ' read parameters for dipole models>>>'

      idx = GetParameterIndex('SIG0')
      idx = iExtraParamMinuit(idx)
      sig0 = pars(idx)

      idx = GetParameterIndex('XLAM')
      idx = iExtraParamMinuit(idx)
      xlam = pars(idx)
      
      idx = GetParameterIndex('X0')
      idx = iExtraParamMinuit(idx)
      x0   = pars(idx)
      
      idx = GetParameterIndex('XM')
      idx = iExtraParamMinuit(idx)
      xm   = pars(idx)

      idx = GetParameterIndex('cBGK')
      idx = iExtraParamMinuit(idx)
      cBGK   = pars(idx)

      idx = GetParameterIndex('eBGK')
      idx = iExtraParamMinuit(idx)
      eBGK   = pars(idx)


      write(*,*) ' <<<the parameters have been read'

      end


      subroutine SetDipoleType
C---------------------------------------
C  Set type of dipole fit.
C  This subroutine ius called by subroutine read_steer.
C---------------------------------------
      implicit none
#include "steering.inc"
#include "dipole.inc"
C---------------------------------

      open (51,file='steering.txt',status='old')
      read (51,NML=DipoleType,END=171,ERR=172)
 171   continue      
      close (51)

      if (TypeOfDipoleFit.eq.'NONE') then
         DipoleModel = 0
      elseif (TypeOfDipoleFit.eq.'GBW') then
         DipoleModel = 1
      elseif (TypeOfDipoleFit.eq.'IIM') then
         DipoleModel = 2
      elseif (TypeOfDipoleFit.eq.'DGLAP+GBW') then
         DipoleModel = 3
      elseif (TypeOfDipoleFit.eq.'DGLAP+IIM') then
         DipoleModel = 4
      elseif (TypeOfDipoleFit.eq.'BGK') then
         DipoleModel = 5
      else
         goto 173
      endif

      goto 173
 172   continue
      print '(''Namelist &DipoleType NOT found'')'
      DipoleModel = 0

 173   continue

      end

C 7/06/2016 A.L: fits with and without saturation
            subroutine SetDipCsType
C---------------------------------------
C Set type of dipole fit.
C This subroutine is called by subroutine read_steer.
C---------------------------------------
      implicit none
#include "steering.inc"
#include "dipole.inc"
C---------------------------------

C Set default to avoid extra printout:
      TypeOfDipCs = 'linear'

      open (51,file='steering.txt',status='old')
      read (51,NML=DipCsType,END=181,ERR=182)
 181   continue      
      close (51)
      

      if (TypeOfDipCs.eq.'saturation') then
         DipCsModel = 1
      elseif (TypeOfDipCs.eq.'linear') then
         DipCsModel = 0
      else
         write(*,*) 'TypeOfDipCs parameter not recognised. Using linear.'
         DipCsModel = 0
         goto 183
      endif

      goto 183
 182   continue
      print '(''Namelist &DipCsType NOT found'')'
      DipCsModel = 0

 183   continue

      end
c
c      
      subroutine dipoleBGK(IDataSet)
C---------------------------------------------------
C
C Dipole model which uses gluon from QCDNUM
C and valence quarks from DGLAP
C
C---------------------------------------------------
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      integer IDataSet
      integer kflag
C----------------------------

C Addition of valence quars:
      call SaveRestorePdfs(0)   ! save all info

            
      call SaveRestorePdfs(1)   ! reset sea, gluon
      kflag = 0
C      call SumRules(kflag)
      call Evolution            ! evolve PDFs without sea, gluon.
C     Do dis:
      
      print *,'gluon',parglue(1),parglue(2),parglue(3)

      print *,'Start DGLAP'
      Call GetNCXsection(IDataSet, HFSCHEME)

      call SaveRestorePdfs(2)   ! store them back
      call Evolution            ! evolve PDFs with gluon

c            print *,'Start dipole model:'

      call DipolePrediction(IDataSet)
*     write(*,*) 'tu jestem'
*     stop

                                ! Restore all PDFs
      call SaveRestorePdfs(3)
      end
