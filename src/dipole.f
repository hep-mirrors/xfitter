      subroutine DipolePrediction(IDataSet)
C=======================================================
C Calculate dipole model predictions
C 22 Nov 2011 P. Belov (pavel.belov@desy.de)
C=======================================================
      implicit none
C-------------------------------------------------------
      include 'dipole.inc'
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'qcdnumhelper.inc'
      include 'fcn.inc'
C-------------------------------------------------------
      integer IDataSet

      integer idxQ2, idxX, idxY, i,  idx

      integer NPmax
      parameter(NPmax=1000)

      double precision X(NPmax),Y(NPmax),Q2(NPmax)
      double precision yplus, yminus

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

      write(*,*) 'Dipole prediction:: sig0, xlam, x0, xm'
      write(*,*)  sig0, xlam, x0, xm
      write(*,*) ' dipolemodel = ', dipolemodel

C
C define the type of the dipole model
C
      if (DipoleModel.eq.1.or.DipoleModel.eq.3) then
	 icharm = 1
      elseif (DipoleModel.eq.2.or.DipoleModel.eq.4) then
	 icharm = 0
      endif

C define SFUNCT_B inner parameter
      cm = 1.30d0

C
C Get indexes for Q2, x and y bins:
C
      idxQ2 = GetBinIndex(IDataSet,'Q2')
      idxX  = GetBinIndex(IDataSet,'x')
      idxY = GetBinIndex(IDataSet,'y')

      if (idxQ2.eq.0 .or. idxX.eq.0 .or. idxY.eq.0) then
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
         Y(i)   = AbstractBins(idxY,idx)
         Q2(i)  = AbstractBins(idxQ2,idx)

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
C      print *,'> x,Q2,F2,FL,Xsec: ',X(i),Q2(i),F2,FL,
C     $            F2 - Y(i)**2/(1.0d0+(1.0d0-Y(i))**2) * FL
CC-------------------------------------------------

C Get the double differential cross section
         XSec = F2 - Y(i)**2/(1.0d0+(1.0d0-Y(i))**2) * FL

C
C Store cross-section prediction in the global cross-sections table:
C
C        THEO(idx) =  XSec
C
        THEO(idx) = THEO(idx) + XSec

      enddo

      end



      subroutine LeaveOnlyValenceQuarks
C-------------------------------------------------------
C Leave only the contribution of the valence quarks
C for DGLAP+Dipole model fits.
C This subroutine ius called by subroutine fcn
C-------------------------------------------------------
      implicit none

      include 'steering.inc'
      include 'pdfparam.inc'      

      Ag  = 0.0d0
      parglue(1) = 0.0d0
      Aubar = 0.0d0
      parubar(1) = 0.0d0
      Adbar = 0.0d0
      pardbar(1) = 0.0d0
      Asea    = 0.0d0
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
      include 'dipole.inc'
      include 'extrapars.inc'
      include 'steering.inc'

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

      write(*,*) ' <<<the parameters have been read'

      end


      subroutine SetDipoleType
C---------------------------------------
C  Set type of dipole fit.
C  This subroutine ius called by subroutine read_steer.
C---------------------------------------
      implicit none
      include 'steering.inc'
      include 'dipole.inc'
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
      else
         goto 173
      endif

      goto 173
 172   continue
      print '(''Namelist &DipoleType NOT found'')'
      DipoleModel = 0

 173   continue

      end
