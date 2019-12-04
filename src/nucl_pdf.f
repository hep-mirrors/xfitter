C 07/02/2017 created by Marina Walt, University of Tuebingen
C module to summarize all modification to xFitter meant to implement the treatment of nuclear PDFs

C in order to be able to call the new routines from the old subprograms of xFitter
C nucl_pdf.f needs to be added to the Makefile in /src folder 
C (and commands make clean, make, install to be executed after the modification of the Makefile).

C -------------------------------------------------------- 
C Test subroutine
C -------------------------------------------------------- 

       subroutine nucl_pdf
       
       implicit none
       
       print*,'Testing nucl_pdf.f'
       
       end


C --------------------------------------------------------    
C Subroutine to set nPDFType
C -------------------------------------------------------- 

C Subroutine to set the PDFType (incl. new value 'nuclear').
C The functions of this new subroutine are meant to replace the former Subroutine SetPDFType().
C For that, it is necessary to modify the existing xFitter source code file 'read_steer.f' as follows.
C First, modify main steering parameters namelist, by adding the new parameters:
C       namelist/xFitter/
C     $     Znucleus,Anucleus	         ! new parameters for nPDFType 'nucleus'
C Second, modify the subroutine:
C       Subroutine SetPDFType()
C       call SetnPDFType()      
C       end

C In order to activate the treatment of flexible input parameters A, Z in the existing xFitter code,
C the code in the file 'interface/src/hf_pdf_calls.f' needs to be modified as fallows:
C      data A,Z /207,82/		! comment out
C ! Old program code
C      if(lead) then
C         if (deuteron) then		! comment out
C            A=2			! comment out
C            Z=1			! comment out
C         endif				! comment out
C ! New program code - set nuclear numbers A,Z to the input parameters from the steering.txt file
C      if(nucleus) then
C         A = Anucleus
C         Z = Znucleus
C ----
C            for fitting of nuclear PDFs with nuclear data, division by A is NOT required, but needs to be disabled.
C ! Old program code
C         tmpU  = (Z*PDFSF( 1) + (A-Z)*PDFSF( 2) )/A	
C         tmpD  = (Z*PDFSF( 2) + (A-Z)*PDFSF( 1) )/A	
C         tmpUb = (Z*PDFSF(-1) + (A-Z)*PDFSF(-2) )/A	
C         tmpDb = (Z*PDFSF(-2) + (A-Z)*PDFSF(-1) )/A	
C ! New program code
C         tmpU  = (Z*PDFSF( 1) + (A-Z)*PDFSF( 2) )
C         tmpD  = (Z*PDFSF( 2) + (A-Z)*PDFSF( 1) )
C         tmpUb = (Z*PDFSF(-1) + (A-Z)*PDFSF(-2) )
C         tmpDb = (Z*PDFSF(-2) + (A-Z)*PDFSF(-1) )  

       Subroutine SetnPDFType()
       
C     07/02/2017 new parameters introduced by Marina Walt, university of Tuebingen 
C     (include/steering.inc needs to be modified as follows)
C      logical nucleus             		!> Flag to trigger nuclear PDF	
C     --------
C     ! modify namelist 
C      common/STEERING/
C     &     Q2VAL,starting_scale,strange_frac, Chi2MaxError,
C     ...
C     $     ,npolyval, lead, deuteron, nucleus !> add deuteron and nucleus to common steering namelist
C     ---------
C      double precision Znucleus		!> proton number Z for PDFType 'nucleus'
C      double precision Anucleus		!> nucleon number A for PDFType 'nucleus'
C      common/CPdfStyle/PDFStyle,PDFType,
C     $     Znucleus,Anucleus     
       
       implicit none
#include "steering.inc"
     
       double precision Ainput, Zinput
       
C -- part the of original source code             
       if (PDFType.eq.'proton'.or. PDFType.eq.'PROTON') then
        nucleus = .false.
        print *,'Fitting for PROTON PDFs, PDFType=', PDFType
       elseif (PDFType.eq.'lead'.or. PDFType.eq.'LEAD') then
         lead = .true. 
         deuteron = .false.
         print *,'Fitting for LEAD PDFs, PDFType=', PDFType
       elseif (PDFType.eq.'DEUTERON'.or. PDFType.eq.'deuteron') then
         lead = .true. 
         deuteron = .true. 
         print *,'Fitting for DEUTERON PDFs, PDFType=', PDFType 
C -- part the of original source code          


       elseif (PDFType.eq.'nucleus'.or. PDFType.eq.'NUCLEUS') then
        nucleus = .true.
        if (Anucleus.eq.0) then
         print '(''Error: nucleon number Anucleus is missing in the steering file'')'
         call HF_stop
        else
         Ainput = Anucleus
        endif
        if (Znucleus.eq.0) then
         print '(''Error: proton number Znucleus is missing in the steering file'')'
C         call HF_stop		! Z-number can be blank, if one wants to treat a neutron (A=1, Z=0)
        else
         Zinput = Znucleus
        endif
        print *,'Fitting for NUCLEAR PDFs, PDFType=', PDFType

       else
        call hf_errlog(300920131,
     $   'F: Unsupported PDFType used!')
       endif
       
       end
       
       
C -------------------------------------------------------- 
C Subroutine to read nuclear data
C 21/02/2017 created by Marina Walt, University of Tuebingen
C -------------------------------------------------------- 

C additional input data are necessary for treatment of nPDFs

       Subroutine ReadnuclDataFile

C       double precision function ReadnuclDataFile()

C to test the new subroutine, the existing subroutine 'ReadDataFile(CFile)' in 'read_data.f' needs to be modified as follows
C ! existing code
C      DATASETInfoDimension(NDATASETS) = NInfo		! existing code
C      do i=1,NInfo					! existing code
C         DATASETInfoNames(i,NDATASETS) = CInfo(i)	! existing code
C         DATASETInfo(i,NDATASETS) =      DataInfo(i)	! existing code
C      enddo						! existing code
C ! new code line          
C         call ReadnuclDataFile			

       
       implicit none       
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "systematics.inc"
#include "theorexpr.inc"
#include "scales.inc"
#include "couplings.inc"
#include "for_debug.inc"

       integer i,k
       integer IDataSet

C   19/11/2018 xfitter-2.0.0        
C       parameter (ninfoMax=100)  !>PARAMETER attribute of 'ninfomax' conflicts with PARAMETER attribute at (1)
C       integer ninfomax !> Symbol 'ninfomax' already has basic type of INTEGER
C       integer  NInfo   !> Symbol 'ninfo' already has basic type of INTEGER
C       double precision DataInfo(ninfoMax)  !> Symbol 'datainfo' already has basic type of REAL
C       character *80 CInfo(ninfoMax) !>Symbol 'cinfo' already has basic type of CHARACTER
       
C output:
       double precision A1,Z1,A2,Z2

C functions:       
       integer GetInfoIndex
       
C Namelist definition:
       namelist/Data/
     $     NInfo,DataInfo,CInfo
          
          
C For reference purpose only...
C Extra info:
C       DATASETInfoDimension(NDATASETS) = NInfo
C       do i=1,NInfo
C         DATASETInfoNames(i,NDATASETS) = CInfo(i)
C         DATASETInfo(i,NDATASETS) =      DataInfo(i)
C       enddo      

C ------- 
       IDataSet = NDATASETS

C      set nuclear number A (variable called 'A1')      
       if(DATASETInfo( GetInfoIndex(IDataSet, 'A1'), IDataSet).gt.0) then
        A1=DATASETInfo( GetInfoIndex(IDataSet, 'A1'), IDataSet)
C ----- debugging
        if (Debug) then
          print *,'A1 provided for NDATASETS'
          print *,NDATASETS
          print *,'is'
          print *,A1
        endif
C -----         
       else
C ----- set to default (proton)        
        A1=1.0
C ----- debugging
        if (Debug) then
          print *,'A1 NOT provided for NDATASETS'
          print *,NDATASETS
          print *,'default value set'
          print *,A1
        endif  
       endif     

C      set proton number Z (variable called 'Z1')
       if(DATASETInfo( GetInfoIndex(IDataSet, 'Z1'), IDataSet).gt.0) then
        Z1=DATASETInfo( GetInfoIndex(IDataSet, 'Z1'), IDataSet)
C ----- debugging
        if (Debug) then
          print *,'Z1 provided for NDATASETS'
          print *,NDATASETS
          print *,'is'
          print *,Z1
        endif
       else
C ----- set to default (proton)        
        Z1=1.0
C ----- debugging
        if (Debug) then
          print *,'Z1 NOT provided for NDATASETS'
          print *,NDATASETS
          print *,'default value set'
          print *,Z1
        endif
       endif        

C      OPTIONAL: set nuclear number A of the second nucleus (variable called 'A2')        
       if(DATASETInfo( GetInfoIndex(IDataSet, 'A2'), IDataSet).gt.0) then
        A2=DATASETInfo( GetInfoIndex(IDataSet, 'A2'), IDataSet)
C ----- debugging
        if (Debug) then
          print *,'A2 provided for NDATASETS'
          print *,NDATASETS
          print *,'is'
          print *,A2
        endif
       endif
       
C      OPTIONAL: set proton number Z of the second nucleus (variable called 'Z2')         
       if(DATASETInfo( GetInfoIndex(IDataSet, 'Z2'), IDataSet).gt.0) then
        Z2=DATASETInfo( GetInfoIndex(IDataSet, 'Z2'), IDataSet)
C ----- debugging
        if (Debug) then
          print *,'Z2 provided for NDATASETS'
          print *,NDATASETS
          print *,'is'
          print *,Z2
        endif
       endif        
       
C       return
       
       end       

       
C --------------------------------------------------------
C Subroutine to consider sumrules for a nucleus with CTEQ parametrisation
C 15/03/2017 created by Marina Walt, University of Tuebingen (modified copy of function SumRuleCTEQ)
C --------------------------------------------------------

      double precision function SumRuleCTEQ_nucl(n,acteq)
C---------------------------------------------------------------
C Sum-rule integral for 
C  
C UF = a0*E**(a3*x)*(1 - x)**a2*x**(a1 + n)*(1 + E**a4*x + E**a5*x**2)
C
C parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
C---------------------------------------------------------------
      implicit none
      
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "pdfparam.inc"
#include "for_debug.inc"

      double precision Anucl1, Anucl2, Anucl
      integer k,j, kmax      
      
      integer n
      double precision acteq(1:27), atempcteq(1:27), acteqinput(1:27)
      double precision YF_nucl, temp
      double precision DGammF,HypG1F1r, HypG1F1,SumRulenTUJU
      
C functions:       
      integer GetInfoIndex
       
C C++ functions:
      external testintegration
      external integratecncteq
      double precision xmin, xmax, res, err
       
C-----------------------------------------------------
      
C      print *,'debugging SumRulesCTEQ_nucl, before if-loop: acteq'
C      print *,acteq
      
      acteqinput = acteq
      
      if (nucleus) then
       Anucl1=DATASETInfo( GetInfoIndex(idxADataSet, 'A1'), idxADataSet)
       
       if (ratio.eq.1.0) then		! value of the global parameter ratio is set inside the FCN routine
        if (ratiostep.eq.1) then
         Anucl=Anucl1
        elseif (ratiostep.eq.2) then
         Anucl2=DATASETInfo( GetInfoIndex(idxADataSet, 'A2'), idxADataSet)
         Anucl=Anucl2
        else
         print '(''Error: ratiostep not defined!'')'
         call HF_stop
        endif
       else
        Anucl=Anucl1
       endif
       

       if (Anucl.gt.0) then
       
C --   10/10/2017
        if (nCTEQframework) then
         kmax=6
        elseif (nTUJU) then
         kmax=5
        else
         kmax=9
        endif 
C --            
       
        do k=2,kmax
         j=2*k-1
         atempcteq(k) = acteq(k)+acteq(j+9)*(1-Anucl**((-1)*acteq(j+9+1)))
         acteq(k) = atempcteq(k)
        enddo
       endif
      endif
      
      xmin = 0.0
      xmax = 1.0
      
      if (nCTEQ) then
        call integratecncteq(n,acteq, xmin, xmax, res, err)
        SumRuleCTEQ_nucl = res
        if (Debug) then
          print *,'c++ integration result'
          print *,res
        endif  
      
      elseif (nTUJU) then
C        numerical integration routine (C++)      
C        call integratecnTUJU(n,acteq, xmin, xmax, res, err)
C        SumRuleCTEQ_nucl = res

C        For debugging purpose only, to compare with analytical integration routine
C        temp = SumRulenTUJU(n,acteq)
C        print *,'check analytical integration,num-res,an-res',SumRuleCTEQ_nucl,temp

C        18/04/2018 Marina Walt, University of Tuebingen
C        analytical integration routine, (x=1).
        SumRuleCTEQ_nucl = SumRulenTUJU(n,acteq)

      
      else

        YF_nucl = (
     &  DGammF(1 + acteq(3) )*DGammF(1 + acteq(2) + n)*(
     &   Hypg1F1(1 + acteq(2) + n,2 + acteq(2) + acteq(3) + n,
     $     acteq(4)) + 
     &   (1 + acteq(2) + n)*DGammF(2 + acteq(2) + acteq(3) + n)*
     &   ( exp(acteq(5))*
     &    Hypg1F1R(2 + acteq(2) + n,3 + acteq(2)
     $     + acteq(3) + n,acteq(4)) + 
     &     exp(acteq(6))*(2 + acteq(2) + n)*
     &    Hypg1F1R(3 + acteq(2) + n,4 + acteq(2)
     $     + acteq(3) + n,acteq(4))
     &   )
     &  )
     & )/DGammF(2 + acteq(2) + acteq(3) + n)

        SumRuleCTEQ_nucl = YF_nucl
      endif
      
      acteq = acteqinput

      end
   

   
C---------------------------------------------------------------
C 18/04/2018 Marina Walt, University of Tuebingen
C Sum-rule integral for 
C UF = a1**x**(a2 + n)*(1 - x)**a3*(1 + a4*x + a5*x**2)
C parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
C---------------------------------------------------------------   
   
   
      double precision function SumRulenTUJU(n,acteq)

      implicit none 
      integer n
      double precision acteq(1:5)
      double precision YF_nucl
      double precision DGammF
C-----------------------------------------------------
      YF_nucl = (
     &  DGammF(1 + acteq(3) )*DGammF(2 + acteq(2) + n)/((1+acteq(2)+n)*DGammF(2+acteq(2)+acteq(3)+n))
     $     + acteq(4)*DGammF(1 + acteq(3) )*DGammF(3 + acteq(2) + n)/((2+acteq(2)+n)*DGammF(3+acteq(2)+acteq(3)+n))
     $     + acteq(5)*DGammF(1 + acteq(3) )*DGammF(4 + acteq(2) + n)/((3+acteq(2)+n)*DGammF(4+acteq(2)+acteq(3)+n))  
     &  )

      SumRulenTUJU = YF_nucl

      end   
   
   
C -------------------------------------------------------- 
C Subroutine to store calculated cross sections in a local file
C 10/04/2017 created by Marina Walt, University of Tuebingen
C -------------------------------------------------------- 

      subroutine WriteCSforDebugging

      implicit none

#include "steering.inc"
#include "pdfparam.inc"
#include "ntot.inc"
#include "datasets.inc"
#include "indata.inc"
#include "systematics.inc"
#include "theo.inc"

      integer i,j,index,k,reacindx
      
      double precision currEcharge
      
      integer idxQ2, idxX, idxY, idxSigma

C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      
      character*48 name
      character*3 charI
      

      do i=1,ndatasets
          idxX = GetBinIndex(i,'x')
          idxY = GetBinIndex(i,'y')
          idxQ2 = GetBinIndex(i,'Q2') 
          reacindx = 0
          currEcharge = 0
          
          write(charI,'(i0)'), i
          
          name='./'//TRIM(OutDirName)//'/crosssections_'//TRIM(charI)//'.dat'
            
            
          open(90,file=name)
             
          write(90,*) '  #ColumnName = ''reaction index'', ''x'', ''Q2'', ''Daten(x-section)'',' //
     & ' ''THEO '' '
          
          do j=1,NDATAPOINTS(i)
             index = DATASETIDX(i,j)
                
             write(90,'(1X,i5,1X,4(e11.5,1X),i4)') 
     $              i,
     $              AbstractBins(idxX,index),
     $              AbstractBins(idxQ2,index),
     $              DATEN(index),
     $              THEO(index)
             
          enddo
          close(90)
          
      enddo
  111  format(1X, F10.3, 2X, F12.6, 2X, 3(F12.6,2X))
      close(90)

C      RETURN
      end subroutine
       
       

      double precision function nCTEQparaDU(x,a)
C----------------------------------------------------
C
C 16/08/2017 Marina Walt, University of Tuebingen
C nCTEQ like parametrisation of Dbar/Ubar as used in nCTEQ15 Paper 'PRD93, 085037 (2016)', equation (2.5)
C
C-----------------------------------------------------
      implicit none
      
C ---------------------------    
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "pdfparam.inc"
#include "for_debug.inc"

      double precision Anucl1,Anucl2,Anucl
      integer k,j      

C functions:       
      integer GetInfoIndex
C ---------------------------             
      
      double precision x,a(1:27)
      double precision UF      
      
      
      UF = a(1)*x**a(2)*(1-x)**a(3) + (1+a(4)*x)*(1-x)**a(5)
       

      nCTEQparaDU = UF

      end

      subroutine HF_Get_PDFs_nucl(x,q2,PDFSF)
C----------------------------------------------------------------------
C Interface to PDF 
C
C  Input:  x, Q2 values
C  Output: 13 PDF values
C----------------------------------------------------------------------
      implicit none
C----------------------------------------------------------------------
#include "steering.inc"
#include "fcn.inc"
#include "ntot.inc"
#include "datasets.inc"
#include "indata.inc"
#include "pdfparam.inc"

      double precision x,q2
     $     ,pdfsf(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF)

      integer i
      double precision Anucl,Znucl, tmpU,tmpD,tmpUb,tmpDb
C     13/02/2017 modified by Marina Walt, University of Tuebingen
C      data A,Z /207,82/	! comment out

      double precision FSNSXQ
      
C -----------------------------   
C 10/03/2017
       integer IDataSet

C functions:       
       integer GetInfoIndex 
      
C----------------------------------------------------------------------
      if (PDFStyle.eq.'LHAPDFNATIVE') then
         if (ExtraPdfs) then
C photon is present !
            call evolvePDFphoton(x, sqrt(q2), pdfsf, pdfsf(7))
         else
            call evolvePDF(x, sqrt(q2), pdfsf)
         endif
         return
      endif

      if ( ExtraPdfs ) then
C QED evolution:
C         call FPDFXQ(iPDFSET,x,q2,PDFSF,ICheck_QCDNUM)
         call ALLFXQ(iPDFSET,x,q2,PDFSF,0,ICHECK_QCDNUM)         
         PDFSF(N_CHARGE_PDF+N_NEUTRAL_PDF) = FSNSXQ( iPDFSET,13,x,q2,ICheck_QCDNUM)
c         print *,ipdfset,x,q2, PDFSF(N_CHARGE_PDF+1), PDFSF(0)
         return
      endif

      if (x.ge.1.D0) then
         do i=-6,6
            PDFSF(i) = 0
         enddo

         if (ICheck_QCDNUM.gt.0) then
            Call HF_errlog(12042201,'W:HF_GET_PDFS X value >= 1.0') 
         endif
         Return
      endif

      if ( CachePDFs ) then
C Cache PDF calls:
         call GetCachedPDFs(ifcncount,x,q2,PDFSF)
      else
C Get PDFs directly:
         Call HF_Get_PDFs_UnCached(x,q2,PDFSF)
      endif
  
C---- switch for lead PDF: Combine to form nuclear pdf; scale by A
C----     For full cross setion on lead, multiply by A

C     27/03/2017 modified by Marina Walt, University of Tuebingen
C                modification to activate the treatment of flexible, dataset dependent input parameters A,Z.

      if(nucleus) then
C --- 30/03/2017 Marina Walt, modification to consider experimental data provided for a ratio sigma(A1)/sigma(A2)
         print *,'debug hf_pdf_calls_nucl.f'

         if (ratio.eq.1.0) then			! value of the global parameter ratio is set inside the FCN routine
          if (ratiostep.eq.1) then
           Anucl=DATASETInfo( GetInfoIndex(idxADataSet, 'A1'), idxADataSet)
           Znucl=DATASETInfo( GetInfoIndex(idxADataSet, 'Z1'), idxADataSet)
          elseif (ratiostep.eq.2) then
           Anucl=DATASETInfo( GetInfoIndex(idxADataSet, 'A2'), idxADataSet)
           Znucl=DATASETInfo( GetInfoIndex(idxADataSet, 'Z2'), idxADataSet)
          else
           print '(''hf_pfd_calls Error: ratiostep not defined!'')'
           call HF_stop
          endif
         else
          Anucl=DATASETInfo( GetInfoIndex(idxADataSet, 'A1'), idxADataSet)
          Znucl=DATASETInfo( GetInfoIndex(idxADataSet, 'Z1'), idxADataSet)
          if (Anucl.eq.0) then 
           Anucl = Anucleus         
           Znucl = Znucleus
          endif
         endif
C ---         
         tmpU  = (Znucl*PDFSF( 1) + (Anucl-Znucl)*PDFSF( 2) )/Anucl
         tmpD  = (Znucl*PDFSF( 2) + (Anucl-Znucl)*PDFSF( 1) )/Anucl
         tmpUb = (Znucl*PDFSF(-1) + (Anucl-Znucl)*PDFSF(-2) )/Anucl
         tmpDb = (Znucl*PDFSF(-2) + (Anucl-Znucl)*PDFSF(-1) )/Anucl
         
         PDFSF( 1) = tmpU
         PDFSF( 2) = tmpD
         PDFSF(-1) = tmpUb
         PDFSF(-2) = tmpDb       
      endif

C----------------------------------------------------------------------
      end

      
      
C>!   24/10/2018 Marina Walt, University of Tuebingen
C                New function to calculate isocorrections on the calculated Xsection in order to reflect exp data
      double precision function isocorrectionfunc(IDataSet,x,Q2)

      implicit none
#include "ntot.inc"
#include "couplings.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "qcdnumhelper.inc"
#include "pdfparam.inc"
#include "for_debug.inc"
#include "steering.inc"
      
      integer IDataSet
    
      double precision x,Q2
      
      double precision EMCflag, SLACflag, NMCflag             
      double precision NMCA, NMCB, F2nF2p     
C Functions:
      integer GetInfoIndex

            
C ------ 
            
              
      if (ratiostep.eq.1) then
                   ADIS = DATASETInfo( GetInfoIndex(IDataSet, 'A1'), IDataSet)
                   ZDIS = DATASETInfo( GetInfoIndex(IDataSet, 'Z1'), IDataSet)
      elseif (ratiostep.eq.2) then
                   ADIS = DATASETInfo( GetInfoIndex(IDataSet, 'A2'), IDataSet)
                   ZDIS = DATASETInfo( GetInfoIndex(IDataSet, 'Z2'), IDataSet)
      endif

      SLACflag = DATASETInfo( GetInfoIndex(IDataSet, 'SLAC'), IDataSet)
      EMCflag = DATASETInfo( GetInfoIndex(IDataSet, 'EMC'), IDataSet)
      NMCflag = DATASETInfo( GetInfoIndex(IDataSet, 'NMC'), IDataSet)
              
      if (SLACflag.eq.1.0) then     !>modifications to consider exp SLAC data provided with isoscalar corrections, as per arXiv:1612.05741v1
                   isocorrectionfunc = (ADIS/2.0)*(1+(1-0.8*x))/(ZDIS+(ADIS-ZDIS)*(1-0.8*x))
      elseif (EMCflag.eq.1.0) then
                   isocorrectionfunc = (ADIS/2.0)*(1+(0.92-0.86*x))/(ZDIS+(ADIS-ZDIS)*(0.92-0.86*x))
      elseif (NMCflag.eq.1.0) then          !> isoscalar correction as per paper from Amaudruz et al. (NMC), Nucl. Phys. B 371 (1992) 3, Appendix B)
                   NMCA=0.979 - 1.692*x+2.797*x**2-4.313*x**3+3.075*x**4
                   NMCB=-0.171*x+0.244*x**2
                   F2nF2p=NMCA*(Q2/20)**NMCB*(1+((x**2)/Q2))              
                   isocorrectionfunc = (ADIS/2.0)*(1+F2nF2p)/(ZDIS+(ADIS-ZDIS)*(F2nF2p)) 
       else
                   isocorrectionfunc = 1.0                   
       endif

C ------ 24/10/2017  
           
      end

      
      
