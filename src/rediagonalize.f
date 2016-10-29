C---------------------------------------------
!> Re-diagonalize PDF eigenvectors 
!> @param Ndata number of data points
!> @param Nsyst number of systematic sources
C---------------------------------------------
      Subroutine rediagonalize(Ndata,Nsyst)
      implicit none      
#include "ntot.inc"
#include "steering.inc"
C#include "datasets.inc"
#include "indata.inc"
#include "systematics.inc"
#include "covar.inc"
#include "theo.inc"

      integer Ndata,NSyst
      integer i,j,k,nsysloc
      double precision theo_err2_loc(Ndata)
      double precision Eigenvalues(Nsyst)  
      double precision RotBeta(Nsyst,Ndata)  ! dynamic 
      double precision totpdf(Ndata),totpdftest(Ndata)
c      double precision, allocatable :: test(:,:)

      logical lreset 
C-------------------------------------

      nsysloc = nsys*2 - nsys  ! that is ugly fix of the Fortran optimization problem
      print *,'sys',nsysloc
C      print *,'ndata=',ndata
      call read_outdirnml
      call read_lhapdfnml

      lreset = .false.
      do k=1,ndata
         theo_err2_loc(k) = (theo_fix(k)/100.)**2
     $        *(theo_stat(k)**2+theo_unc(k)**2)

c         print *,'ho',sqrt(theo_err2_loc(k))
         if (theo_err2_loc(k).eq.0) then
            call hf_errlog(13122014,
     $           'W: Rediagonalize: Recet weights to unity')
            lreset = .true.
         endif
      enddo
      if (lreset) then
         do k=1,ndata
            theo_err2_loc(k) = 1.0
         enddo
      endif

      do i=1,nsyst
         do j=1,nsyst
            cov(i,j) = 0.
            if (i.eq.j) then
               cov(i,j) = 1.0
            endif
            do k=1,ndata
               cov(i,j) = cov(i,j) + beta(i,k)*beta(j,k)
     $              *theo_fix(k)**2/10000./theo_err2_loc(k)
            enddo
         enddo
      enddo

      call MyDSYEVD(Nsyst,Cov,NTot,Eigenvalues)

      print '(''Eigenvalues:'')'
      do i=nsyst,1,-1
         print '(i4,F10.2)',i,Eigenvalues(i)
      enddo

      print '(''Creating pdf_rotation.dat file ...'')'

      open (51,file=trim(OutDirName)//'/pdf_rotation.dat'
     $     ,status='unknown')
      write (51,'(''LHAPDF set='',A32)') 
     $     trim(adjustl(LHAPDFSET))
      
      write (51,'(i4)') NSyst
      do i=NSyst,1,-1
         write (51,'(i5,200F10.6)') i,( Cov(j,i),j=1,NSyst)
      enddo
      close (51)

      print *,' '
      print '(''Fraction of uncertainty^2 in rotated vectors'')'
      do k=1,ndata
         totpdf(k) = 0.
         totpdftest(k) = 0.
         do i=1,nsyst
            rotbeta(i,k) = 0.
            do j=1,nsyst
               rotbeta(i,k) = rotbeta(i,k) + 
     $              beta(j,k)*cov(j,i)               
            enddo            
            totpdf(k) = totpdf(k) + beta(i,k)**2
            totpdftest(k) = totpdftest(k) + rotbeta(i,k)**2
         enddo
         totpdf(k) = sqrt(totpdf(k))
         totpdftest(k) = sqrt(totpdftest(k))

c         print *,totpdf(k),totpdftest(k)

         print '(i4, 100F5.1)',k,(rotbeta(i,k)**2
     $        /totpdf(k)**2*100.,i=nsyst,1,-1)
      enddo
      print *,' '
      

      do k=1,ndata
         do j=1,nsyst
            beta(nsyst-j+1,k) = rotbeta(j,k)/100.
         enddo
      enddo
c      print *,nsyst,nsysloc
      nsys = 0
c      print *,nsyst,nsysloc      

      call WriteTheoryFiles(nsysloc,theo_fix,.true.)

C-------------------------------------
      end
