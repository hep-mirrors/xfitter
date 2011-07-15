      subroutine qcdnum_ini

*     ------------------------------------------------
*     QCDNUM initialisation
*     ------------------------------------------------

      implicit none

      INCLUDE 'steering.inc'
      INCLUDE 'thresholds.inc'
      INCLUDE 'couplings.inc'
      INCLUDE 'datasets.inc'

c set-up of the constants
      integer nc,nf,nfe
      double precision cf,ffnscg,tr
      data nc,nf,nfe,cf,ffnscg,tr /3,6,3,1.33333d0,3.d0,0.5d0/

      double precision zeta2,zeta3,zeta4
      data zeta2 /1.6449340668482264365D0/
      data zeta3 /1.2020569031595942854D0/
      data zeta4 /1.0823232337111381916D0/
  
      integer iord
      integer nxin, nxout, nqin, nqout, iosp
 

C RT parameters:
      double precision alphaS0in,alambdain,flavorin,qsctin,qsdtin
      integer iordin,inullin


      double precision xmin(5)

      integer  iwt(5)

      integer NQGRID, NXGRID

      PARAMETER (NQGRID=11)
      PARAMETER (NXGRID=5)
      

      double precision  QARR(NQGRID),WGT(NQGRID)
      data iosp/2/                   !x grid, lin/qua/spli
c      DATA WGT/2.d0,  1.d0,  1.d0,  1.d0,  1.d0, 1.d0,
      DATA WGT/1.d0,  1.d0,  1.d0,  1.d0, 1.d0,
     $     1.d0, 1.d0, 2.d0, 2.d0, 1.d0, 2.d0/
c      DATA QARR/0.25, 1., 1.1, 1.2,1.6, 1.8, 1.9,
      DATA QARR/1., 1.1, 1.2,1.6, 1.8, 1.9,
     $     1.95,2.89,22.5625,30625.,500000000./
c--   linear x grid
      data xmin/9.9d-7,0.01d0,0.10d0,0.40d0,0.70d0/

      data iwt/1,2,4,8,16/

      double precision qmas, hqmass
      dimension qmas(3), hqmass(3)


      double precision as0,r20
      data as0/0.364/, r20/2.D0/!, nfin/0/ !alphas, NNLO, VFNS
      integer I,ndum,ierr

      double precision a,b
      double precision aq2,bq2,qs0,qt

      integer id1,id2,iq0,iqb,iqc,iqt

      integer nflav
      integer nw,nwords,nx

C Functions:
      integer iqfrmq
      double precision asfunc
C---------------------------------------------------------------------------------------

      q0 = starting_scale
      qc = HF_THRE(1)**2
      qb = HF_THRE(2)**2

      call qcinit(6,' ')        !initialize
      call setord(I_FIT_ORDER)         !LO, NLO, NNLO


      call gxmake(xmin,iwt,nxgrid,200,nx,iosp)        !x-grid
      call gqmake(qarr,wgt,nqgrid,120,nqout)          !mu2-grid
      iq0 =iqfrmq(q0)
      iqc =iqfrmq(qc)
      iqb =iqfrmq(qb)



!SG: hardwire top mass
!      qt = 171.3**2.
      qt = 200.**2
      iqt =iqfrmq(qt)

c      print*,'q0,  qc, qb, qt, iq0,iqc,iqb,iqt', q0,  qc, qb, qt, iq0,iqc,iqb,iqt
c      print*,'vvvv', qfrmiq(iq0), qfrmiq(iqc), qfrmiq(iqb),qfrmiq(iqt)

      if (vfnsINDX.lt.4) then
         call setcbt(0,iqc,iqb,iqt) !thesholds in the ffns
 !! careful in ffns  alphas should be 0.1046         
         call setalf(dble(rtalphas), Mz*Mz) !input alphas
      endif

      call readwt(22,'unpolarised.wgt',id1,id2,nw,ierr)
      if(ierr.ne.0) then
         call fillwt(0,id1,id2,nw) !calculate weights
cv         call dmpwgt(1,22,'unpolarised.wgt')
      else 
         print*,' ERRRROR in read unpolarised wieght', ierr
      endif
      write(6,'(/'' weight: words used ='',I10)') nw     

      call zmreadw(22,'zmstf.wgt',nwords,ierr)
      if(ierr.ne.0) then
         call zmfillw(nwords)
cv         call zmdumpw(22,'zmstf.wgt')
      else 
         print*,' ERRRROR in read zmstf weight', ierr
      endif      
      write(6,'(/'' ZMSTF: words used ='',I10)') nwords      

! setting of the evolution parameters


c      call SETABR(1.D0,0.D0)  ! mur scale variation
c      call ZMDEFQ2(1.D0,0.D0) ! muf scale variation

      if ((mod(HFSCHEME,10).eq.2).or.(vfnsINDX.eq.2)) then
         alambdaIn=0.307
         qs0=1.d0
         alphas0in =asfunc(qs0,nflav,ierr)
         qsdtIn = 4.d0 * qc
         qsctIn = 4.d0 * qb
         flavorIn = 3
         iordIn = 1
         inullIn=0
C-
         call RT_Set_Input(alphaS0in,alambdain,flavorin,qsctin,qsdtin,iordin,inullin)
         call WATE96
      elseif((HFSCHEME.eq.4).or.(vfnsINDX.ge.4)) then

         hqmass(1) = hf_mass(1)
         hqmass(2) = hf_mass(2)
         hqmass(3) = 999.0D0


         call setcbt(3,iqc,iqb,999) !thesholds in the ffns
 !! careful in ffns  alphas should be 0.1046         
         call setalf(dble(rtalphas), Mz*Mz) !input alphas

cv now set the scale
C--- aq=1.0 and bq=0 sets the heavy quarks factorisation scale=Q2
         aq2       = 1.D0
cv         bq2       = 4.d0*qc
         bq2       =  0.D0
C--   Weights HQSTF this maintains an up to date weight file in cwd
         call hqreadw(22,'hqstf.wgt',nwords,ierr)
C--   Check that the parameters read in are those intended      
         if(ierr.eq.0) then
            call hqparms(qmas,a,b)
            if(qmas(1).ne.hqmass(1)) ierr = 1
            if(qmas(2).ne.hqmass(2)) ierr = 1
            if(qmas(3).ne.hqmass(3)) ierr = 1
            if(a.ne.aq2)             ierr = 1
            if(b.ne.bq2)             ierr = 1
            print*,'ERROR!!!!!!!! in hqreadw!', ierr
         endif

         if(ierr.ne.0) then
            print*,'here setting the factorisation scale...'
            call hqfillw(3,hqmass,aq2,bq2,nwords)
cv            call hqdumpw(22,'hqstf.wgt')
         endif      
         write(6,'(/'' HQSTF: words used ='',I10)') nwords      

cv now settings for serghey alechin's code!
cv     set-up of the grids
c         q2min=1.5d0
cv         q2max=1.1d4
c         q2max=1.1d5
c         q2ini=9d0
c         xbmin=1d-5
c         nxpgrid=100
c         nxmgrid=100
c         nspgrid=30
c         nsmgrid=60
c         xbmax=0.999d0
c         x1=0.2
c         xlog1=log(x1)
c         xlog2=log(1-x1)
c         qs0=1.d0


cv this is passed to BMSN code!
c        q20alphas=(91.2d0)**2

cvasfunc(q20alphas,nflav,ierr)
c     c-quark mass
c        rmass(8)=HF_THRE(1)
c     b-quark mass
c        rmass(10)=HF_THRE(2)

      endif
      
      print*,'exit qcdnum_Ini'
      return
      end
