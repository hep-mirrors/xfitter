      subroutine qcdnum_ini

*     ------------------------------------------------
*     QCDNUM initialisation
*     ------------------------------------------------

      implicit double precision(a-h,o-z)

      INCLUDE 'steering.inc'
      INCLUDE 'thresholds.inc'
      INCLUDE 'couplings.inc'
      INCLUDE 'datasets.inc'
cv FFNS
c      INCLUDE 'CONSTCOM.'
c      INCLUDE 'APSCOM5.'
c set-up of the constants
      data nc,nf,nfe,cf,ffnscg,tr /3,6,3,1.33333d0,3.d0,0.5d0/
      data zeta2 /1.6449340668482264365D0/
      data zeta3 /1.2020569031595942854D0/
      data zeta4 /1.0823232337111381916D0/
  
      integer iord
      integer nxin, nxout, nqin, nqout, iosp
 

      COMMON/INPUT/alphaS0,alambda,flavor,qsct,qsdt,iord,inull
      dimension xmin(5)
      integer  iwt(5)
      PARAMETER (NQGRID=10)
      PARAMETER (NXGRID=5)

      DIMENSION QARR(NQGRID),WGT(NQGRID)
      data iosp/2/                   !x grid, lin/qua/spli
c      DATA WGT/2.d0,  1.d0,  1.d0,  1.d0,  1.d0, 1.d0,
      DATA WGT/1.d0,  1.d0,  1.d0,  1.d0, 1.d0,
     $     1.d0, 2.d0, 2.d0, 1.d0, 2.d0/
c      DATA QARR/0.25, 1., 1.1, 1.2,1.6, 1.8, 1.9,
      DATA QARR/1., 1.1, 1.2,1.6, 1.8, 1.9,
     $     1.95,22.5625,30625.,500000000./
c--   linear x grid
      data xmin/9.9d-7,0.01d0,0.10d0,0.40d0,0.70d0/

      data iwt/1,2,4,8,16/

      double precision qmas, hqmass
      dimension qmas(3), hqmass(3)


      data as0/0.364/, r20/2.D0/!, nfin/0/ !alphas, NNLO, VFNS
      integer I,ndum,ierr



      q0 = starting_scale
      qc = HF_THRE(1)**2
      qb = HF_THRE(2)**2

      call qcinit(6,' ')        !initialize
      call setord(I_FIT_ORDER)         !LO, NLO, NNLO


      call gxmake(xmin,iwt,nxgrid,200,nx,iosp)        !x-grid
      call gqmake(qarr,wgt,nqgrid,120,nqout)              !mu2-grid
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
      write(lunout,'(/'' weight: words used ='',I10)') nw     

      call zmreadw(22,'zmstf.wgt',nwords,ierr)
      if(ierr.ne.0) then
         call zmfillw(nwords)
cv         call zmdumpw(22,'zmstf.wgt')
      else 
         print*,' ERRRROR in read zmstf wieght', ierr
      endif      
      write(lunout,'(/'' ZMSTF: words used ='',I10)') nwords      

! setting of the evolution parameters


c      call SETABR(1.D0,0.D0)  ! mur scale variation
c      call ZMDEFQ2(1.D0,0.D0) ! muf scale variation

      if ((HFSCHEME.eq.2).or.(vfnsINDX.eq.2)) then
         alambda=0.307
         qs0=1.d0
         alphas0=asfunc(qs0,nflav,ierr)
         qsdt = 4.d0 * qc
         qsct = 4.d0 * qb
         flavor = 3
         iord = 1
         inull=0
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
            print*,'ERRRRRRRROR!!!!!!!! in hqreadw!', ierr
         endif

         if(ierr.ne.0) then
            print*,'here setting the factorisation scale...'
            call hqfillw(3,hqmass,aq2,bq2,nwords)
cv            call hqdumpw(22,'hqstf.wgt')
         endif      
         write(lunout,'(/'' HQSTF: words used ='',I10)') nwords      

cv now settings for serghey alechin's code!
cv     set-up of the grids
         q2min=1.5d0
cv         q2max=1.1d4
         q2max=1.1d5
         q2ini=9d0
         xbmin=1d-5
         nxpgrid=100
         nxmgrid=100
         nspgrid=30
         nsmgrid=60
         xbmax=0.999d0
         x1=0.2
         xlog1=log(x1)
         xlog2=log(1-x1)
         qs0=1.d0


cv this is passed to BMSN code!
        q20alphas=(91.2d0)**2
        falphas0=0.1135d0

cvasfunc(q20alphas,nflav,ierr)
c     c-quark mass
c        rmass(8)=HF_THRE(1)
c     b-quark mass
c        rmass(10)=HF_THRE(2)

      endif
      
      print*,'exit qcdnum_Ini'
      return
      end
