! maximum number of data sets and points per set
      integer, parameter :: mset=4, mpint=200
      integer, parameter :: xnd=40, qnd=15, nbk=42

! grid name
      character*100 :: fname

! storing total number of sets
      integer :: cset

! storing reference Lambda scale, no of points, fname per set
      real(8) :: Nscale(mset) ! hardwired currently
      integer :: pset(mset)
      character*100 :: ffnm(mset)
! storing pp or ppbar, and mur/muf
      integer :: sig(mset)
      real(8) :: srof(mset)

! storing reference x-Q-mu scale per point
      real(8) :: xgd(mset,mpint,xnd)
      real(8) :: qgd(mset,mpint,qnd)
      real(8) :: smu(mset,mpint)
      integer :: nxgd(mset,mpint)
      integer :: nqgd(mset,mpint)

! storing grid weight per point
      real(8) :: ciwgt(10,xnd,xnd*qnd,5,mset,mpint)
