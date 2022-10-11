      module inputdata
      USE nrtype
!      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(8)
!     Constants
      real(dp), parameter ::  mu=2.d-5,g=9.81d0,df2=2.d0,df100=100.d0,theta=0.015d0,xacc=1.d-15
      real(dp), parameter ::  denatm=1.225d0,dengas=0.38d0,pn=2.5d0,dz=0.01d0

! extremes for zbrent root search in probfunction subroutine
!      real(dp), parameter :: x1tmp=-50.d0,x2tmp=-1.d-10,x3tmp=1.d-10,x4tmp=50.d0!,tol=1.d-15


! input data common throughout the code
      real(dp) :: dens,dm,dens1,d1mm,sigma1,psi1,nclass1,dens2,d2mm,sigma2,psi2,nclass2,probt,zlam,zlams,c0,ks
      
      end module inputdata
