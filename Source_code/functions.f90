      function func(x)
      USE inputdata; USE nrtype
      implicit none
	INTERFACE
		FUNCTION tstud(x)
                USE inputdata; USE nrtype
		IMPLICIT NONE
		REAL(dp) :: tstud
		REAL(dp), INTENT(IN) :: x
		END FUNCTION tstud

		FUNCTION cum(x)
                USE inputdata; USE nrtype
		IMPLICIT NONE
		REAL(dp) :: cum
		REAL(dp), INTENT(IN) :: x
		END FUNCTION cum

		FUNCTION qtrap(a,b)
        	USE nrtype; USE nrutil, ONLY : nrerror
         	REAL(dp), INTENT(IN) :: a,b
        	REAL(dp) :: qtrap
		END FUNCTION qtrap

	END INTERFACE
		REAL(dp), INTENT(IN) :: x
      REAL(dp) :: fx,alfa,usqnrm,d1m,d2m,ushavg,ushmax,ushmin,ztot,ztavg,ztmax,ztmin,zshr,&
     & z0avg,z0max,z0min,pnsavg,pnsmax,pnsmin,z0,pns,px,sigsim,mudstr,mxdstr,mndstr,mm,musim,&
     & pnstmp,func,cd1,ftemp,cavg,cmin,cmax
      real(dp) :: denavg,dfmax,dfmin,uavg,umax,umin,densp
      real(dp) :: t1,t2
      integer nfunc,nfunc1
      common/nf/ nfunc
      common/c1/ fx,alfa,usqnrm
      common/tlaycom/ d1m,d2m,cd1
      common/forprof2/ ushavg,ushmax,ushmin
      common/prof/ densp,ztot,ztavg,ztmax,ztmin,z0avg,z0max,z0min,pnsavg,pnsmax,pnsmin
      common/cnewt/ z0,zshr,pns,pnstmp
      common/prob/ px,sigsim,mudstr,mxdstr,mndstr,mm,musim
      select case (nfunc)
      case (1) !twocomponent
      func=(0.69d0*g*(d1m**3)*x*(1.33d0*dens1-1.33d0*x))/((mu**2)*(((g*(psi1**1.6d0)*&
     &(d1m**3)*x*(dens1-x))/(mu**2))**1.0412))
      case (2) !twolayer
      func=(4.d0*g*(dens2-(theta*g*dm*dens)/(x+theta*g*dm)))/&
     &(3.d0*x*((theta*g*dm*dens)/(x+theta*g*dm)))
      case (3) !twolayer
      func=(0.69d0*g*(d2m**3)*x*(1.33d0*dens2-1.33d0*x))/((mu**2)*&
     &(((g*(psi2**1.6d0)*(d2m**3)*x*(dens2-x))/(mu**2))**1.0412))
      case (4) !testt
      func=alfa*2.d0-tstud(x)
      case (5) !twocomponent
      func=(g*((3.d0*cd1*dens2)/(g*d1m)+(4.d0*(dens2-dens1))/(x)))/(3.d0*dens1)
      case (6) !twocomponent
      func=(0.69d0*g*(d2m**3)*x*(1.33d0*dens2-1.33d0*x))/((mu**2)*&
     &(((g*(psi2**1.6d0)*(d2m**3)*x*(dens2-x))/(mu**2))**1.0412))
      case (7) !twolayer, twocomponent
      ftemp=qtrap(x,usqnrm)
      func=fx-ftemp
      case (8) !twolayer, twocomponent
      ftemp=qtrap(usqnrm,x)
      func=fx-ftemp
      case (9) !funcv
      func=densp*c0*((zlams/(ztot-zlams))*((ztot-x)/x))**pnstmp+&
     &dengas*(1.d0-(c0*((zlams/(ztot-zlams))*((ztot-x)/x))**pnstmp))
      case (10) !funcv
      func=dengas+(densp-dengas)*c0*((z0/(ztot-z0))*((ztot-x)/x))**pns
      case (11) !profiles
      if(x.gt.ztot) then
      cavg=0.d0
      else
      cavg=c0*((z0avg/(ztavg-z0avg))*((ztavg-x)/x))**pnsavg
      endif
      if(x.le.z0avg) cavg=c0
      denavg=cavg*densp+(1.d0-cavg)*dengas
      uavg=ushavg*((1.d0/0.4d0)*log(x/ks)+8.5d0)
      func=0.5d0*denavg*uavg**2
      case (12)  !profiles
      if(x.gt.ztmax) then
      cmin=0.d0
      else
      cmin=c0*((z0min/(ztmax-z0min))*((ztmax-x)/x))**pnsmin
      endif
      if(x.le.z0min) cmin=c0
      dfmin=cmin*densp+(1.d0-cmin)*dengas
      umax=ushmax*((1.d0/0.4d0)*log(x/ks)+8.5d0)
      func=0.5d0*dfmin*umax**2
      case (13)  !profiles
      if(x.gt.ztmin) then
      cmax=0.d0
      else
      cmax=c0*((z0max/(ztmin-z0max))*((ztmin-x)/x))**pnsmax
      endif
      if(x.le.z0max) cmax=c0
      dfmax=cmax*densp+(1.d0-cmax)*dengas
      umin=ushmin*((1.d0/0.4d0)*log(x/ks)+8.5d0)
      func=0.5d0*dfmax*umin**2
      case (14)  !profiles
      if(x.gt.ztot) then
      func=0.d0
      else
      func=c0*((z0avg/(ztavg-z0avg))*((ztavg-x)/x))**pnsavg
      endif
      if(x.le.z0avg) func=c0
      case (15)  !profiles
      if(x.gt.ztmin) then
      func=0.d0
      else
      func=c0*((z0max/(ztmin-z0max))*((ztmin-x)/x))**pnsmax
      endif
      if(x.le.z0max) func=c0
      case (16) !profiles
      if(x.gt.ztmax) then
      func=0.d0
      else
      func=c0*((z0min/(ztmax-z0min))*((ztmax-x)/x))**pnsmin
      endif
      if(x.le.z0min) func=c0
      case (17)  !probfunction
      func=px-cum(x)
      case (18)  !probfunction
      func=mxdstr**x-mudstr**x-(mudstr**x-mndstr**x)
      end select
      end function func

      function func1(x)
      USE inputdata; USE nrtype
      implicit none
      real(dp) :: x,d1m,d2m,func1,cd1
      integer :: nfunc1
      common/nf1/ nfunc1
      common/tlaycom/ d1m,d2m,cd1
      select case (nfunc1)
      case (1) !twolayer
      func1=(4.d0*g*(dens2-(theta*g*dm*dens)/(x+theta*g*dm)))/&
     &(3.d0*x*((theta*g*dm*dens)/(x+theta*g*dm)))
      case (2) !teocomponent
      func1=(g*((3.d0*cd1*dens2)/(g*d1m)+(4.d0*(dens2-dens1))/(x)))/(3.d0*dens1)
      end select
      end function func1

      function cum(x)
      USE nrtype
      implicit none
      REAL(dp) :: cum
      REAL(dp), INTENT(IN) :: x
      cum=0.5d0+0.5d0*erf(x/sqrt(2.d0))
      end function cum

      function dlog2(xx)
      USE nrtype
      implicit none
      REAL(dp) :: dlog2
      REAL(dp), INTENT(IN) :: xx
      dlog2=dlog(xx)/dlog(2.d0)
      end function dlog2

      function rad(xx)
      USE nrtype
      implicit none
      REAL(dp) :: rad
      REAL(dp), INTENT(IN) :: xx
      rad=((2.d0*3.1416d0)/360.d0)*xx
      end function rad

      function grad(xx)
      USE nrtype
      implicit none
      REAL(dp) :: grad
      REAL(dp), INTENT(IN) :: xx
      grad=(360.d0/(2.d0*3.1416d0))*xx
      end function grad

