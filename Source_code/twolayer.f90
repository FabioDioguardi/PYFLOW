      subroutine twolayer
      USE inputdata; USE nrtype
      implicit none
	INTERFACE

		FUNCTION func(x)
                USE inputdata; USE nrtype
		IMPLICIT NONE
		REAL(dp) :: func
                REAL(dp), INTENT(IN) :: x
		END FUNCTION func

		FUNCTION dlog2(xx)
                USE inputdata; USE nrtype
		IMPLICIT NONE
		REAL(dp) :: dlog2
		REAL(dp), INTENT(IN) :: xx
		END FUNCTION dlog2

		FUNCTION zbrent(x1,x2)
                USE inputdata; USE nrtype; USE nrutil, ONLY : nrerror
                IMPLICIT NONE
	        REAL(dp), INTENT(IN) :: x1,x2
         	REAL(dp) :: zbrent
         	END FUNCTION zbrent

		FUNCTION qtrap(a,b)
        	USE nrtype; USE nrutil, ONLY : nrerror
         	REAL(dp), INTENT(IN) :: a,b
        	REAL(dp) :: qtrap
		END FUNCTION qtrap

       	END INTERFACE
      character(len=20) :: nmodel,cmd1
      REAL(dp) :: fx,d1m,d2m,alfa,usqnrm,tcalc,ttab,dennrm,denmax,denmin,tauavg,taumax,taumin,ushavg,ushmax,ushmin
      REAL(dp) :: d2,cdd2,cdd100,usq2,usq100,cddavg,cdavg,dmodm,dmod,denmod,cddnrm,aleft,aright,cddlft,cddrgt,&
     &usqmin,usqmax,cddmin,cddmax,nfreed,nclass,cd1,cd2,s
      integer :: nfunc,nfunc1,ntest,nrest
      common/c1/ fx,alfa,usqnrm
      common/tlaycom/ d1m,d2m,cd1
      common/tres/ tcalc,ttab
      common/tinp/ nfreed,nclass
      common/nf/ nfunc
      common/nf1/ nfunc1
      common/forprof1/ dennrm,denmax,denmin,tauavg,taumax,taumin
      common/forprof2/ ushavg,ushmax,ushmin
      nmodel='Two layer method'
      nclass=nclass2
      alfa=probt/2.d0
      d2=-dlog2(d2mm)
      d2m=d2mm/1000.d0
      cdd2=(4.d0*(dens2-df2))/(dm*theta*3.d0*(dens-df2))                !Cd/d for flow density = 2 kg/m^3
      cdd100=(4.d0*(dens2-df100))/(dm*theta*3.d0*(dens-df100))          !Cd/d for flow density = 100 kg/m^3
      write(52,350)nmodel,cdd2,cdd100
      write(*,350)nmodel,cdd2,cdd100
      fx=cdd2
      usq2=-(4.d0*g*dm*theta*(dens-dens2))/(3.d0*fx*theta*dm*dens-4.d0*dens2)     !ushear^2 for flow density = 2 kg/m^3
      fx=cdd100
      usq100=-(4.d0*g*dm*theta*(dens-dens2))/(3.d0*fx*theta*dm*dens-4.d0*dens2)   !ushear^2 for flow density = 100 kg/m^3
      write(52,351)usq2,usq100
      write(*,351)usq2,usq100
!     Cd/davg
      nfunc=2
      s=qtrap(usq100,usq2)
      cddavg=(1.d0/(usq2-usq100))*s
!     Cdavg
      nfunc=3
      s=qtrap(df2,df100)
      cdavg=(1.d0/(df100-df2))*s
      write(52,352)cddavg,cdavg
      write(*,352)cddavg,cdavg
!     Model diameter (m and phi units)
      dmodm=cdavg/cddavg
      dmod=-dlog2(dmodm*1000.d0)
      write(52,353)dmodm,dmod
      write(*,353)dmodm,dmod
!     t-Student test for comparing dmod and d
  301 call testt(d2,dmod,sigma2,nrest)
      write(*,303)ttab,tcalc
  303 format(/,'t tabulated =',f7.3,2x,'t calculated =',f7.3,/)
      if(nrest.eq.2) then
  302 write(*,*)'Failed t-test! Stop the calculations?'
      write(*,*)'1: Stop the calculations 2:Reduce significance level'
      read(*,*)ntest
      select case (ntest)
      case (1)
      stop
      case (2)
      write(*,*)'Write the new significance level'
      read(*,*)probt
      alfa=probt/2.d0
      goto 301
      case default
      write(*,*)'Wrong choice'
      goto 302
      end select
      else
      write(*,*)'Test t OK'
      endif
!     Model flow density
      denmod=(3.d0*cddavg*theta*dm*dens-4.d0*dens2)/(3.d0*theta*dm*cddavg-4.d0)
!     Normalized flow density
      dennrm=(d2m*dens2)/((dmodm*(dens2-denmod))/denmod+d2m)
!     Normalized Cd/d
      cddnrm=(4.d0*(dens2-dennrm))/(dm*theta*3.d0*(dens-dennrm))
      write(52,355)denmod,dennrm,cddnrm
      write(*,355)denmod,dennrm,cddnrm
!     Normalized ushear^2
      fx=cddnrm
      usqnrm=-(4.d0*g*dm*theta*(dens-dens2))/(3.d0*fx*theta*dm*dens-4.d0*dens2)
      write(52,356)usqnrm
      write(*,356)usqnrm
      if(usqnrm.ge.usq2.or.usqnrm.le.usq100) then
      write(*,*)'Warning! ush^2 norm. is not in the range ush^2 (2 kg/m^3) - ush^2(100 kg/m^3)'
      write(*,*)'This is due to the very low accuracy level of the T - Student test'
      write(52,*)'Warning! ush^2 norm. is not in the range ush^2 (2 kg/m^3) - ush^2(100 kg/m^3)'
      write(52,*)'This is due to the very low accuracy level of the T - Student test'
      write(*,*)'Type "stop" if you want to stop the calculations'
      write(*,*)'Any other input will continue the calculations'
      read(*,*)cmd1
      if(cmd1.eq.'stop') stop
      endif
!     Cd/d(34% left), Cd/d(34% right)
      nfunc=2
      aleft=qtrap(usq100,usqnrm)
      cddlft=0.68d0*aleft
      aright=qtrap(usqnrm,usq2)
      cddrgt=0.68d0*aright
      write(52,357)aleft,cddlft,aright,cddrgt
      write(*,357)aleft,cddlft,aright,cddrgt
!     ushear^2 max and min
      nfunc1=1
      nfunc=7
      fx=cddlft
      write(52,*)'ush^2 min calculation residuals'
      write(*,*)'ush^2 min calculation residuals'
      usqmin=zbrent(1.d-5,600.d0)
      nfunc=8
      fx=cddrgt
      write(52,*)''
      write(52,*)'ush^2 max calculation residuals'
      write(*,*)''
      write(*,*)'ush^2 max calculation residuals'
      usqmax=zbrent(1.d-5,600.d0)
      write(52,358)usqmin,usqmax
      write(*,358)usqmin,usqmax
!     Cd/d(ushear^2 max) Cd/d(ushear^2 min)
      nfunc=2
      cddmin=func(usqmin)
      cddmax=func(usqmax)
      write(52,359)cddmin,cddmax
      write(*,359)cddmin,cddmax
!     Flow densities as a function of Cd/d(ushear^2 max) Cd/d(ushear^2 min)
      denmax=(3.d0*cddmin*theta*dm*dens-4.d0*dens2)/(3.d0*theta*dm*cddmin-4.d0)
      denmin=(3.d0*cddmax*theta*dm*dens-4.d0*dens2)/(3.d0*theta*dm*cddmax-4.d0)
      write(52,360)denmax,denmin
      write(*,360)denmax,denmin
!     Shear stresses
      tauavg=dennrm*usqnrm
      taumin=denmin*usqmax
      taumax=denmax*usqmin
!     Shear velocities
      ushavg=sqrt(usqnrm)
      ushmin=sqrt(usqmin)
      ushmax=sqrt(usqmax)
  350 format(a20,/,'Cd/d(2 kg/m^3) =',f12.4,1x,'Cd/d(100 kg/m^3) =',f12.4,/)
  351 format('ush^2(2kg/m^3) =',f8.4,1x,'ush^2(100 kg/m^3) =',f8.4,/)
  352 format('Cd/d avg =',f12.4,1x,'Cd avg =',f10.4,/)
  353 format('dmod (m) =',f8.4,1x,'dmod (phi) =',f8.4,/)
  355 format('Denmod =',f8.4,1x,'Den norm =',f8.4,1x,'Cd/d norm =',f10.4,/)
  356 format('ush^2 norm =',f8.4,/)
  357 format('Aleft =',f10.4,1x,'Cd/d left =',f10.4,1x,'Aright =',f10.4,1x,'Cd/d right =',f10.4,/)
  358 format('ush^2 min =',f8.4,1x,'ush^2 max =',f8.4,/)
  359 format('Cd/d(ush^2 min) =',f12.4,1x,'Cd/d(ush^2 max) =',f12.4,/)
  360 format('Den(Cd/d(ush^2 min)) =',f8.4,1x,'Den(Cd/d(ush^2 max)) =',f8.4,/)
      end subroutine twolayer
