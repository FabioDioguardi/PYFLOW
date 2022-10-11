      subroutine profiles
      USE inputdata; USE nrtype
      implicit none
	INTERFACE
                FUNCTION func(x)
                USE inputdata; USE nrtype
		IMPLICIT NONE
		REAL(dp) :: func
		REAL(dp), INTENT(IN) :: x
		END FUNCTION func

                FUNCTION grad(xx)
                USE inputdata; USE nrtype
                implicit none
		REAL(dp) :: grad
                real(dp), INTENT(IN) :: xx
		END FUNCTION grad

                FUNCTION rad(xx)
                USE inputdata; USE nrtype
                implicit none
		REAL(dp) :: rad
                real(dp), INTENT(IN) :: xx
		END FUNCTION rad

		FUNCTION qtrap(a,b)
        	USE nrtype; USE nrutil, ONLY : nrerror
         	REAL(dp), INTENT(IN) :: a,b
        	REAL(dp) :: qtrap
		END FUNCTION qtrap
		
	END INTERFACE
      real(dp), dimension(20) :: zdynpr,pzav1,pzmax1,pzmin1,pzavg,&
      pzmax,pzmin,czav1,czmax1,czmin1,czavg,czmax,czmin,zc
      real(dp), dimension(2) :: x
      real(dp), dimension(3) :: arr
      real(dp) :: dennrm,denmax,denmin,tauavg,taumax,taumin,ushavg,ushmax,ushmin
      real(dp) :: z,ztot,ztavg,ztmax,ztmin,z0avg,z0max,z0min,pnsavg,pnsmax,pnsmin
      real(dp) :: z0,zshr,pns,pnstmp
      real(dp) :: p10av1,p10mx1,p10mn1,zsfavg,zsfmax,zsfmin,c2av1,c2max1,c2min1,&
      p10avg,p10max,p10min,c2avg,c2max,c2min
      real(dp) :: den,cavg,cmax,cmin,denavg,densp,dfmax,dfmin,pdynav,pdynmx,pdynmn,&
     &s,senx,sloapp,uavg,umax,umin
      integer :: nfunc,nfunc1,nnewt,i,j,ipr,ic,ndynpr,nc
      LOGICAL check
      character(42) :: avgcon,maxcon,mincon
      character(37) :: avgpd,maxpd,minpd
      character(30) :: avgvel,maxvel,minvel
      character(33) :: avgden,maxden,minden
      character(10) :: height,cmd
      common/nf/ nfunc
      common/nf1/ nfunc1
      common/forprof1/ dennrm,denmax,denmin,tauavg,taumax,taumin
      common/forprof2/ ushavg,ushmax,ushmin
      common/prof/ densp,ztot,ztavg,ztmax,ztmin,z0avg,z0max,z0min,pnsavg,pnsmax,pnsmin
      common/prof1/ p10av1,p10mx1,p10mn1,c2av1,c2max1,c2min1,pzav1,pzmax1,pzmin1,czav1,czmax1,czmin1
      common/prof2/ zsfavg,zsfmax,zsfmin
      common/prof3/ zdynpr,p10avg,p10max,p10min,c2avg,c2max,c2min,pzavg,pzmax,pzmin,&
     &czavg,czmax,czmin,zc,ipr,ic,ndynpr,nc
      common/cnewt/ z0,zshr,pns,pnstmp
      common/cnewt2/ den,nnewt
      avgcon='50th percentile particle concentration (-)'
      maxcon='84th percentile particle concentration (-)'
      mincon='16th percentile particle concentration (-)'
      avgpd='50th percentile dynamic pressure (Pa)'
      maxpd='84th percentile dynamic pressure (Pa)'
      minpd='16th percentile dynamic pressure (Pa)'
      avgvel='50th percentile velocity (m/s)'
      maxvel='84th percentile velocity (m/s)'
      minvel='16th percentile velocity (m/s)'
      avgden='50th percentile density (kg/m^3)'
      maxden='84th percentile density (kg/m^3)'
      minden='16th percentile density (kg/m^3)'
      height='Height (m)'      
!     Average. maximum and minimum shear flow particle concentration
      densp=dens2
      cavg=(dennrm-dengas)/(densp-dengas)
      cmax=(denmax-dengas)/(densp-dengas)
      cmin=(denmin-dengas)/(densp-dengas)
      write(52,363)cavg,cmax,cmin
      write(*,363)cavg,cmax,cmin
  363 format('C 50th =',e12.5,1x,'C 84th =',e12.5,1x,'C 16th =',e12.5,/)
!     Average, maximum and minimum shear flow thickness
      ztavg=zlam/cavg
      ztmax=zlam/cmin
      ztmin=zlam/cmax
      write(52,364)ztavg,ztmax,ztmin
      write(*,364)ztavg,ztmax,ztmin
  364 format('Htot 50th =',f8.3,1x,'Htot 84th =',f8.3,1x,'Htot 16th =',f8.3,/)
      nnewt=1
 !    Pnsusp,avg and slope angle
      ztot=ztavg
      den=dennrm
      x(1)=1.d0
      x(2)=1.d0
      write(52,*)'Pns avg and zsf avg calculation residuals'
      write(*,*)'Pns avg and zsf avg calculation residuals'
      call newt(x,check)
      pnsavg=x(1)
      zsfavg=x(2)
      z0avg=zlams
      senx=tauavg/((den-denatm)*g*zsfavg)
      sloapp=grad(asin(senx))
      write(52,366)pnsavg,z0avg,sloapp
      write(*,366)pnsavg,z0avg,sloapp
  366 format('Pnsusp avg =',f6.3,1x,'z0avg =',f9.6,1x,'Slope (ø) =',f6.3,/)
      nnewt=2
!     Pnsusp,max and z0max
      den=denmax
      zsfmin=taumin/((den-denatm)*g*sin(rad(sloapp)))
      zshr=zsfmin
      ztot=ztmin
      x(1)=0.1d0
      x(2)=0.00001d0
      write(52,*)'Pns max and z0 max calculation residuals'
      write(*,*)'Pns max and z0 max calculation residuals'
      call newt(x,check)
      pnsmax=x(1)
      z0max=x(2)
      write(52,367)pnsmax,z0max
      write(*,367)pnsmax,z0max
  367 format('Pnsusp max =',f6.3,1x,'z0max =',f9.6,/)
!     Pnsusp,min and z0min
      den=denmin
      zsfmax=taumax/((den-denatm)*g*sin(rad(sloapp)))
      zshr=zsfmax
      ztot=ztmax
      x(1)=0.1d0
      x(2)=0.0001d0
      write(52,*)'Pns min and z0 min calculation residuals'
      write(*,*)'Pns min and z0 min calculation residuals'
      call newt(x,check)
      pnsmin=x(1)
      z0min=x(2)
      write(52,368)pnsmin,z0min
      write(*,368)pnsmin,z0min
  368 format('Pnsusp min =',f6.3,1x,'z0min =',f9.6,/)
!     PROFILES
      z=zlam
!     Average profiles
      write(53,170)height,avgcon,maxcon,mincon
      write(54,171)height,avgpd,maxpd,minpd
      write(55,172)height,avgvel,maxvel,minvel
      write(56,173)height,avgden,maxden,minden
  170 format(a10,3(2x,a42))
  171 format(a10,3(2x,a37))
  172 format(a10,3(2x,a30))
  173 format(a10,3(2x,a33))
  102 if(z.gt.ztavg) then
      cavg=0.d0
      else
      cavg=c0*((z0avg/(ztavg-z0avg))*((ztavg-z)/z))**pnsavg
      endif
      if(z.le.z0avg) cavg=c0
      denavg=cavg*densp+(1.d0-cavg)*dengas
      uavg=ushavg*((1.d0/0.4d0)*log(z/ks)+8.5d0)
      pdynav=0.5d0*denavg*uavg**2
!     Maximum profiles
      if(z.gt.ztmin) then
      cmax=0.d0
      else
      cmax=c0*((z0max/(ztmin-z0max))*((ztmin-z)/z))**pnsmax
      endif
      if(z.le.z0max) cmax=c0
      dfmax=cmax*densp+(1.d0-cmax)*dengas
      umax=ushmax*((1.d0/0.4d0)*log(z/ks)+8.5d0)
!     Minimum profiles
      if(z.gt.ztmax) then
      cmin=0.d0
      else
      cmin=c0*((z0min/(ztmax-z0min))*((ztmax-z)/z))**pnsmin
      endif
      if(z.le.z0min) cmin=c0
      dfmin=cmin*densp+(1.d0-cmin)*dengas
      umin=ushmin*((1.d0/0.4d0)*log(z/ks)+8.5d0)
      pdynmx=0.5d0*dfmin*umax**2
      pdynmn=0.5d0*dfmax*umin**2
      write(53,174)z,cavg,cmax,cmin
      write(54,175)z,pdynav,pdynmx,pdynmn
      write(55,176)z,uavg,umax,umin
      write(56,177)z,denavg,dfmax,dfmin
  174 format(1x,f8.3,12x,e12.5,2(33x,e12.5))
  175 format(1x,f8.3,12x,e12.4,2(27x,e12.4))
  176 format(1x,f8.3,12x,f7.3,2(26x,f7.3))
  177 format(1x,f8.3,12x,f8.3,2(27x,f8.3))
      z=z+dz
      if(z.gt.ztavg) goto 101
      goto 102
!     Specific heights dynamic pressure
      write(*,*)'####################################################'
  101 write(*,*)'Average dynamic pressure will be calculated for 10 m'
  500 write(*,*)'Do you need to calculate over other heights?'
      write(*,*)'1:Yes 2:No'
      read(*,*)ndynpr
      if(ndynpr.ne.1.and.ndynpr.ne.2) then
      write(*,*)'Wrong choice'
      goto 500
      else
      if(ndynpr.eq.1) then
      write(*,*)'Write below the heights for the dynamic pressure'
      do i=1,1000
  505 write(*,*)'z (m) = ?'
      read(*,*)zdynpr(i)
      if(zdynpr(i).le.zlam) then
      write(*,*)'Warning. z < layer thickness'
      goto 505
      else
      ipr=i
      write(*,*)'Type "stop" if you want to end the height input'
      write(*,*)'Type anything else if you want to continue the input'
      read(*,*)cmd
      if(cmd.eq.'stop') goto 501
      endif
      enddo
      endif
      endif
  501 cmd='ciao'
      nfunc=11
      s=qtrap(z0avg,10.d0)
      p10av1=(1.d0/(10.d0-z0avg))*s
      do j=1,ipr
      s=qtrap(z0avg,zdynpr(j))
      pzav1(j)=(1.d0/(zdynpr(j)-z0avg))*s
      enddo
      nfunc=12
      s=qtrap(z0min,10.d0)
      p10mx1=(1.d0/(10.d0-z0min))*s
      do j=1,ipr
      s=qtrap(z0min,zdynpr(j))
      pzmax1(j)=(1.d0/(zdynpr(j)-z0min))*s
      enddo
      nfunc=13
      s=qtrap(z0max,10.d0)
      p10mn1=(1.d0/(10.d0-z0max))*s
      do j=1,ipr
      s=qtrap(z0max,zdynpr(j))
      pzmin1(j)=(1.d0/(zdynpr(j)-z0max))*s
      enddo
!     Particle concentration at 2 m and at user requested heights
      write(*,*)'Particle concentration will be calculated for z=2 m'
  502 write(*,*)'Do you need to calculate over other heights?'
      write(*,*)'1:Yes 2:No'
      read(*,*)nc
      if(nc.ne.1.and.nc.ne.2) then
      write(*,*)'Wrong choice'
      goto 502
      else
      if(nc.eq.1) then
      write(*,*)'Write below the heights for the particle concentration'
      do i=1,1000
  504 write(*,*)'z (m) = ?'
      read(*,*)zc(i)
      if(zc(i).le.zlam) then
      write(*,*)'Warning. z < layer thickness'
      goto 504
      else
      ic=i
      write(*,*)'Type "stop" if you want to end the height input'
      write(*,*)'Type anything else if you want to continue the input'
      read(*,*)cmd
      if(cmd.eq.'stop') goto 503
      endif
      enddo
      endif
      endif
  503 nfunc=14
      c2av1=func(2.d0)
      do j=1,ic
      czav1(j)=func(zc(j))
      enddo
      nfunc=15
      c2max1=func(2.d0)
      do j=1,ic
      czmax1(j)=func(zc(j))
      enddo
      nfunc=16
      c2min1=func(2.d0)
      do j=1,ic
      czmin1(j)=func(zc(j))
      enddo
!     Real average, maximum and minimum solutions
      arr(1)=p10av1
      arr(2)=p10mx1
      arr(3)=p10mn1
      call piksrt(3,arr)
      p10max=arr(3)
      p10avg=arr(2)
      p10min=arr(1)
      do j=1,ipr
      arr(1)=pzav1(j)
      arr(2)=pzmax1(j)
      arr(3)=pzmin1(j)
      call piksrt(3,arr)
      pzmax(j)=arr(3)
      pzavg(j)=arr(2)
      pzmin(j)=arr(1)
      enddo
      arr(1)=c2av1
      arr(2)=c2max1
      arr(3)=c2min1
      call piksrt(3,arr)
      c2max=arr(3)
      c2avg=arr(2)
      c2min=arr(1)
      do j=1,ic
      arr(1)=czav1(j)
      arr(2)=czmax1(j)
      arr(3)=czmin1(j)
      call piksrt(3,arr)
      czmax(j)=arr(3)
      czavg(j)=arr(2)
      czmin(j)=arr(1)
      enddo
      end subroutine profiles
      
      SUBROUTINE piksrt(n,arr)
      USE nrtype
      implicit none
      INTEGER :: n,i,j
      real(dp), dimension(n) :: arr
      real(dp) :: a
      do j=2,n
      a=arr(j)
      do i=j-1,1,-1
      if(arr(i).le.a) goto 10
      arr(i+1)=arr(i)
      enddo
      i=0
   10 arr(i+1)=a
      enddo
      end subroutine piksrt
