      USE inputdata; USE nrtype
      implicit none
      real(dp), dimension(20) :: zdynpr,pzav1,pzmax1,pzmin1,pzavg,pzmax,pzmin,&
     &czav1,czmax1,czmin1,czavg,czmax,czmin,zc
      real(dp) :: dennrm,denmax,denmin,tauavg,taumax,taumin,ushavg,ushmax,ushmin,fx,usqnrm
      real(dp) :: p10av1,p10mx1,p10mn1,c2av1,c2max1,c2min1
      real(dp) :: zsfavg,zsfmax,zsfmin,p10avg,p10max,p10min,c2avg,c2max,c2min
      real(dp) :: alfa,dentmp,dtemp,sigtmp,psitmp,ncltmp
      real(dp) :: tcalc,ttab
      real(dp) :: densp,ztot,ztavg,ztmax,ztmin,z0avg,z0max,z0min,pnsavg,pnsmax,pnsmin
      integer :: i,ipr,ic,ndynpr,nc,nmod,ninp,ndata
      common/c1/ fx,alfa,usqnrm
      common/tres/ tcalc,ttab
      common/forprof/ dennrm,denmax,denmin,tauavg,taumax,taumin,ushavg,ushmax,ushmin
      common/prof/ densp,ztot,ztavg,ztmax,ztmin,z0avg,z0max,z0min,pnsavg,pnsmax,pnsmin
      common/prof1/ p10av1,p10mx1,p10mn1,c2av1,c2max1,c2min1,pzav1,pzmax1,pzmin1,czav1,czmax1,czmin1
      common/prof2/ zsfavg,zsfmax,zsfmin
      common/prof3/ zdynpr,p10avg,p10max,p10min,c2avg,c2max,c2min,pzavg,pzmax,pzmin&
     &,czavg,czmax,czmin,zc,ipr,ic,ndynpr,nc
      write(*,*)'###PROGRAM PYFLOW (2013) by Dioguardi Fabio###'
      write(*,*)'###Based on Dellino et al. (2008) model ####'
      write(*,*)''
   20 write(*,*)'Choose a method'
      write(*,*)'1: Two layer method'
      write(*,*)'2: Two components method'
      read(*,*)nmod                                                     ! Model choice
      if(nmod.ne.1.and.nmod.ne.2) then
      write(*,*)'Wrong choice'
      goto 20
      else
  120 write(*,*)'Input data'                                            ! How to provide input data: input file or keyboard
      write(*,*)'1:File input 2:keyboard'
      read(*,*)ninp
      select case (ninp)
      case (1)
      open(51,file='input')
      read(51,*)dens,dm,dens1,d1mm,sigma1,psi1,nclass1,dens2,d2mm,sigma2,& ! input data
     &psi2,nclass2,probt,zlam,zlams,c0,ks
      case (2)
      if(nmod.eq.1) then
      write(*,*)'Density of the entrained particle (kg/m^3)'
      read(*,*)dens
      write(*,*)'Diameter of the entrained particle (m)'
      read(*,*)dm
      write(*,*)'Particle density (kg/m^3)'
      read(*,*)dens2
      write(*,*)'Particle equivalent diameter of the median size (mm)'
      read(*,*)d2mm
      write(*,*)'Sorting particle grainsize distribution (phi)'
      read(*,*)sigma2
      write(*,*)'Particle shape factor (-)'
      read(*,*)psi2
      write(*,*)'Classes number particle grainsize distribution (-)'
      read(*,*)nclass2
      write(*,*)'Significance level t-test'
      read(*,*)probt
      write(*,*)'Layer thickness (m)'
      read(*,*)zlam
      write(*,*)'Sublayer thickness (m) (write 0 if it is not known)'
      read(*,*)zlams
      write(*,*)'Particle concentration in the layer (-)'
      read(*,*)c0
      write(*,*)'Substrate roughness (m)'
      read(*,*)ks
      else
      write(*,*)'Density particle 1 (kg/m^3)'
      read(*,*)dens1
      write(*,*)'Particle 1 equivalent diameter of the median size (mm)'
      read(*,*)d1mm
      write(*,*)'Sorting particle grainsize distribution 1 (phi)'
      read(*,*)sigma1
      write(*,*)'Particle 1 shape factor (-)'
      read(*,*)psi1
      write(*,*)'Classes number particle grainsize distribution 1 (-)'
      read(*,*)nclass1
      write(*,*)'Density particle 2 (kg/m^3)'
      read(*,*)dens2
      write(*,*)'Particle 2 equivalent diameter of the median size (mm)'
      read(*,*)d2mm
      write(*,*)'Sorting particle grainsize distribution 2(phi)'
      read(*,*)sigma2
      write(*,*)'Particle 2 shape factor (-)'
      read(*,*)psi2
      write(*,*)'Classes number particle grainsize distribution 2 (-)'
      read(*,*)nclass2
      write(*,*)'Significance level t-test'
      read(*,*)probt
      write(*,*)'Layer thickness (m)'
      read(*,*)zlam
      write(*,*)'Sublayer thickness (m) (write 0 if it is not known)'
      read(*,*)zlams
      write(*,*)'Particle concentration in the layer (-)'
      read(*,*)c0
      write(*,*)'Substrate roughness (m)'
      read(*,*)ks
      endif
      case default
      write(*,*)'Wrong choice'
      goto 120
      end select
      alfa=probt/2.d0
!     Check: by default phase 2 is represented by less dense particles. If this is not the case on input, here the program exchanges phase 1 with phase 2
      if(dens2.gt.dens1.and.nmod.eq.2) then
      dentmp=dens1
      dtemp=d1mm
      sigtmp=sigma1
      psitmp=psi1
      ncltmp=nclass1
      dens1=dens2
      d1mm=d2mm
      sigma1=sigma2
      psi1=psi2
      nclass1=nclass2
      dens2=dentmp
      d2mm=dtemp
      sigma2=sigtmp
      psi2=psitmp
      nclass2=ncltmp
      endif
      if(zlams.eq.0.d0) zlams=d2mm/1000.d0                               ! If sublayer is not known (thus user set 0), it is set equal to the less dense particle diameter
      open(50,file='results.dat')
      open(52,file='log.dat')
      open(53,file='conc_profile.dat')
      open(54,file='pdyn_profile.dat')
      open(55,file='vel_profile.dat')
      open(56,file='dens_profile.dat')
      write(*,*)'Data summary'
      write(50,*)'Data summary'
      if(nmod.eq.1) then
      write(*,*)'Two layer method'                                      ! Data input summary
      write(*,190)dens,dm,dens2,d2mm,sigma2,psi2,nclass2,probt,zlam,zlams,c0,ks
      else
      write(*,*)'Two component method'
      write(*,191)dens1,d1mm,sigma1,psi1,nclass1,dens2,d2mm,sigma2,psi2,nclass2,probt,zlam,zlams,c0,ks
      endif
  130 write(*,*)'Correct data?'
      write(*,*)'1:Yes 2:No'
      read(*,*)ndata
      select case (ndata)
      case (1)
      if(nmod.eq.1) then
      write(50,*)'Two layer method'
      write(50,190)dens,dm,dens2,d2mm,sigma2,psi2,nclass2,probt,zlam,zlams,c0,ks
  190 format('Density of the entrained particle (kg/m^3)',7x,f10.3,/,&
     &'Diameter of the entrained particle (m)',11x,f10.3,/,&
     &'Particle density (kg/m^3)',24x,f10.3,/,&
     &'Particle equivalent diameter of the median size (mm)',f8.4,/,&
     &'Sorting particle grainsize distribution (phi)',6x,f8.3,/,&
     &'Particle shape factor (-)',26x,f8.3,/,&
     &'Classes number particle grainsize distribution (-)',2x,f4.0,/&
     &'Significance level t-test (-)',22x,f8.3,/,&
     &'Layer thickness (m)',34x,f8.5,/,&
     &'Sublayer thickness (m)',31x,f8.5,/,&
     &'Particle concentration in the layer (-)',12x,f8.3,/,&
     &'Substrate roughness (m)',30x,f8.5,/)
      else
      write(50,*)'Two component method'
      write(50,191)dens1,d1mm,sigma1,psi1,nclass1,dens2,d2mm,sigma2,psi2,&
     &nclass2,probt,zlam,zlams,c0,ks
  191 format('Density particle 1 (kg/m^3)',25x,f10.3,/,&
     &'Particle 1 equivalent diameter of the median size (mm)',f8.3,/,&
     &'Sorting particle grainsize distribution 1 (phi)',7x,f8.3,/,&
     &'Particle 1 shape factor (-)',27x,f8.3,/,&
     &'Classes number particle grainsize distribution 1 (-)',3x,f4.0,/,&
     &'Density particle 2 (kg/m^3)',25x,f10.3,/,&
     &'Particle 2 equivalent diameter of the median size (mm)',f8.3,/,&
     &'Sorting particle grainsize distribution 2 (phi)',7x,f8.3,/,&
     &'Particle 2 shape factor (-)',27x,f8.3,/,&
     &'Classes number particle grainsize distribution 2 (-)',3x,f4.0,/,&
     &'Significance level t-test (-)',25x,f8.3,/,&
     &'Layer thickness (m)',37x,f8.5,/,&
     &'Sublayer thickness (m)',34x,f8.5,/,&
     &'Particle concentration in the layer (-)',15x,f8.3,/,&
     &'Substrate roughness (m)',33x,f8.5,/)
      endif
      case (2)
      close (51)
      goto 120
      case default
      write(*,*)'Wrong choice'
      goto 130
      end select
! CALL TWOLAYER OR TWOCOMPONENT DEPENDING ON THE USER'S CHOICE #######################
      select case (nmod)
      case (1)
!     Two layer model
      call twolayer
      case (2)
!     Two componentS model
      call twocomponent
      end select
! CALL PROFILES FOR CALCULATING VERTICAL PROFILES OF FLUID-DYNAMIC VARIABLES ######################
      call profiles
      endif
! Write results
      write(*,*)''
      write(*,*)'Results'
      write(50,*)''
      write(50,*)'Results'
      write(50,200)dennrm,denmax,denmin,ztavg,ztmax,ztmin,zsfavg,zsfmax,&
     &zsfmin,ushavg,ushmax,ushmin,tauavg,taumax,taumin,p10avg,p10max,&
     &p10min,c2avg,c2max,c2min
      write(*,200)dennrm,denmax,denmin,ztavg,ztmax,ztmin,zsfavg,zsfmax,&
     &zsfmin,ushavg,ushmax,ushmin,tauavg,taumax,taumin,p10avg,p10max,&
     &p10min,c2avg,c2max,c2min
  200 format('Average density (kg/m^3)                            '&
     &,f10.3,/,&
     &'Maximum density (kg/m^3)                            ',f10.3,/,&
     &'Minimum density (kg/m^3)                            ',f10.3,/,&
     &'Average total flow thickness (m)                    ',f10.3,/,&
     &'Maximum total flow thickness (m)                    ',f10.3,/,&
     &'Minimum total flow thickness (m)                    ',f10.3,/,&
     &'Average shear flow thickness (m)                    ',f10.3,/,&
     &'Maximum shear flow thickness (m)                    ',f10.3,/,&
     &'Minimum shear flow thickness (m)                    ',f10.3,/,&
     &'Average shear velocity (m/s)                        ',f10.3,/,&
     &'Maximum shear velocity (m/s)                        ',f10.3,/,&
     &'Minimum shear velocity (m/s)                        ',f10.3,/,&
     &'Average velocity shear stress (Pa)                  ',f10.3,/,&
     &'Maximum velocity shear stress (Pa)                  ',f10.3,/,&
     &'Minimum velocity shear stress (Pa)                  ',f10.3,/,&
     &'Average specific 10m dynamic pressure (Pa)          ',f10.3,/,&
     &'Maximum specific 10m dynamic pressure (Pa)          ',f10.3,/,&
     &'Minimum specific 10m dynamic pressure (Pa)          ',f10.3,/,&
     &'Average 2m particle concentration (-)               ',e10.5,/,&
     &'Maximum 2m particle concentration (-)               ',e10.5,/,&
     &'Minimum 2m particle concentration (-)               ',e10.5,/)
      write(*,*)'### User requested outputs ###'
      write(50,*)'### User requested outputs ###'
      write(52,*)'### User requested outputs ###'
      write(52,201)p10av1,p10mx1,p10mn1,c2av1,c2max1,c2min1
  201 format('50th percentile specific 10m dynamic pressure (Pa)  '&
     &,f10.3,/,&
     &'84th percentile specific 10m dynamic pressure (Pa)  ',f10.3,/,&
     &'16th percentile specific 10m dynamic pressure (Pa)  ',f10.3,//,&
     &'50th percentile 2m particle concentration (-)       ',e10.5,/,&
     &'84th percentile 2m particle concentration (-)       ',e10.5,/,&
     &'16th percentile 2m particle concentration (-)       ',e10.5,/)
      if(ndynpr.eq.2.and.nc.eq.2) write(*,*)'No user requested outputs'
      if(ndynpr.eq.2) goto 210
      do i=1,ipr
      write(*,202)zdynpr(i),pzavg(i),pzmax(i),pzmin(i)
      write(50,202)zdynpr(i),pzavg(i),pzmax(i),pzmin(i)
      write(*,203)zdynpr(i),pzav1(i),pzmax1(i),pzmin1(i)
      write(52,203)zdynpr(i),pzav1(i),pzmax1(i),pzmin1(i)
      enddo
  210 if(nc.eq.2) goto 211
      do i=1,ic
      write(*,204)zc(i),czavg(i),czmax(i),czmin(i)
      write(50,204)zc(i),czavg(i),czmax(i),czmin(i)
      write(52,205)zc(i),czav1(i),czmax1(i),czmin1(i)
      write(*,205)zc(i),czav1(i),czmax1(i),czmin1(i)
      enddo
  202 format('z =',f6.2,/,&
     &'Average specific dynamic pressure (Pa)            ',f10.3,/,&
     &'Maximum specific dynamic pressure (Pa)            ',f10.3,/,&
     &'Minimum specific dynamic pressure (Pa)            ',f10.3,//)
  203 format('z =',f6.2,/,&
     &'50th percentile specific dynamic pressure (Pa)    ',f10.3,/,&
     &'84th percentile specific dynamic pressure (Pa)    ',f10.3,/,&
     &'16th percentile specific dynamic pressure (Pa)    ',f10.3,//)
  204 format('z =',f6.2,/,&
     &'Average particle concentration            ',e10.5,/,&
     &'Maximum particle concentration            ',e10.5,/,&
     &'Minimum particle concentration            ',e10.5,//)
  205 format('z =',f6.2,/,&
     &'50th percentile particle concentration    ',e10.5,/,&
     &'84th percentile particle concentration    ',e10.5,/,&
     &'16th percentile particle concentration    ',e10.5,//)
  211 write(*,*)''
      write(*,*)'Test t-Student summary'
      write(50,*)''
      write(50,*)'Test t-Student summary'
      write(50,206)alfa,ttab,tcalc
      write(*,206)alfa,ttab,tcalc
  206 format('Significance level   ',f10.6,/,&
     &'Theoretical t value ',f8.3,/,&
     &'Calculated t value  ',f8.3,//)
      write(*,*)'###### PROBABILITY FUNCTIONS ######'
      write(50,*)'###### PROBABILITY FUNCTIONS ######'
! PROBABILITY FUNCTIONS CALCULATION FOR DYNAMIC PRESSURE AND PARTICLE CONCENTRATION##################################
      call probfunction
      end
