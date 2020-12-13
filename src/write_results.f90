      subroutine write_results
      USE inoutdata; USE nrtype
      implicit none
      integer :: i
      write(*,*)''
      write(*,*)'Results'
      write(50,*)''
      write(50,*)'Results'
      if(only_deprates) goto 212
      write(50,200)dennrm,denmax,denmin,ztavg,ztmax,ztmin,zsfavg,zsfmax,&
     &zsfmin,ushavg,ushmax,ushmin,tauavg,taumax,taumin,pnsavg,pnsmax,pnsmin,p10avg,p10max,&
     &p10min,c2avg,c2max,c2min!,rtot(1),rtot(2),rtot(3),tdep(1),tdep(3),tdep(2)
      write(*,200)dennrm,denmax,denmin,ztavg,ztmax,ztmin,zsfavg,zsfmax,&
     &zsfmin,ushavg,ushmax,ushmin,tauavg,taumax,taumin,pnsavg,pnsmax,pnsmin,p10avg,p10max,&
     &p10min,c2avg,c2max,c2min!,rtot(1),rtot(2),rtot(3),tdep(1),tdep(3),tdep(2)

      write(*,*)''
      write(*,*)'Test t-Student summary'
      write(50,*)''
      write(50,*)'Test t-Student summary'
      write(50,206)alfa,ttab,tcalc
      write(*,206)alfa,ttab,tcalc

      write(*,*)'### User requested outputs ###'
      write(50,*)'### User requested outputs ###'
      write(52,*)'### User requested outputs ###'
      write(52,201)p10av1,p10mx1,p10mn1,c2av1,c2max1,c2min1
      if(usr_z_dynpr.eqv..FALSE..and.usr_z_c.eqv..FALSE.) write(*,*)'No user requested outputs'
      if(usr_z_dynpr.eqv..FALSE.) goto 210
      do i=1,ipr
      write(*,202)zdynpr(i),pzavg(i),pzmax(i),pzmin(i)
      write(50,202)zdynpr(i),pzavg(i),pzmax(i),pzmin(i)
      write(*,203)zdynpr(i),pzav1(i),pzmax1(i),pzmin1(i)
      write(52,203)zdynpr(i),pzav1(i),pzmax1(i),pzmin1(i)
      enddo
  210 if(usr_z_c.eqv..FALSE.) goto 211
      do i=1,ic
      write(*,204)zc(i),czavg(i),czmax(i),czmin(i)
      write(50,204)zc(i),czavg(i),czmax(i),czmin(i)
      write(52,205)zc(i),czav1(i),czmax1(i),czmin1(i)
      write(*,205)zc(i),czav1(i),czmax1(i),czmin1(i)
      enddo

  211 if(.not.deprates) goto 213

  212 write(*,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
      write(50,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
      write(52,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
      do i=1,kmax
      write(*,207)i
      write(50,207)i
      write(50,208)zlam_final(i),rtot_susp(i),tdep_susp(i),rho_flow(i),ush_flow(i),ctot_flow(i),ctot_dep(i),ctot_susp(i)
      write(*,208)zlam_final(i),rtot_susp(i),tdep_susp(i),rho_flow(i),ush_flow(i),ctot_flow(i),ctot_dep(i),ctot_susp(i)
      write(50,215)srw(i),qtot(i),sqratio(i)
      write(*,215)srw(i),qtot(i),sqratio(i)
      if(zlam_massive.eq.undefined) cycle
      write(50,209)rtot_massive(i),tdep_massive(i),ctot_massive(i)
      write(*,209)rtot_massive(i),tdep_massive(i),ctot_massive(i)
      write(*,214)rtot_susp(i)+rtot_massive(i),tdep_susp(i)+tdep_massive(i)
      write(50,214)rtot_susp(i)+rtot_massive(i),tdep_susp(i)+tdep_massive(i)
      enddo

  213 write(*,*)'###### PROBABILITY FUNCTIONS ######'
      write(50,*)'###### PROBABILITY FUNCTIONS ######'
      
      
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
     &'Average suspension Rouse number (-)                 ',f10.3,/,&
     &'Maximum suspension Rouse number (-)                 ',f10.3,/,&
     &'Minimum suspension Rouse number (-)                 ',f10.3,/,&
     &'Average specific 10m dynamic pressure (Pa)          ',f10.3,/,&
     &'Maximum specific 10m dynamic pressure (Pa)          ',f10.3,/,&
     &'Minimum specific 10m dynamic pressure (Pa)          ',f10.3,/,&
     &'Average 2m particle concentration (-)               ',e10.5,/,&
     &'Maximum 2m particle concentration (-)               ',e10.5,/,&
     &'Minimum 2m particle concentration (-)               ',e10.5,/,/)!&
!     &'Average deposition rate (kg/(m^2*s))                ',f10.3,/,&
!     &'Maximum deposition rate (kg/(m^2*s))                ',f10.3,/,&
!     &'Minimum deposition rate (kg/(m^2*s))                ',f10.3,/,&
!     &'Average deposition time (d)                         ',f10.3,/,&
!     &'Maximum deposition time (d)                         ',f10.3,/,&
!     &'Minimum deposition time (d)                         ',f10.3,/)
  201 format('50th percentile specific 10m dynamic pressure (Pa)  '&
     &,f10.3,/,&
     &'84th percentile specific 10m dynamic pressure (Pa)  ',f10.3,/,&
     &'16th percentile specific 10m dynamic pressure (Pa)  ',f10.3,//,&
     &'50th percentile 2m particle concentration (-)       ',e10.5,/,&
     &'84th percentile 2m particle concentration (-)       ',e10.5,/,&
     &'16th percentile 2m particle concentration (-)       ',e10.5,/)
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
  206 format('Significance level   ',f10.6,/,&
     &'Theoretical t value ',f8.3,/,&
     &'Calculated t value  ',f8.3,//)
  207 format(/,'SOLUTION ',i1,/)
  208 format(/,'******SUSPENSION LOAD******',/, &
      'Laminated layer thickness (m) = ',e10.3,/, &
      'Deposition rate of turbulent suspension (kg m^-2 s^-1) = ',e10.3,/, &
      'Deposition time of turbulent suspension (s) = ',e10.3,/, &
      'Flow density (kg m-3) = ',e10.3,/, &
      'Flow shear velocity (m s-1) = ',e10.3,/, &
      'Total particle concentration of turbulent suspension = ',e10.3,/, &
      'Total particle concentration of turbulent suspension at deposition = ',e10.3,/, &
      'Total particle concentration of turbulent suspension remaining at suspension = ',e10.3,/)
  209 format('***WASH LOAD CONTRIBUTING TO THE FINE MASSIVE LAYER***',/, &
      'Deposition rate of wash load (kg m^-2 s^-1) = ',e10.3,/, &
      'Deposition time of wash load(s) = ',e10.3,/, &
      'Total particle concentration of wash load= ',e10.3,//)
  214 format('***TOTAL DEPOSITION RATE AND TIME***',/, &
      'Total deposition rate (kg m^-2 s^-1) = ',e10.3,/, &
      'Total deposition time (s) = ',e10.3,//)
  215 format('Volumetric sedimentation rate per unit width (m^2 s^-1) = ',e10.3,/, &
      'Bedload transportation rate per unit width (m^2 s^-1) = ',e10.3,/, &
      'Srw/Qb ratio (-) = ',f8.3,//)
      end subroutine write_results
