      subroutine write_results
      USE inoutdata; USE nrtype
      implicit none
      integer :: i
	  open(fout, file=trim(output_dir)//trim(path_sep)//'results.dat', position="append", status="old", action="write")
      write(*,*)''
      write(*,*)'Results'
      write(fout,*)''
      write(fout,*)'Results'
      if(only_deprates) then
	      write(*,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
		  write(fout,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
		  write(flog,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
		  do i=1,kmax
			  write(*,207)i
			  write(fout,207)i
			  write(fout,208)zlam_final(i),rtot_susp(i),tdep_susp(i),ctot_flow(i),ctot_dep(i),ctot_susp(i)
			  write(*,208)zlam_final(i),rtot_susp(i),tdep_susp(i),ctot_flow(i),ctot_dep(i),ctot_susp(i)
			  write(fout,211)srw(i),qtot(i),sqratio(i)
			  write(*,211)srw(i),qtot(i),sqratio(i)
			  if(zlam_massive.eq.undefined) cycle
			  write(fout,209)rtot_massive(i),tdep_massive(i),ctot_massive(i)
			  write(*,209)rtot_massive(i),tdep_massive(i),ctot_massive(i)
			  write(*,210)rtot_susp(i)+rtot_massive(i),tdep_susp(i)+tdep_massive(i)
			  write(fout,210)rtot_susp(i)+rtot_massive(i),tdep_susp(i)+tdep_massive(i)
		  enddo
	  endif
      write(fout,200)dennrm,denmax,denmin,z0avg,z0max,z0min,ztavg,ztmax,ztmin,zsfavg,zsfmax,&
     &zsfmin,ushavg,ushmax,ushmin,tauavg,taumax,taumin,pnsavg,pnsmax,pnsmin,p10avg,p10max,&
     &p10min,c2avg,c2max,c2min
      write(*,200)dennrm,denmax,denmin,z0avg,z0max,z0min,ztavg,ztavg,ztmax,ztmin,zsfavg,zsfmax,&
     &zsfmin,ushavg,ushmax,ushmin,tauavg,taumax,taumin,pnsavg,pnsmax,pnsmin,p10avg,p10max,&
     &p10min,c2avg,c2max,c2min
	  if(rhogavg.ne.undefined.and.rhogmax.ne.undefined.and.rhogmin.ne.undefined) then
		write(fout, 212) rhogavg, rhogmax, rhogmin, tavg, tmax, tmin
		write(*, 212) rhogavg, rhogmax, rhogmin, tavg, tmax, tmin
	  endif
      write(*,*)''
      write(*,*)'Test t-Student summary'
      write(fout,*)''
      write(fout,*)'Test t-Student summary'
      write(fout,206)alfa,ttab,tcalc
      write(*,206)alfa,ttab,tcalc

      write(*,*)'### User requested outputs ###'
      write(fout,*)'### User requested outputs ###'
      write(flog,*)'### User requested outputs ###'
      write(flog,201)p10av1,p10mx1,p10mn1,c2av1,c2max1,c2min1
	  write(fout,201)p10av1,p10mx1,p10mn1,c2av1,c2max1,c2min1
      if(usr_z_dynpr.eqv..FALSE..and.usr_z_c.eqv..FALSE..and.usr_z_t.eqv..FALSE.) write(*,*)'No user requested outputs'
      if(usr_z_dynpr) then
		  do i=1,ipr
			  write(*,202)zdynpr(i),pzavg(i),pzmax(i),pzmin(i)
			  write(fout,202)zdynpr(i),pzavg(i),pzmax(i),pzmin(i)
			  write(*,203)zdynpr(i),pzav1(i),pzmax1(i),pzmin1(i)
			  write(flog,203)zdynpr(i),pzav1(i),pzmax1(i),pzmin1(i)
		  enddo
	  endif
      if(usr_z_c) then
		  do i=1,ic
			  write(*,204)zc(i),czavg(i),czmax(i),czmin(i)
			  write(fout,204)zc(i),czavg(i),czmax(i),czmin(i)
			  write(flog,205)zc(i),czav1(i),czmax1(i),czmin1(i)
			  write(*,205)zc(i),czav1(i),czmax1(i),czmin1(i)
		  enddo
	  endif
      if(usr_z_t.or.calc_t_mix) then
		  do i=1,itemp
			  write(*,213)zt(i),tzav(i),tzmax(i),tzmin(i)
			  write(fout,213)zt(i),tzav(i),tzmax(i),tzmin(i)
			  write(flog,214)zt(i),tzav1(i),tzmax1(i),tzmin1(i)
			  write(*,214)zt(i),tzav(i),tzmax(i),tzmin(i)
		  enddo
	  endif
  
      if(deprates) then
	      write(*,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
		  write(fout,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
		  write(flog,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
		  do i=1,kmax
			  write(*,207)i
			  write(fout,207)i
			  write(fout,208)zlam_final(i),rtot_susp(i),tdep_susp(i),ctot_flow(i),ctot_dep(i),ctot_susp(i)
			  write(*,208)zlam_final(i),rtot_susp(i),tdep_susp(i),ctot_flow(i),ctot_dep(i),ctot_susp(i)
			  write(fout,211)srw(i),qtot(i),sqratio(i)
			  write(*,211)srw(i),qtot(i),sqratio(i)
			  if(zlam_massive.eq.undefined) cycle
			  write(fout,209)rtot_massive(i),tdep_massive(i),ctot_massive(i)
			  write(*,209)rtot_massive(i),tdep_massive(i),ctot_massive(i)
			  write(*,210)rtot_susp(i)+rtot_massive(i),tdep_susp(i)+tdep_massive(i)
			  write(fout,210)rtot_susp(i)+rtot_massive(i),tdep_susp(i)+tdep_massive(i)
		  enddo
	  endif

      write(*,*)'###### PROBABILITY FUNCTIONS ######'
      write(fout,*)'###### PROBABILITY FUNCTIONS ######'
      close(fout)
      
  200 format('Average density (kg/m^3)                            '&
     &,f10.3,/,&
     &'Maximum density (kg/m^3)                            ',f10.3,/,&
     &'Minimum density (kg/m^3)                            ',f10.3,/,&
     &'Average viscous sublayer thickness (m)              ',f10.5,/,&
     &'Maximum viscous sublayer thickness (m)              ',f10.5,/,&
     &'Minimum viscous sublayer thickness (m)              ',f10.5,/,&	 
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
     &'Minimum 2m particle concentration (-)               ',e10.5)!&
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
      'Total particle concentration of turbulent suspension = ',e10.3,/, &
      'Total particle concentration of turbulent suspension at deposition = ',e10.3,/, &
      'Total particle concentration of turbulent suspension remaining at suspension = ',e10.3,/)
  209 format('***WASH LOAD CONTRIBUTING TO THE FINE MASSIVE LAYER***',/, &
      'Deposition rate of wash load (kg m^-2 s^-1) = ',e10.3,/, &
      'Deposition time of wash load(s) = ',e10.3,/, &
      'Total particle concentration of wash load= ',e10.3,//)
  210 format('***TOTAL DEPOSITION RATE AND TIME***',/, &
      'Total deposition rate (kg m^-2 s^-1) = ',e10.3,/, &
      'Total deposition time (s) = ',e10.3,//)
  211 format('Volumetric sedimentation rate per unit width (m^2 s^-1) = ',e10.3,/, &
      'Bedload transportation rate per unit width (m^2 s^-1) = ',e10.3,/, &
      'Srw/Qb ratio (-) = ',f8.3,//)
  212 format('Average gas density (kg/m^3)                        '&
     &,f10.3,/,&
     &'Maximum gas density (kg/m^3)                        ',f10.3,/,&
     &'Minimum gas density (kg/m^3)                        ',f10.3,/,&
     &'Average flow temperature (K)                        ',f10.3,/,&
     &'Maximum flow temperature (K)                        ',f10.3,/,&	 
	 &'Minimum flow temperature (K)                        ',f10.3,//)
  213 format('z =',f6.2,/,&
     &'Average flow temperature (K)            ',f10.3,/,&
     &'Maximum flow temperature (K)            ',f10.3,/,&
     &'Minimum flow temperature (K)            ',f10.3,//)
  214 format('z =',f6.2,/,&
     &'50th percentile flow temperature (K)            ',f10.3,/,&
     &'84th percentile flow temperature (K)            ',f10.3,/,&
     &'16th percentile flow temperature (K)            ',f10.3,//)
      end subroutine write_results
