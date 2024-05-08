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
	else
		write(fout,200)dennrm,denmin,denmax,z0avg,z0min,z0max,ztavg,ztmax,ztmin,zsfavg,zsfmax,&
	   &zsfmin,ushavg,ushmax,ushmin,tauavg,taumax,taumin,pnsavg,pnsmin,pnsmax,p10avg,p10max,&
	   &p10min,c2avg,c2min,c2max,c2dpavg,c2dpmin,c2dpmax
		write(*,200)dennrm,denmin,denmax,z0avg,z0min,z0max,ztavg,ztmax,ztmin,zsfavg,zsfmax,&
	   &zsfmin,ushavg,ushmax,ushmin,tauavg,taumax,taumin,pnsavg,pnsmin,pnsmax,p10avg,p10max,&
	   &p10min,c2avg,c2min,c2max,c2dpavg,c2dpmin,c2dpmax
		if(tavg.ne.undefined.and.tmax.ne.undefined.and.tmin.ne.undefined) then
			write(fout, 212) rhogavg, rhogmax, rhogmin, tavg, tmin, tmax
			write(*, 212) rhogavg, rhogmax, rhogmin, tavg, tmin, tmax
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
		if(usr_z_dynpr.eqv..FALSE..and.usr_z_c.eqv..FALSE..and.usr_z_t.eqv..FALSE.) write(*,*)'No user requested outputs'
		if(usr_z_dynpr) then
			do i=1,20
				if(zdynpr(i).eq.undefined) cycle
				write(*,202)zdynpr(i),pzavg(i),pzmax(i),pzmin(i)
				write(fout,202)zdynpr(i),pzavg(i),pzmax(i),pzmin(i)
			enddo
		endif
		if(usr_z_c) then
			do i=1,20
				if(zc(i).eq.undefined) cycle
				write(*,204)zc(i),czavg(i),czmin(i),czmax(i)
				write(fout,204)zc(i),czavg(i),czmin(i),czmax(i)
			enddo
			do i=1,20
				if(zc(i).eq.undefined) cycle
				write(*,205)zc(i),czdpavg(i),czdpmin(i),czdpmax(i)
				write(fout,205)zc(i),czdpavg(i),czdpmin(i),czdpmax(i)
			enddo
		endif
		if(usr_z_t.and.calc_t_mix) then
			do i=1,20
				if(zt(i).eq.undefined) cycle
				write(*,213)zt(i),tzavg(i),tzmax(i),tzmin(i)
				write(fout,213)zt(i),tzavg(i),tzmax(i),tzmin(i)
			enddo
		endif
	  
		if(deprates) then
			write(*,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
			write(fout,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
			write(flog,*)'###DEPOSITION RATE AND TIME CALCULATIONS'
			do i=1,kmax
				write(flog,*)i
				if(i.eq.1) then
					write(*,*)'50th percentile solution'
					write(fout,*)'50th percentile solution'
					write(flog,*)'50th percentile solution'
				else if(i.eq.2) then
					write(*,*)'84th percentile solution'
					write(fout,*)'84th percentile solution'
					write(flog,*)'84th percentile solution'
				else
					write(*,*)'16th percentile solution'
					write(fout,*)'16th percentile solution'
					write(flog,*)'16th percentile solution'
				endif
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
	endif

    write(*,*)'###### PROBABILITY FUNCTIONS ######'
    write(fout,*)'###### PROBABILITY FUNCTIONS ######'
    close(fout)
      
200 format('50th percentile solution flow density (kg/m^3)                         ',f10.3,/,&
     &'84th percentile solution flow density (kg/m^3)                         ',f10.3,/,&
     &'16th percentile solution flow density (kg/m^3)                         ',f10.3,/,&
     &'50th percentile solution viscous sublayer thickness (m)                ',f10.5,/,&
     &'84th percentile solution viscous sublayer thickness (m)                ',f10.5,/,&
     &'16th percentile solution viscous sublayer thickness (m)                ',f10.5,/,&	 
     &'50th percentile solution total flow thickness (m)                      ',f10.3,/,&
     &'84th percentile solution total flow thickness (m)                      ',f10.3,/,&
     &'16th percentile solution total flow thickness (m)                      ',f10.3,/,&
     &'50th percentile solution shear flow thickness (m)                      ',f10.3,/,&
     &'84th percentile solution shear flow thickness (m)                      ',f10.3,/,&
     &'16th percentile solution shear flow thickness (m)                      ',f10.3,/,&
     &'50th percentile solution shear velocity (m/s)                          ',f10.3,/,&
     &'84th percentile solution shear velocity (m/s)                          ',f10.3,/,&
     &'16th percentile solution shear velocity (m/s)                          ',f10.3,/,&
     &'50th percentile solution velocity shear stress (Pa)                    ',f10.3,/,&
     &'84th percentile solution velocity shear stress (Pa)                    ',f10.3,/,&
     &'16th percentile solution velocity shear stress (Pa)                    ',f10.3,/,&
     &'50th percentile solution suspension Rouse number (-)                   ',f10.3,/,&
     &'84th percentile solution suspension Rouse number (-)                   ',f10.3,/,&
     &'16th percentile solution suspension Rouse number (-)                   ',f10.3,/,&
     &'50th percentile solution 10m depth-averaged dynamic pressure (Pa)      ',f10.3,/,&
     &'84th percentile solution 10m depth-averaged dynamic pressure (Pa)      ',f10.3,/,&
     &'16th percentile solution 10m depth-averaged dynamic pressure (Pa)      ',f10.3,/,&
     &'50th percentile solution 2m particle concentration (-)                 ',e10.5,/,&
     &'84th percentile solution 2m particle concentration (-)                 ',e10.5,/,&
     &'16th percentile solution 2m particle concentration (-)                 ',e10.5,/,&
     &'50th percentile solution 2m depth-averaged particle concentration (-)  ',e10.5,/,&
     &'84th percentile solution 2m depth-averaged particle concentration (-)  ',e10.5,/,&
     &'16th percentile solution 2m depth-averaged particle concentration (-)  ',e10.5)!&
  201 format('50th percentile specific 10m dynamic pressure (Pa)       ',f10.3,/,&
     &'84th percentile specific 10m dynamic pressure (Pa)  ',f10.3,/,&
     &'16th percentile specific 10m dynamic pressure (Pa)  ',f10.3,//,&
     &'50th percentile 2m particle concentration (-)       ',e10.5,/,&
     &'84th percentile 2m particle concentration (-)       ',e10.5,/,&
     &'16th percentile 2m particle concentration (-)       ',e10.5,/)
  202 format('z =',f6.2,/,&
     &'50th percentile solution depth-averaged dynamic pressure (Pa)            ',f10.3,/,&
     &'84th percentile solution depth-averaged dynamic pressure (Pa)            ',f10.3,/,&
     &'16th percentile solution depth-averaged dynamic pressure (Pa)            ',f10.3,//)
  204 format('z =',f6.2,/,&
     &'50th percentile solution particle concentration            ',e10.5,/,&
     &'84th percentile solution particle concentration            ',e10.5,/,&
     &'16th percentile solution particle concentration            ',e10.5,//)
  205 format('z =',f6.2,/,&
     &'50th percentile solution depth-averaged particle concentration            ',e10.5,/,&
     &'84th percentile solution depth-averaged particle concentration            ',e10.5,/,&
     &'16th percentile solution depth-averaged particle concentration            ',e10.5,//)
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
  212 format('50th percentile solution gas density (kg/m^3)                        '&
     &,f10.3,/,&
     &'84th percentile solution gas density (kg/m^3)                        ',f10.3,/,&
     &'16th percentile solution gas density (kg/m^3)                        ',f10.3,/,&
     &'50th percentile solution flow temperature (K)                        ',f10.3,/,&
     &'84th percentile solution flow temperature (K)                        ',f10.3,/,&	 
	 &'16th percentile solution flow temperature (K)                        ',f10.3,//)
  213 format('z =',f6.2,/,&
     &'50th percentile solution flow temperature (K)            ',f10.3,/,&
     &'84th percentile solution flow temperature (K)            ',f10.3,/,&
     &'16th percentile solution flow temperature (K)            ',f10.3,//)
  214 format('z =',f6.2,/,&
     &'50th percentile flow temperature (K)            ',f10.3,/,&
     &'84th percentile flow temperature (K)            ',f10.3,/,&
     &'16th percentile flow temperature (K)            ',f10.3,//)
      end subroutine write_results
