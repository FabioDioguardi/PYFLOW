   subroutine profiles(convergence)
   USE inoutdata; USE nrtype
   implicit none
    
   INTERFACE
	  FUNCTION func(x)
		 USE inoutdata; USE nrtype
		 IMPLICIT NONE
		 REAL(dp) :: func
		 REAL(dp), INTENT(IN) :: x
	  END FUNCTION func

	  FUNCTION grad(xx)
		 USE inoutdata; USE nrtype
		 implicit none
		 REAL(dp) :: grad
		 real(dp), INTENT(IN) :: xx
	  END FUNCTION grad

	  FUNCTION rad(xx)
		 USE inoutdata; USE nrtype
		 implicit none
		 REAL(dp) :: rad
		 real(dp), INTENT(IN) :: xx
	  END FUNCTION rad

	  FUNCTION qsimp(a, b)
		 USE nrtype; USE nrutil, ONLY: nrerror
		 REAL(dp), INTENT(IN) :: a, b
		 REAL(dp) :: qsimp
	  END FUNCTION qsimp
	
   END INTERFACE
   
   real(dp), dimension(3) :: x
   real(dp), dimension(5) :: x_gas
   real(dp), dimension(3) :: arr
   real(dp) :: z, u
   real(dp) :: cavg, cmax, cmin, denavg, dfmax, dfmin, pdynav, pdynmx, pdynmn,&
   s, senx, uavg, umax, umin, dz0, avg_z0, avg_pns, avg_rhog, avg_zt, avg_zsf
   integer :: i, j, i_x, n_solutions_z0
   LOGICAL :: check, convergence
   character(42) :: avgcon, maxcon, mincon
   character(37) :: avgpd, maxpd, minpd
   character(30) :: avgvel, maxvel, minvel
   character(33) :: avgden, maxden, minden
   character(32) :: avgtemp, maxtemp, mintemp
   character(10) :: height, cmd
   avgcon = '50th percentile particle concentration (-)'
   maxcon = '84th percentile particle concentration (-)'
   mincon = '16th percentile particle concentration (-)'
   avgpd = '50th percentile dynamic pressure (Pa)'
   maxpd = '84th percentile dynamic pressure (Pa)'
   minpd = '16th percentile dynamic pressure (Pa)'
   avgvel = '50th percentile velocity (m/s)'
   maxvel = '84th percentile velocity (m/s)'
   minvel = '16th percentile velocity (m/s)'
   avgden = '50th percentile density (kg/m^3)'
   maxden = '84th percentile density (kg/m^3)'
   minden = '16th percentile density (kg/m^3)'
   avgtemp = '50th percentile temperature (K)'
   maxtemp = '84th percentile temperature (K)'
   mintemp = '16th percentile temperature (K)'
   height = 'Height (m)'
   calc_t_mix = .false.
   convergence = .true.
   rho_air = p_air / (r_air * t_air) !For calculating Tmix
   rho_gas = p_air / (r_gas * t_gas) !For calculating Tmix
   epsdz0 = 1.d-1
   !     Average, maximum and minimum shear flow particle concentration
   if (model .eq. 'TWOLAYERS') then
	  densp = rhos(1, 0)
   else
	  densp = rhos(2, 0)
   end if
   zlams_or = zlams
   if (zlams .eq. UNDEFINED) then
	  !Set zlam = 1 cm, the typical laminae thickness
	  write (*, *) 'WARNING! Command ZLAMS missing in input.dat'
	  write (*, *) 'Setting ZLAMS=1 cm'
	  write (flog, *) 'WARNING! Command ZLAMS missing in input.dat'
	  write (flog, *) 'Setting ZLAMS=1 cm'
	  z0temp = 0.01d0
   else
	  z0temp = zlams
   end if
   if (dengas .eq. undefined) then
	  calc_t_mix = .true.
	  ! 50th percentile solution
	  z0 = z0temp
	  dz0 = epsdz0 * z0
	  den = dennrm
	  zsfavg = tauavg/((den - rho_air)*g*sin(rad(slope_ground)))
	  nnewt = 3
	  zshr = zsfavg
	  write (flog, *) 'Pns avg, rho_g avg and ztot avg calculation residuals'
	  write (*, *) 'Pns avg, rho_g avg and ztot avg calculation residuals'
	  n_solutions_z0 = 0
	  avg_z0 = 0.d0
	  avg_pns = 0.d0
	  avg_rhog = 0.d0
	  avg_zt = 0.d0
	  do while(z0.gt.1.d-10) 
		 write(flog,*)'Solution with z0 = ', z0
		 write(*,*)'Solution with z0 = ', z0
		 x(1) = pnsavgguess
		 x(2) = rhogavgguess
		 x(3) = ztavgguess
		 call newt(x, check, 3)
		 ! Double check the solutions obtained from newt
		 do i_x = 1, size(x)
			if (x(i_x).lt.0.d0.or.x(2).gt.rho_air) check = .true.
			! Spurious convergence case, i.e. check=false but non realistic solutions obtained
			if(check) exit
		 enddo	
		 if(check) then
			! No convergence case
			if(zlams.eq.UNDEFINED) then
			   z0 = z0 - dz0
			   write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
			   check = .false.
			else
			   z0 = 0.01d0
			   dz0 = epsdz0 * z0
			   write(flog,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   zlams = UNDEFINED
			   cycle
			endif
		 else
			avg_z0 = avg_z0 + z0
			avg_pns = avg_pns + x(1)
			avg_rhog = avg_rhog + x(2)
			avg_zt = avg_zt + x(3)
			n_solutions_z0 = n_solutions_z0 + 1
			if(zlams_or.eq.UNDEFINED) then
				z0 = z0 - dz0
				cycle
			else
				if(zlams.eq.UNDEFINED) then
					z0 = z0 - dz0
					cycle
				else
					exit
				endif
			endif
		 endif
	  enddo
	  if(check.or.n_solutions_z0.eq.0) then
		 write(flog,*)'Unable to find a solution for Pns avg, rho_g avg and ztot avg'
		 write(*,*)'Unable to find a solution for Pns avg, rho_g avg and ztot avg'
		 convergence = .false.
		 if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
		 return
	  endif
	  z0avg = avg_z0 / n_solutions_z0
	  pnsavg = avg_pns / n_solutions_z0
	  rhogavg = avg_rhog / n_solutions_z0
	  ztavg = avg_zt / n_solutions_z0
	  write(*,360) pnsavg, rhogavg, ztavg, z0avg
	  write(flog,360) pnsavg, rhogavg, ztavg, z0avg
360   format('Pnsusp avg =', f6.3, 1x, 'rhogavg =', f9.6, 1x, 'ztavg =', f8.3, 1x, 'z0avg =', f9.6,/)
	  cavg = (dennrm - rhogavg)/(densp - rhogavg)
	  c_gas_avg = (rhogavg - rho_air) / (rho_gas - rho_air)
	  c_air_avg = 1.d0 - c_gas_avg
	  r_mix = c_gas_avg * r_gas + c_air_avg * r_air
	  t_mix_avg = p_air / (r_mix * rhogavg)
	  cgastemp = c_gas_avg
	  cairtemp = c_air_avg
	  nfunc = 21
	  tavg = func(cavg)
	  n_solutions_z0 = 0
	  avg_z0 = 0.d0
	  avg_pns = 0.d0
	  avg_rhog = 0.d0
	  avg_zt = 0.d0
	  ! 84th percentile solution
	  zlams = zlams_or
	  z0 = z0temp
	  dz0 = epsdz0 * z0
	  den = denmin
	  zsfmax = taumax/((den - rho_air)*g*sin(rad(slope_ground)))
	  zshr = zsfmax
	  write (flog, *) 'Pns min, rho_g max and ztot max calculation residuals'
	  write (*, *) 'Pns min, rho_g max and ztot max calculation residuals'
	  do while(z0.gt.1.d-10) 
		 write(flog,*)'Solution with z0 = ', z0
		 write(*,*)'Solution with z0 = ', z0
		 x(1) = pnsminguess
		 x(2) = rhogmaxguess
		 x(3) = ztmaxguess
		 call newt(x, check, 3)
		 ! Double check the solutions obtained from newt
		 do i_x = 1, size(x)
			if (x(i_x).lt.0.d0.or.x(2).gt.rho_air) check = .true.
			! Spurious convergence case, i.e. check=false but non realistic solutions obtained
			if(check) exit
		 enddo	
		 if(check) then
			! No convergence case
			if(zlams.eq.UNDEFINED) then
			   z0 = z0 - dz0
			   write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
			   check = .false.
			else
			   z0 = 0.01d0
			   dz0 = epsdz0 * z0
			   write(flog,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   zlams = UNDEFINED
			   cycle
			endif
		 else
			avg_z0 = avg_z0 + z0
			avg_pns = avg_pns + x(1)
			avg_rhog = avg_rhog + x(2)
			avg_zt = avg_zt + x(3)
			n_solutions_z0 = n_solutions_z0 + 1
			if(zlams_or.eq.UNDEFINED) then
				z0 = z0 - dz0
				cycle
			else
				if(zlams.eq.UNDEFINED) then
					z0 = z0 - dz0
					cycle
				else
					exit
				endif
			endif
		 endif
	  enddo
	  if(check.or.n_solutions_z0.eq.0) then
		 write(flog,*)'Unable to find a solution for Pns min, rho_g max and ztot max'
		 write(*,*)'Unable to find a solution for Pns min, rho_g max and ztot max'
		 convergence = .false.
		 if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
		 return
	  endif
	  z0min = avg_z0 / n_solutions_z0
	  pnsmin = avg_pns / n_solutions_z0
	  rhogmax = avg_rhog / n_solutions_z0
	  ztmax = avg_zt / n_solutions_z0
	  write(*,361) pnsmin, rhogmax, ztmax, z0min
	  write(flog,361) pnsmin, rhogmax, ztmax, z0min
361   format('Pnsusp min =', f6.3, 1x, 'rhogmax =', f9.6, 1x, 'ztmax =', f8.3, 1x, 'z0min =', f9.6,/)		
	  cmin = (denmin - rhogmax)/(densp - rhogmax) ! FABIO: ho cambiato in cmin per essere consistente
	  c_gas_max = (rhogmax - rho_air) / (rho_gas - rho_air)
	  c_air_max = 1.d0 - c_gas_max
	  r_mix = c_gas_max * r_gas + c_air_max * r_air
	  t_mix_min = p_air / (r_mix * rhogmax) 
	  cgastemp = c_gas_max
	  cairtemp = c_air_max
	  nfunc = 21
	  tmin = func(cmin) 
	  n_solutions_z0 = 0
	  avg_z0 = 0.d0
	  avg_pns = 0.d0
	  avg_rhog = 0.d0
	  avg_zt = 0.d0
	  ! 16th percentile solution
	  zlams = zlams_or
	  z0 = z0temp
	  dz0 = epsdz0 * z0
	  den = denmax
	  zsfmin = taumin/((den - rho_air)*g*sin(rad(slope_ground)))
	  zshr = zsfmin		
	  write (flog, *) 'Pns max, rho_g min and ztot min calculation residuals'
	  write (*, *) 'Pns max, rho_g min and ztot min calculation residuals'
	  do while(z0.gt.1.d-10)
		 write(flog,*)'Solution with z0 = ', z0
		 write(*,*)'Solution with z0 = ', z0
		 x(1) = pnsmaxguess
		 x(2) = rhogminguess
		 x(3) = ztminguess
		 call newt(x, check, 3)
		 ! Double check the solutions obtained from newt
		 do i_x = 1, size(x)
			if (x(i_x).lt.0.d0.or.x(2).gt.rho_air) check = .true.
			! Spurious convergence case, i.e. check=false but non realistic solutions obtained
			if(check) exit
		 enddo	
		 if(check) then
			! No convergence case
			if(zlams.eq.UNDEFINED) then
			   z0 = z0 - dz0
			   write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
			   check = .false.
			else
			   z0 = 0.01d0
			   dz0 = epsdz0 * z0
			   write(flog,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   zlams = UNDEFINED
			   cycle
			endif
		 else
			avg_z0 = avg_z0 + z0
			avg_pns = avg_pns + x(1)
			avg_rhog = avg_rhog + x(2)
			avg_zt = avg_zt + x(3)
			n_solutions_z0 = n_solutions_z0 + 1
			if(zlams_or.eq.UNDEFINED) then
				z0 = z0 - dz0
				cycle
			else
				if(zlams.eq.UNDEFINED) then
					z0 = z0 - dz0
					cycle
				else
					exit
				endif
			endif
		 endif
	  enddo
	  if(check.or.n_solutions_z0.eq.0) then
		 write(flog,*)'Unable to find a solution for Pns max, rho_g min and ztot min'
		 write(*,*)'Unable to find a solution for Pns max, rho_g min and ztot min'
		 convergence = .false.
		 if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
		 return
	  endif
	  z0max = avg_z0 / n_solutions_z0
	  pnsmax = avg_pns / n_solutions_z0
	  rhogmin = avg_rhog / n_solutions_z0
	  ztmin = avg_zt / n_solutions_z0
	  write(*,362) pnsmax, rhogmin, ztmin, z0max
	  write(flog,362) pnsmax, rhogmin, ztmin, z0max
362   format('Pnsusp max =', f6.3, 1x, 'rhogmin =', f9.6, 1x, 'ztmin =', f8.3, 1x, 'z0max =', f9.6,/)	
	  cmax = (denmax - rhogmin)/(densp - rhogmin)
	  c_air_min = 1.d0 - c_gas_min
	  r_mix = c_gas_min * r_gas + c_air_min * r_air
	  t_mix_max = p_air / (r_mix * rhogmin)
	  cgastemp = c_gas_min
	  cairtemp = c_air_min
	  nfunc = 21
	  tmax = func(cmax)	 
	  n_solutions_z0 = 0
	  avg_z0 = 0.d0
	  avg_pns = 0.d0
   else
	  cavg = (dennrm - dengas)/(densp - dengas)
	  cmax = (denmax - dengas)/(densp - dengas)
	  cmin = (denmin - dengas)/(densp - dengas)
	  write (flog, 363) cavg, cmax, cmin
	  write (*, 363) cavg, cmax, cmin
363	  format('C 50th =', e12.5, 1x, 'C 84th =', e12.5, 1x, 'C 16th =', e12.5,/)
	  !     Average, maximum and minimum total flow thickness
	  ztavg = zlam/cavg
	  ztmax = zlam/cmin
	  ztmin = zlam/cmax
	  write (flog, 364) ztavg, ztmax, ztmin
	  write (*, 364) ztavg, ztmax, ztmin
364   format('Htot 50th =', f8.3, 1x, 'Htot 84th =', f8.3, 1x, 'Htot 16th =', f8.3,/)
	  if (slope_ground .eq. undefined) then
		 avg_z0 = 0.d0
		 avg_pns = 0.d0
		 avg_zsf = 0.d0
		 n_solutions_z0 = 0
		 nnewt = 1
		 !    Pnsusp,avg and slope angle
		 ztot = ztavg
		 den = dennrm
		 z0 = z0temp
		 dz0 = epsdz0 * z0
		 do while(z0.gt.1.d-10)
			write(flog,*)'Solution with z0guess = ', z0
			write(*,*)'Solution with z0guess = ', z0
			x(1) = pnsavgguess
			x(2) = zsfavgguess
			write (flog, *) 'Pns avg and zsf avg calculation residuals'
			write (*, *) 'Pns avg and zsf avg calculation residuals'
			call newt(x, check, 2)
		 ! Double check the solutions obtained from newt
			do i_x = 1, size(x)
			   if (x(i_x).lt.0.d0) check = .true.
			! Spurious convergence case, i.e. check=false but non realistic solutions obtained
			   if(check) exit
			enddo	
			if(check) then
			! No convergence case
			   if(zlams.eq.UNDEFINED) then
				  z0 = z0 - dz0
				  write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
				  write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
				  check = .false.
			   else
				   z0 = 0.01d0
				   dz0 = epsdz0 * z0
				   write(flog,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
				   write(*,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
				   zlams = UNDEFINED
				   cycle
			   endif
			else
				avg_pns = avg_pns + x(1)
				avg_zsf = avg_zsf + x(2)
				avg_z0 = avg_z0 + z0
				n_solutions_z0 = n_solutions_z0 + 1
				if(zlams_or.eq.UNDEFINED) then
					z0 = z0 - dz0
					cycle
				else
					if(zlams.eq.UNDEFINED) then
						z0 = z0 - dz0
						cycle
					else
						exit
					endif
				endif
			endif
		 enddo
		 if(check.or.n_solutions_z0.eq.0) then
			write(flog,*)'Unable to find a solution for z0avg and Pnsavg'
			write(*,*)'Unable to find a solution for z0avg and Pnsavg'
			convergence = .false.
			if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
			return
		 endif	
		 pnsavg = avg_pns / n_solutions_z0
		 zsfavg = avg_zsf / n_solutions_z0
		 z0avg = avg_z0 / n_solutions_z0		 
		 senx = tauavg/((den - rho_air)*g*zsfavg)
		 slope_ground = grad(asin(senx))
		 write (flog, 366) pnsavg, z0avg, slope_ground
		 write (*, 366) pnsavg, z0avg, slope_ground
366 	 format('Pnsusp avg =', f6.3, 1x, 'z0avg =', f9.6, 1x, 'Slope (ø) =', f6.3,/)
	  else
		 avg_z0 = 0.d0
		 avg_pns = 0.d0
		 n_solutions_z0 = 0
		 nnewt = 2
		 !     Pnsusp,avg and z0avg
		 z0 = z0temp
		 dz0 = epsdz0 * z0
		 den = dennrm
		 zsfavg = tauavg/((den - rho_air)*g*sin(rad(slope_ground)))
		 zshr = zsfavg
		 ztot = ztavg
		 write (flog, *) 'Pns avg and z0 avg calculation residuals'
		 write (*, *) 'Pns avg and z0 avg calculation residuals'
		 do while(z0.gt.1.d-10)
			write(flog,*)'Solution with z0guess = ', z0
			write(*,*)'Solution with z0guess = ', z0
			x(1) = pnsavgguess
			x(2) = z0
			call newt(x, check, 2)
		 ! Double check the solutions obtained from newt
			do i_x = 1, size(x)
			   if (x(i_x).lt.0.d0) check = .true.
			! Spurious convergence case, i.e. check=false but non realistic solutions obtained
			   if(check) exit
			enddo	
			if(check) then
			! No convergence case
			   if(zlams.eq.UNDEFINED) then
				  z0 = z0 - dz0
				  write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
				  write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
				  check = .false.
			   else
				   z0 = 0.01d0
				   dz0 = epsdz0 * z0
				   write(flog,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
				   write(*,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
				   zlams = UNDEFINED
				   cycle
			   endif
			else
				avg_pns = avg_pns + x(1)
				avg_z0 = avg_z0 + x(2)
				n_solutions_z0 = n_solutions_z0 + 1
				if(zlams_or.eq.UNDEFINED) then
					z0 = z0 - dz0
					cycle
				else
					if(zlams.eq.UNDEFINED) then
						z0 = z0 - dz0
						cycle
					else
						exit
					endif
				endif
			endif
		 enddo
		 pnsavg = avg_pns / n_solutions_z0
		 z0avg = avg_z0 / n_solutions_z0
		 if(check.or.n_solutions_z0.eq.0) then
			write(flog,*)'Unable to find a solution for z0avg and Pnsavg'
			write(*,*)'Unable to find a solution for z0avg and Pnsavg'
			convergence = .false.
			if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
			return
		 endif
		 write (flog, 365) pnsavg, z0avg
		 write (*, 365) pnsavg, z0avg
365 	 format('Pnsusp avg =', f6.3, 1x, 'z0avg =', f9.6,/)
	  end if
	  nnewt = 2
	  !     Pnsusp,max and z0max
	  n_solutions_z0 = 0
	  avg_z0 = 0.d0
	  avg_pns = 0.d0
	  zlams = zlams_or
	  z0 = z0temp
	  dz0 = epsdz0 * z0
	  den = denmax
	  zsfmin = taumin/((den - rho_air)*g*sin(rad(slope_ground)))
	  zshr = zsfmin
	  ztot = ztmin
	  write (flog, *) 'Pns max and z0 max calculation residuals'
	  write (*, *) 'Pns max and z0 max calculation residuals'
	  do while(z0.gt.1.d-10)
		 write(flog,*)'Solution with z0guess = ', z0
		 write(*,*)'Solution with z0guess = ', z0
		 x(1) = pnsmaxguess
		 x(2) = z0
		 call newt(x, check, 2)
		 ! Double check the solutions obtained from newt
		 do i_x = 1, size(x)
			if (x(i_x).lt.0.d0) check = .true.
			! Spurious convergence case, i.e. check=false but non realistic solutions obtained
			if(check) exit
		 enddo	
		 if(check) then
			! No convergence case
			if(zlams.eq.UNDEFINED) then
			   z0 = z0 - dz0
			   write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
			   check = .false.
			else
			   z0 = 0.01d0
			   dz0 = epsdz0 * z0
			   write(flog,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   zlams = UNDEFINED
			   cycle
			endif
		 else
			avg_pns = avg_pns + x(1)
			avg_z0 = avg_z0 + x(2)
			n_solutions_z0 = n_solutions_z0 + 1
			if(zlams_or.eq.UNDEFINED) then
				z0 = z0 - dz0
				cycle
			else
				if(zlams.eq.UNDEFINED) then
					z0 = z0 - dz0
					cycle
				else
					exit
				endif
			endif
		 endif
	  enddo
	  pnsmax = avg_pns / n_solutions_z0
	  z0max = avg_z0 / n_solutions_z0
	  if(check.or.n_solutions_z0.eq.0) then
		 write(flog,*)'Unable to find a solution for z0max and Pnsmax'
		 write(*,*)'Unable to find a solution for z0max and Pnsmax'
		 convergence = .false.
		 if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
		 return
	  endif
	  write (flog, 367) pnsmax, z0max
	  write (*, 367) pnsmax, z0max
367    format('Pnsusp max =', f6.3, 1x, 'z0max =', f9.6,/)
	  !     Pnsusp,min and z0min
	  n_solutions_z0 = 0
	  avg_z0 = 0.d0
	  avg_pns = 0.d0
	  zlams = zlams_or
	  z0 = z0temp
	  dz0 = epsdz0 * z0	  
	  den = denmin
	  zsfmax = taumax/((den - rho_air)*g*sin(rad(slope_ground)))
	  zshr = zsfmax
	  ztot = ztmax
	  write (flog, *) 'Pns min and z0 min calculation residuals'
	  write (*, *) 'Pns min and z0 min calculation residuals'
	  do while(z0.gt.1.d-10) 
		 write(flog,*)'Solution with z0guess = ', z0
		 write(*,*)'Solution with z0guess = ', z0
		 x(1) = pnsminguess
		 x(2) = z0
		 call newt(x, check, 2)
		 ! Double check the solutions obtained from newt
		 do i_x = 1, size(x)
			if (x(i_x).lt.0.d0) check = .true.
			! Spurious convergence case, i.e. check=false but non realistic solutions obtained
			if(check) exit
		 enddo	
		 if(check) then
			! No convergence case
			if(zlams.eq.UNDEFINED) then
			   z0 = z0 - dz0
			   write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
			   check = .false.
			else
			   z0 = 0.01d0
			   dz0 = epsdz0 * z0
			   write(flog,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   write(*,*)'Warning. Unable to converge with the provided z0. Setting z0 = ', z0
			   zlams = UNDEFINED
			   cycle
			endif
		 else
			avg_pns = avg_pns + x(1)
			avg_z0 = avg_z0 + x(2)
			n_solutions_z0 = n_solutions_z0 + 1
			if(zlams_or.eq.UNDEFINED) then
				z0 = z0 - dz0
				cycle
			else
				if(zlams.eq.UNDEFINED) then
					z0 = z0 - dz0
					cycle
				else
					exit
				endif
			endif
		 endif
	  enddo
	  pnsmin = avg_pns / n_solutions_z0
	  z0min = avg_z0 / n_solutions_z0
	  if(check.or.n_solutions_z0.eq.0) then
		 write(flog,*)'Unable to find a solution for z0min and Pnsmin'
		 write(*,*)'Unable to find a solution for z0min and Pnsmin'
		 convergence = .false.
		 if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
		 return
	  endif	  
	  write (flog, 368) pnsmin, z0min
	  write (*, 368) pnsmin, z0min
368   format('Pnsusp min =', f6.3, 1x, 'z0min =', f9.6,/)
	  rhogavg = dengas
	  rhogmax = dengas
	  rhogmin = dengas
   endif
   ! PROFILES
   open(53, file=trim(output_dir)//trim(path_sep)//'concentration_profile.dat')
   open(54, file=trim(output_dir)//trim(path_sep)//'dynamic_pressure_profile.dat')
   open(55, file=trim(output_dir)//trim(path_sep)//'velocity_profile.dat')
   open(56, file=trim(output_dir)//trim(path_sep)//'density_profile.dat')
   open(57, file=trim(output_dir)//trim(path_sep)//'temperature_profile.dat')
   z = zlam
   ! Average profiles
   write(53, 170) height, avgcon, maxcon, mincon
   write(54, 171) height, avgpd, maxpd, minpd
   write(55, 172) height, avgvel, maxvel, minvel
   write(56, 173) height, avgden, maxden, minden
   write(57, 174) height, avgtemp, maxtemp, mintemp
170 format(a10, 3(2x, a42))
171 format(a10, 3(2x, a37))
172 format(a10, 3(2x, a30))
173 format(a10, 3(2x, a33))
174 format(a10, 3(2x, a32))
   do while(z.le.ztavg)
	  if (z .gt. ztavg) then
		 cavg = 0.d0
	  else
		 cavg = c0*((z0avg/(ztavg - z0avg))*((ztavg - z)/z))**pnsavg
	  end if
	  denavg = cavg*densp + (1.d0 - cavg)*rhogavg
	  uavg = ushavg*((1.d0/kvk)*log(z/ks) + 8.5d0)
	  if (z .le. z0avg) cavg = c0
	  if (calc_t_mix) then
		 c_gas_avg = (rhogavg - rho_air) / (rho_gas - rho_air)
		 c_air_avg = 1.d0 - c_gas_avg
		 r_mix = c_gas_avg * r_gas + c_air_avg * r_air
		 t_mix_avg = p_air / (r_mix * rhogavg)
		 cgastemp = c_gas_avg
		 cairtemp = c_air_avg
		 nfunc = 21
		 t_mix_avg = func(cavg)
	  endif
	  pdynav = 0.5d0*denavg*uavg**2
	  !     Maximum profiles
	  if (z .gt. ztmin) then
		 cmax = 0.d0
	  else
		 cmax = c0*((z0max/(ztmin - z0max))*((ztmin - z)/z))**pnsmax
	  end if
	  dfmax = cmax*densp + (1.d0 - cmax)*rhogmin  ! rhogmin corresponds to the 84th percentile solution of Pdyn, since it comes from umax
	  umax = ushmax*((1.d0/kvk)*log(z/ks) + 8.5d0)
	  if (z .le. z0max) cmax = c0
	  if (calc_t_mix) then
		 c_gas_max = (rhogmax - rho_air) / (rho_gas - rho_air)
		 c_air_max = 1.d0 - c_gas_max
		 r_mix = c_gas_max * r_gas + c_air_max * r_air
		 t_mix_max = p_air / (r_mix * rhogmax) 
		 cgastemp = c_gas_min
		 cairtemp = c_air_min
		 nfunc = 21
		 t_mix_max = func(cmax)
	  endif
	  !     Minimum profiles
	  if (z .gt. ztmax) then
		 cmin = 0.d0
	  else
		 cmin = c0*((z0min/(ztmax - z0min))*((ztmax - z)/z))**pnsmin
	  end if
	  dfmin = cmin*densp + (1.d0 - cmin)*rhogmax ! rhogmax corresponds to the 16th percentile solution of Pdyn, since it comes from umin
	  umin = ushmin*((1.d0/kvk)*log(z/ks) + 8.5d0)
	  if (z .le. z0min) cmin = c0
	  if (calc_t_mix) then
		 c_gas_min = (rhogmin - rho_air) / (rho_gas - rho_air)
		 if(c_gas_min.gt.1.d0) then
			c_gas_min = 1.d0
			rho_gas = rhogmin
		 endif
		 c_air_min = 1.d0 - c_gas_min
		 r_mix = c_gas_min * r_gas + c_air_min * r_air
		 t_mix_min = p_air / (r_mix * rhogmin)
		 cgastemp = c_gas_max
		 cairtemp = c_air_max
		 nfunc = 21
		 t_mix_min = func(cmin)		
	  endif
	  pdynmx = 0.5d0*dfmin*umax**2
	  pdynmn = 0.5d0*dfmax*umin**2
	  write (53, 175) z, cavg, cmin, cmax
	  write (54, 176) z, pdynav, pdynmx, pdynmn
	  write (55, 177) z, uavg, umax, umin
	  write (56, 178) z, denavg, dfmin, dfmax
	  if (calc_t_mix) write (57, 178) z, t_mix_avg, t_mix_min, t_mix_max
175 format(1x, f8.3, 12x, e12.5, 2(33x, e12.5))
176 format(1x, f8.3, 12x, e12.4, 2(27x, e12.4))
177 format(1x, f8.3, 12x, f7.3, 2(26x, f7.3))
178 format(1x, f8.3, 12x, f8.3, 2(27x, f8.3))
	  z = z + dz
   enddo
   nfunc = 11
   !     Specific heights dynamic pressure
   if (ztmin .lt. 10.d0) then
	  write (*, *) 'Warning! Minimum total flow thickness is less than 10 m'
	  write (*, *) 'Calculating specific dynamic pressure over z less than ', ztmin
	  write (*, *) 'Otherwise free atmosphere will be taken into account'
	  write (*, *) ''
	  write (flog, *) 'Warning! Minimum total flow thickness is less than 10m'
	  write (flog, *) 'Calculating specific dynamic pressure over z less than ', ztmin
	  write (flog, *) 'Otherwise free atmosphere will be taken into account'
	  write (flog, *) ''
	  s = qsimp(z0avg, ztmin)
	  p10avg = (1.d0/(ztmin - z0avg))*s
	  nfunc = 12
	  s = qsimp(z0min, ztmin)
	  p10max = (1.d0/(ztmin - z0min))*s
	  nfunc = 13
	  s = qsimp(z0max, ztmin)
	  p10min = (1.d0/(ztmin - z0max))*s
   else
	  s = qsimp(z0avg, 10.d0)
	  p10avg = (1.d0/(10.d0 - z0avg))*s
	  nfunc = 12
	  s = qsimp(z0min, 10.d0)
	  p10max = (1.d0/(10.d0 - z0min))*s
	  nfunc = 13
	  s = qsimp(z0max, 10.d0)
	  p10min = (1.d0/(10.d0 - z0max))*s
   endif

   if (usr_z_dynpr) then
	  ipr = 0
	  do j = 1, 20
		if (zdynpr(j) .gt. ztmin) then
			zdynpr(j) = UNDEFINED
			pzavg(j) = UNDEFINED
			pzmax(j) = UNDEFINED
			pzmin(j) = UNDEFINED
			cycle
		else
			ipr = ipr + 1
			nfunc = 11
			s = qsimp(z0avg, zdynpr(j))
			pzavg(j) = (1.d0/(zdynpr(j) - z0avg))*s
			nfunc = 12
			s = qsimp(z0min, zdynpr(j))
			pzmax(j) = (1.d0/(zdynpr(j) - z0min))*s
			nfunc = 13
			s = qsimp(z0max, zdynpr(j))
			pzmin(j) = (1.d0/(zdynpr(j) - z0max))*s
		endif
	  enddo
	endif

   if (calc_t_mix) then
	  if (usr_z_t) then
		itemp = 0
		do j = 1, 20
			if (zt(j) .gt. ztmin) then
				zt(j) = UNDEFINED
				tzavg(j) = UNDEFINED
				tzmax(j) = UNDEFINED
				tzmin(j) = UNDEFINED
				cycle
			else
				itemp = itemp + 1
				nfunc = 22
				zttemp = ztavg
				z0temp = z0avg
				pnstemp = pnsavg
				cgastemp = c_gas_avg
				cairtemp = c_air_avg
				s = qsimp(z0temp, zt(j))
				tzavg(j) = (1.d0/(zt(j) - z0temp))*s
				zttemp = ztmin
				z0temp = z0max
				pnstemp = pnsmax
				cgastemp = c_gas_min
				cairtemp = c_air_min
				s = qsimp(z0temp, zt(j))
				tzmax(j) = (1.d0/(zt(j) - z0temp))*s
				zttemp = ztmax
				z0temp = z0min
				pnstemp = pnsmin
				cgastemp = c_gas_max
				cairtemp = c_air_max				
				s = qsimp(z0temp, zt(j))
				tzmin(j) = (1.d0/(zt(j) - z0temp))*s
			endif
		enddo
	  end if
   end if
	
   nfunc = 14
   c2avg = func(2.d0)
   s = qsimp(z0avg, 2.d0)
   c2dpavg = (1.d0/(2.d0 - z0avg))*s
   nfunc = 16
   c2min = func(2.d0)  ! 84th percentile general solution
   s = qsimp(z0min, 2.d0)
   c2dpmin = (1.d0/(2.d0 - z0max))*s
   nfunc = 15
   c2max = func(2.d0)  ! 16th percentile general solution
   s = qsimp(z0max, 2.d0)
   c2dpmax = (1.d0/(2.d0 - z0min))*s

   if (usr_z_c) then
		ic = 0
		do j = 1, 20
			if (zc(j) .gt. ztmin) then
				zc(j) = UNDEFINED
				czavg(j) = UNDEFINED
				czmax(j) = UNDEFINED
				czmin(j) = UNDEFINED
				czdpavg(j) = UNDEFINED
				czdpmax(j) = UNDEFINED
				czdpmin(j) = UNDEFINED
				cycle
			else
				ic = ic + 1
				nfunc = 14
				czavg(j) = func(zc(j))
				s = qsimp(z0avg, zc(j))
				czdpavg(j) = (1.d0/(zc(j) - z0avg))*s
				nfunc = 16
				czmin(j) = func(zc(j))  ! 84th percentile general solution
				s = qsimp(z0max, zc(j))
				czdpmin(j) = (1.d0/(zc(j) - z0max))*s
				czmax(j) = func(zc(j))  ! 16th percentile general solution
				s = qsimp(z0min, zc(j))
				czdpmax(j) = (1.d0/(zc(j) - z0min))*s
			endif
		enddo
   end if

   if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
   end subroutine profiles