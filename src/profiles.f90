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
	s, senx, uavg, umax, umin, dz0, epsdz0
	REAL :: rtemp(100),r(1001)
	integer :: i, j, i_r, i_x, mult
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
    mult = 0
    do i_r = 1, size(r), 2
        if(i_r.eq.1) then
            r(i_r) = 1.d0
            mult = 1
        else
            mult = mult + 1
            r(i_r) = -float(mult - 1)
            r(i_r - 1) = float(mult - 1)
        endif
    enddo
   !     Average, maximum and minimum shear flow particle concentration
	if (model .eq. 'TWOLAYERS') then
		densp = rhos(1, 0)
	else
		densp = rhos(2, 0)
	end if
	if (zlams .eq. UNDEFINED) then
		!Set zlam = 1 cm, the typical laminae thickness
		write (*, *) 'WARNING! Command ZLAMS missing in input.dat'
		write (*, *) 'Setting ZLAMS=1 cm'
		write (flog, *) 'WARNING! Command ZLAMS missing in input.dat'
		write (flog, *) 'Setting ZLAMS=1 cm'
		zlams = 0.01d0
		z0temp = zlams/10.d0
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
		do i_r = 1, size(r)
			x(1) = pnsavgguess
			x(2) = rhogavgguess
			x(3) = ztavgguess
			call newt(x, check, 3)
			! Double check the solutions obtained from newt
			do i_x = 1, size(x)
				if (x(i_x).lt.0.d0) check = .true.
				! Spurious convergence case, i.e. check=false but non realistic solutions obtained
				if(check) exit
			enddo	
			if(check) then
				! No convergence case
				z0 = z0 + r(i_r) * dz0
				!if(z0.lt.0.d0) z0=z0temp
				if(z0.le.0.d0.or.z0.gt.2*zlams) cycle
				write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
				write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
				if(i_r.ne.size(r)) check = .false.
				continue
			else
				exit
			endif
		enddo
		if(check) then
			write(flog,*)'Unable to find a solution for Pns avg, rho_g avg and ztot avg'
			write(*,*)'Unable to find a solution for Pns avg, rho_g avg and ztot avg'
			convergence = .false.
			if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
			return
			!call exit_pyflow
		endif
		z0avg = z0
		pnsavg = x(1)
		rhogavg = x(2)
		ztavg = x(3)	
		write(*,360) pnsavg, rhogavg, ztavg, z0avg
		write(flog,360) pnsavg, rhogavg, ztavg, z0avg
360     format('Pnsusp avg =', f6.3, 1x, 'rhogavg =', f9.6, 1x, 'ztavg =', f8.3, 1x, 'z0avg =', f9.6,/)
		cavg = (dennrm - rhogavg)/(densp - rhogavg)
		c_gas_avg = (rhogavg - rho_air) / (rho_gas - rho_air)
		c_air_avg = 1.d0 - c_gas_avg
		r_mix = c_gas_avg * r_gas + c_air_avg * r_air
		t_mix_avg = p_air / (r_mix * rhogavg)
		cgastemp = c_gas_avg
		cairtemp = c_air_avg
		nfunc = 21
		tavg = func(cavg)
		! 84th percentile solution
		z0 = z0temp
		dz0 = epsdz0 * z0
		den = denmin
		zsfmax = taumax/((den - rho_air)*g*sin(rad(slope_ground)))
		zshr = zsfmax
		write (flog, *) 'Pns min, rho_g max and ztot max calculation residuals'
		write (*, *) 'Pns min, rho_g max and ztot max calculation residuals'
		do i_r = 1, size(r)
			x(1) = pnsminguess
			x(2) = rhogmaxguess
			x(3) = ztmaxguess
			call newt(x, check, 3)
			! Double check the solutions obtained from newt
			do i_x = 1, size(x)
				if (x(i_x).lt.0.d0) check = .true.
				! Spurious convergence case, i.e. check=false but non realistic solutions obtained
				if(check) exit
			enddo		  
			if(check) then
				! No convergence case
				z0 = z0 + r(i_r) * dz0
				!if(z0.lt.0.d0) z0=z0temp
				if(z0.le.0.d0.or.z0.gt.2*zlams) cycle
				write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
				write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
				if(i_r.ne.size(r)) check = .false.
				continue
			else
				exit
			endif	
		enddo
		if(check) then
			write(flog,*)'Unable to find a solution for Pns min, rho_g max and ztot max'
			write(*,*)'Unable to find a solution for Pns min, rho_g max and ztot max'
			convergence = .false.
			if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
			return
			!call exit_pyflow
		endif
		z0min = z0
		pnsmin = x(1)
		rhogmax = x(2)
		ztmax = x(3)
		write(*,361) pnsmin, rhogmax, ztmax, z0min
		write(flog,361) pnsmin, rhogmax, ztmax, z0min
361     format('Pnsusp min =', f6.3, 1x, 'rhogmax =', f9.6, 1x, 'ztmax =', f8.3, 1x, 'z0min =', f9.6,/)		
		cmin = (denmin - rhogmax)/(densp - rhogmax) ! FABIO: ho cambiato in cmin per essere consistente
		c_gas_max = (rhogmax - rho_air) / (rho_gas - rho_air)
		c_air_max = 1.d0 - c_gas_max
		r_mix = c_gas_max * r_gas + c_air_max * r_air
		t_mix_min = p_air / (r_mix * rhogmax) 
		cgastemp = c_gas_max
		cairtemp = c_air_max
		nfunc = 21
		tmin = func(cmin) 
		! 16th percentile solution
		z0 = z0temp
		dz0 = epsdz0 * z0
		den = denmax
		zsfmin = taumin/((den - rho_air)*g*sin(rad(slope_ground)))
		zshr = zsfmin		
		write (flog, *) 'Pns max, rho_g min and ztot min calculation residuals'
		write (*, *) 'Pns max, rho_g min and ztot min calculation residuals'
		do i_r = 1, size(r)
			x(1) = pnsmaxguess
			x(2) = rhogminguess
			x(3) = ztminguess
			call newt(x, check, 3)
			! Double check the solutions obtained from newt
			do i_x = 1, size(x)
				if (x(i_x).lt.0.d0) check = .true.
				! Spurious convergence case, i.e. check=false but non realistic solutions obtained
				if(check) exit
			enddo	
			if(check) then
			! No convergence case
				z0 = z0 + r(i_r) * dz0
				!if(z0.lt.0.d0) z0=z0temp
				if(z0.le.0.d0.or.z0.gt.2*zlams) cycle
				write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0
				write(*,*)'Warning. Unable to converge. Setting z0 = ', z0
				if(i_r.ne.size(r)) check = .false.
				continue
			else
				exit
			endif
		enddo
		if(check) then
			write(flog,*)'Unable to find a solution for Pns max, rho_g min and ztot min'
			write(*,*)'Unable to find a solution for Pns max, rho_g min and ztot min'
			convergence = .false.
			if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
			return
			!call exit_pyflow
		endif
		z0max = z0
		pnsmax = x(1)
		rhogmin = x(2)
		ztmin = x(3)
		write(*,362) pnsmax, rhogmin, ztmin, z0max
		write(flog,362) pnsmax, rhogmin, ztmin, z0max
362     format('Pnsusp max =', f6.3, 1x, 'rhogmin =', f9.6, 1x, 'ztmin =', f8.3, 1x, 'z0max =', f9.6,/)	
		cmax = (denmax - rhogmin)/(densp - rhogmin)
		c_air_min = 1.d0 - c_gas_min
		r_mix = c_gas_min * r_gas + c_air_min * r_air
		t_mix_max = p_air / (r_mix * rhogmin)
		cgastemp = c_gas_min
		cairtemp = c_air_min
		nfunc = 21
		tmax = func(cmax)	 
	else
		cavg = (dennrm - dengas)/(densp - dengas)
		cmax = (denmax - dengas)/(densp - dengas)
		cmin = (denmin - dengas)/(densp - dengas)
		write (flog, 363) cavg, cmax, cmin
		write (*, 363) cavg, cmax, cmin
363 	format('C 50th =', e12.5, 1x, 'C 84th =', e12.5, 1x, 'C 16th =', e12.5,/)
		!     Average, maximum and minimum total flow thickness
		ztavg = zlam/cavg
		ztmax = zlam/cmin
		ztmin = zlam/cmax
		write (flog, 364) ztavg, ztmax, ztmin
		write (*, 364) ztavg, ztmax, ztmin
364 	  	format('Htot 50th =', f8.3, 1x, 'Htot 84th =', f8.3, 1x, 'Htot 16th =', f8.3,/)
		if (slope_ground .eq. undefined) then
			nnewt = 1
			!    Pnsusp,avg and slope angle
			ztot = ztavg
			den = dennrm
			x(1) = pnsavgguess
			x(2) = zsfavgguess
			write (flog, *) 'Pns avg and zsf avg calculation residuals'
			write (*, *) 'Pns avg and zsf avg calculation residuals'
			call newt(x, check, 2)
			pnsavg = x(1)
			zsfavg = x(2)
			z0avg = zlams
			senx = tauavg/((den - rho_air)*g*zsfavg)
			slope_ground = grad(asin(senx))
			write (flog, 366) pnsavg, z0avg, slope_ground
			write (*, 366) pnsavg, z0avg, slope_ground
366 	    format('Pnsusp avg =', f6.3, 1x, 'z0avg =', f9.6, 1x, 'Slope (ø) =', f6.3,/)
		else
			nnewt = 2
			!     Pnsusp,avg and z0avg
			den = dennrm
			zsfavg = tauavg/((den - rho_air)*g*sin(rad(slope_ground)))
			zshr = zsfavg
			ztot = ztavg
			dz0 = epsdz0 * z0avgguess
			do i_r = 1, size(r)
				x(1) = pnsavgguess
				x(2) = z0avgguess
				write (flog, *) 'Pns avg and z0 avg calculation residuals'
				write (*, *) 'Pns avg and z0 avg calculation residuals'
				call newt(x, check, 2)			
				! Double check the solutions obtained from newt
				do i_x = 1, size(x)
					if (x(i_x).lt.0.d0) check = .true.
					! Spurious convergence case, i.e. check=false but non realistic solutions obtained
					if(check) exit
				enddo	
				if(check) then
				! No convergence case
					z0avgguess = z0avgguess + r(i_r) * dz0
					!if(z0avgguess.lt.0.d0) z0avgguess = x(2)
					if(z0avgguess.le.0.d0.or.z0.gt.2*zlams) cycle
					write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0avgguess
					write(*,*)'Warning. Unable to converge. Setting z0 = ', z0avgguess
					if(i_r.ne.size(r)) check = .false.
					continue
				else
					exit
				endif				
			enddo
			pnsavg = x(1)
			z0avg = x(2)
			write (flog, 365) pnsavg, z0avg
			write (*, 365) pnsavg, z0avg
365 	    format('Pnsusp avg =', f6.3, 1x, 'z0avg =', f9.6,/)
		end if
		nnewt = 2
		!     Pnsusp,max and z0max
		den = denmax
		zsfmin = taumin/((den - rho_air)*g*sin(rad(slope_ground)))
		zshr = zsfmin
		ztot = ztmin
		dz0 = epsdz0 * z0maxguess
		write (flog, *) 'Pns max and z0 max calculation residuals'
		write (*, *) 'Pns max and z0 max calculation residuals'
		do i_r = 1, size(r)
			x(1) = pnsmaxguess
			x(2) = z0maxguess
			call newt(x, check, 2)
			! Double check the solutions obtained from newt
			do i_x = 1, size(x)
				if (x(i_x).lt.0.d0) check = .true.
				! Spurious convergence case, i.e. check=false but non realistic solutions obtained
				if(check) exit
			enddo	
			if(check) then
			! No convergence case
				z0maxguess = z0maxguess + r(i_r) * dz0
				!if(z0maxguess.lt.0.d0) z0maxguess = x(2)
				if(z0maxguess.le.0.d0.or.z0.gt.2*zlams) cycle
				write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0maxguess
				write(*,*)'Warning. Unable to converge. Setting z0 = ', z0maxguess
				if(i_r.ne.size(r)) check = .false.
				continue
			else
				exit
			endif
		enddo
		pnsmax = x(1)
		z0max = x(2)
		write (flog, 367) pnsmax, z0max
		write (*, 367) pnsmax, z0max
367 	format('Pnsusp max =', f6.3, 1x, 'z0max =', f9.6,/)
		!     Pnsusp,min and z0min
		den = denmin
		zsfmax = taumax/((den - rho_air)*g*sin(rad(slope_ground)))
		zshr = zsfmax
		ztot = ztmax
		dz0 = epsdz0 * z0minguess
		do i_r = 1, size(r)
			x(1) = pnsminguess
			x(2) = z0minguess
			write (flog, *) 'Pns min and z0 min calculation residuals'
			write (*, *) 'Pns min and z0 min calculation residuals'
			call newt(x, check, 2)
			! Double check the solutions obtained from newt
			do i_x = 1, size(x)
				if (x(i_x).lt.0.d0) check = .true.
				! Spurious convergence case, i.e. check=false but non realistic solutions obtained
				if(check) exit
			enddo	
			if(check) then
			! No convergence case
				z0minguess = z0minguess + r(i_r) * dz0
				!if(z0minguess.lt.0.d0) z0minguess = x(2)
				if(z0minguess.le.0.d0.or.z0.gt.2*zlams) cycle
				write(flog,*)'Warning. Unable to converge. Setting z0 = ', z0minguess
				write(*,*)'Warning. Unable to converge. Setting z0 = ', z0minguess
				if(i_r.ne.size(r)) check = .false.
				continue
			else
				exit
			endif
		enddo
		pnsmin = x(1)
		z0min = x(2)
		write (flog, 368) pnsmin, z0min
		write (*, 368) pnsmin, z0min
368   	format('Pnsusp min =', f6.3, 1x, 'z0min =', f9.6,/)
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
175 	format(1x, f8.3, 12x, e12.5, 2(33x, e12.5))
176 	format(1x, f8.3, 12x, e12.4, 2(27x, e12.4))
177 	format(1x, f8.3, 12x, f7.3, 2(26x, f7.3))
178 	format(1x, f8.3, 12x, f8.3, 2(27x, f8.3))
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
		nfunc = 11
		do j = 1, 20
			if (zdynpr(j) .eq. UNDEFINED) exit
			ipr = ipr + 1
			s = qsimp(z0avg, zdynpr(j))
			pzavg(j) = (1.d0/(zdynpr(j) - z0avg))*s
		end do
		nfunc = 12
		do j = 1, 20
			if (zdynpr(j) .eq. UNDEFINED) exit
			s = qsimp(z0min, zdynpr(j))
			pzmax(j) = (1.d0/(zdynpr(j) - z0min))*s
		end do
		nfunc = 13
		do j = 1, 20
			if (zdynpr(j) .eq. UNDEFINED) exit
			s = qsimp(z0max, zdynpr(j))
			pzmin(j) = (1.d0/(zdynpr(j) - z0max))*s
		end do
	end if

	if (calc_t_mix) then
		if (usr_z_t) then
			itemp = 0
			nfunc = 22
			do j = 1, 20
				if (zt(j) .eq. UNDEFINED) exit
				itemp = itemp + 1
				zttemp = ztavg
				z0temp = z0avg
				pnstemp = pnsavg
				cgastemp = c_gas_avg
				cairtemp = c_air_avg
				s = qsimp(z0temp, zt(j))
				tzav(j) = (1.d0/(zt(j) - z0temp))*s
			end do
			do j = 1, 20
				if (zt(j) .eq. UNDEFINED) exit
				zttemp = ztmin
				z0temp = z0max
				pnstemp = pnsmax
				cgastemp = c_gas_min
				cairtemp = c_air_min
				s = qsimp(z0temp, zt(j))
				tzmax(j) = (1.d0/(zt(j) - z0temp))*s
			end do
			do j = 1, 20
				if (zt(j) .eq. UNDEFINED) exit
				zttemp = ztmax
				z0temp = z0min
				pnstemp = pnsmin
				cgastemp = c_gas_max
				cairtemp = c_air_max				
				s = qsimp(z0temp, zt(j))
				tzmin(j) = (1.d0/(zt(j) - z0temp))*s
			end do
		end if
	end if
	
	nfunc = 14
	!c2av1 = func(2.d0)
	c2avg = func(2.d0)
	s = qsimp(z0avg, 2.d0)
	c2dpavg = (1.d0/(2.d0 - z0avg))*s
	nfunc = 16
	!c2max1 = func(2.d0)
	c2min = func(2.d0)  ! 84th percentile general solution
	s = qsimp(z0min, 2.d0)
	c2dpmin = (1.d0/(2.d0 - z0max))*s
	nfunc = 15
	c2max = func(2.d0)  ! 16th percentile general solution
	s = qsimp(z0max, 2.d0)
	c2dpmax = (1.d0/(2.d0 - z0min))*s
	!c2min1 = func(2.d0)

	if (usr_z_c) then
		ic = 0
		nfunc = 14
		do j = 1, 20
			if (zc(j) .eq. UNDEFINED) exit
			ic = ic + 1
			!czav1(j) = func(zc(j))
			czavg(j) = func(zc(j))
			s = qsimp(z0avg, zc(j))
			czdpavg = (1.d0/(zc(j) - z0avg))*s
		end do
		nfunc = 16
		do j = 1, ic
			if (zc(j) .eq. UNDEFINED) exit
			!czmax1(j) = func(zc(j))
			czmin(j) = func(zc(j))  ! 84th percentile general solution
			s = qsimp(z0max, zc(j))
			czdpmin = (1.d0/(zc(j) - z0max))*s
		end do
		nfunc = 15
		do j = 1, ic
			if (zc(j) .eq. UNDEFINED) exit
			!czmin1(j) = func(zc(j))
			czmax(j) = func(zc(j))  ! 16th percentile general solution
			s = qsimp(z0min, zc(j))
			czdpmax = (1.d0/(zc(j) - z0min))*s
		end do
	end if
	!     Real average, maximum and minimum solutions
	! arr(1) = p10av1
	! arr(2) = p10mx1
	! arr(3) = p10mn1
	! call piksrt(3, arr)
	! p10max = arr(3)
	! p10avg = arr(2)
	! p10min = arr(1)
	! do j = 1, ipr
		! arr(1) = pzav1(j)
		! arr(2) = pzmax1(j)
		! arr(3) = pzmin1(j)
		! call piksrt(3, arr)
		! pzmax(j) = arr(3)
		! pzavg(j) = arr(2)
		! pzmin(j) = arr(1)
	! end do
	
	! arr(1) = c2av1
	! arr(2) = c2max1
	! arr(3) = c2min1
	! call piksrt(3, arr)
	! c2max = arr(3)
	! c2avg = arr(2)
	! c2min = arr(1)
	! do j = 1, ic
		! arr(1) = czav1(j)
		! arr(2) = czmax1(j)
		! arr(3) = czmin1(j)
		! call piksrt(3, arr)
		! czmax(j) = arr(3)
		! czavg(j) = arr(2)
		! czmin(j) = arr(1)
	! end do
	
	! if(calc_t_mix) then
		! arr(1) = tavg
		! arr(2) = tmax
		! arr(3) = tmin
		! call piksrt(3, arr)
		! tmax = arr(3)
		! tavg = arr(2)
		! tmin = arr(1)
		! do j = 1, itemp
			! arr(1) = tzav1(j)
			! arr(2) = tzmax1(j)
			! arr(3) = tzmin1(j)
			! call piksrt(3, arr)
			! tzmax(j) = arr(3)
			! tzav(j) = arr(2)
			! tzmin(j) = arr(1)
		! end do
	! end if
	if (zlams_or .eq. UNDEFINED) zlams = UNDEFINED
    end subroutine profiles

    ! SUBROUTINE piksrt(n, arr)
        ! USE nrtype
        ! implicit none
        ! INTEGER :: n, i, j
        ! real(dp), dimension(n) :: arr
        ! real(dp) :: a
        ! do j = 2, n
            ! a = arr(j)
            ! do i = j - 1, 1, -1
               ! if (arr(i) .le. a) goto 10
               ! arr(i + 1) = arr(i)
            ! end do
            ! i = 0
! 10          arr(i + 1) = a
        ! end do
    ! end subroutine piksrt