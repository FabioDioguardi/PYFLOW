      function funcv(x)
         USE inoutdata; USE nrtype
         implicit none
         INTERFACE

            FUNCTION qsimp(a, b)
               USE nrtype; USE nrutil, ONLY: nrerror
               REAL(dp), INTENT(IN) :: a, b
               REAL(dp) :: qsimp
            END FUNCTION qsimp

            FUNCTION cdmodel(dpart, rhop, rho)
               USE inoutdata; USE nrtype
               REAL(dp), INTENT(IN) :: dpart, rhop, rho
               REAL(dp) :: cdmodel
            END FUNCTION cdmodel

         END INTERFACE
         REAL(dp), DIMENSION(:), INTENT(IN) :: x
         REAL(dp), DIMENSION(size(x)) :: funcv
         real(dp) :: s, cg
         real(dp) :: rhoptot, rhospst, ctottemp
         integer :: j
         select case (nnewt)
!     System of equation solved for Pnsusp (pnstmp) and zsf (zshr)

         case (1)
            pnstmp = x(1)
            zshr = x(2)
            funcv(1) = rho_air - (densp*c0*((zlams/(ztot - zlams))*((ztot - zshr)/&
           &zshr))**pnstmp) - (dengas*(1.d0 - (c0*((zlams/(ztot - zlams))*((ztot&
           &- zshr)/zshr))**pnstmp)))
            nfunc = 9
            s = qsimp(zlams, zshr)
            funcv(2) = den - (1.d0/(zshr - zlams))*s
         case (2)
!     System of equations solved for Pnsusp (pns) and z0 (z0)
            pns = x(1)
            z0 = x(2)
            funcv(1) = rho_air - dengas - (densp - dengas)*c0*((z0/(ztot - z0))*((ztot - zshr)/zshr))**pns
            nfunc = 10
            s = qsimp(z0, zshr)
            funcv(2) = den - (1.d0/(zshr - z0))*s
         case (3)
            pns = x(1)
            dengas = x(2)
			if(x(3).le.zshr) then
				write(52,*)'Warning! Non realistic values encountered when calculating ztot'
				write(52,*)'Restarting the loop'
				write(*,*)'Warning! Non realistic values encountered when calculating ztot'
				write(*,*)'Restarting the loop'
				funcv(1) = 10.d0
				funcv(2) = 10.d0
				funcv(3) = 10.d0
				return
			endif
            ztot = x(3)
            funcv(1) = rho_air - dengas - (densp - dengas)*c0*((z0/(ztot - z0))*((ztot - zshr)/zshr))**pns
            nfunc = 10
            s = qsimp(z0, zshr)
            funcv(2) = den - (1.d0/(zshr - z0))*s
            funcv(3) = ztot - zlam/((den - dengas)/(densp - dengas))
            if (isnan(funcv(1)) .or. isnan(funcv(2))) stop
		
		case (4)
			write(123,*) ' I am here! ', r_gas, r_air, rhogavg, p_air
			funcv(1) = x(1) - x(2) * r_gas / r_air
			funcv(2) = rhogavg - ((x(3) * x(2)) / rhogavg) * x(2) - (1.d0 - ((x(3) * x(2)) / rhogavg)) * x(1)
			!funcv(2) = rhogavg - x(3) * x(2) - (1.d0 - x(3)) * x(1)
			funcv(3) = x(4) - ((x(3) * x(2)) / rhogavg) * r_gas - (1.d0 - ((x(3) * x(2)) / rhogavg)) * r_air
			!funcv(3) = x(4) - x(3) * r_gas - (1.d0 - x(3)) * r_air
			funcv(4) = x(2) - rhogavg  * x(4) / r_gas
         end select
      end function funcv

