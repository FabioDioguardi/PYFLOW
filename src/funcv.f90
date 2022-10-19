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
         real(dp) :: s
         real(dp) :: rhoptot, rhospst, ctottemp
         integer :: j
         select case (nnewt)
!     System of equation solved for Pnsusp (pnstmp) and zsf (zshr)

         case (1)
            pnstmp = x(1)
            zshr = x(2)
            funcv(1) = denatm - (densp*c0*((zlams/(ztot - zlams))*((ztot - zshr)/&
           &zshr))**pnstmp) - (dengas*(1.d0 - (c0*((zlams/(ztot - zlams))*((ztot&
           &- zshr)/zshr))**pnstmp)))
            nfunc = 9
            s = qsimp(zlams, zshr)
            funcv(2) = den - (1.d0/(zshr - zlams))*s
         case (2)
!     System of equations solved for Pnsusp (pns) and z0 (z0)
            pns = x(1)
            z0 = x(2)
            funcv(1) = denatm - dengas - (densp - dengas)*c0*((z0/(ztot - z0))*((ztot - zshr)/zshr))**pns
            nfunc = 10
            s = qsimp(z0, zshr)
            funcv(2) = den - (1.d0/(zshr - z0))*s
         case (3)
            pns = x(1)
            dengas = x(2)
            ztot = x(3)
			!write(*,*)pns, dengas, ztot
            funcv(1) = denatm - dengas - (densp - dengas)*c0*((z0/(ztot - z0))*((ztot - zshr)/zshr))**pns
			!write(*,*)funcv(1)
            nfunc = 10
			!write(*,*)z0, zshr, den
            s = qsimp(z0, zshr)
            funcv(2) = den - (1.d0/(zshr - z0))*s
			!write(*,*)funcv(2)
            funcv(3) = ztot - zlam/((den - dengas)/(densp - dengas))
			!write(*,*) pns, dengas, ztot, funcv(1), funcv(2), funcv(3)
            if (isnan(funcv(1)) .or. isnan(funcv(2))) stop

         end select
      end function funcv

