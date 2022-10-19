      function func(x)
         USE inoutdata; USE nrtype
         implicit none
         INTERFACE
            FUNCTION tstud(x)
               USE inoutdata; USE nrtype
               IMPLICIT NONE
               REAL(dp) :: tstud
               REAL(dp), INTENT(IN) :: x
            END FUNCTION tstud

            FUNCTION cum(x)
               USE inoutdata; USE nrtype
               IMPLICIT NONE
               REAL(dp) :: cum
               REAL(dp), INTENT(IN) :: x
            END FUNCTION cum

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

            FUNCTION znorm(x)
               USE inoutdata; USE nrtype
               IMPLICIT NONE
               REAL(dp) :: znorm
               REAL(dp), INTENT(IN) :: x
            END FUNCTION znorm

         END INTERFACE
         REAL(dp), INTENT(IN) :: x
         REAL(dp) :: func, ftemp, cavg, cmin, cmax
         real(dp) :: denavg, dfmax, dfmin, uavg, umax, umin
         select case (nfunc)
         case (1) !twocomponent
            func = cdmodel(d1m, rhos(1, 0), x)
         case (2) !twolayer
            func = (4.d0*g*(rhos(1, 0) - (theta*g*dm_ent*dens_ent)/(x + theta*g*dm_ent)))/&
           &(3.d0*x*((theta*g*dm_ent*dens_ent)/(x + theta*g*dm_ent)))
         case (3) !twolayer
            func = cdmodel(d1m, rhos(1, 0), x)
         case (4) !testt
            func = alfa*2.d0 - tstud(x)
         case (5) !twocomponent
            func = (g*((3.d0*cd1*rhos(2, 0))/(g*d1m) + (4.d0*(rhos(2, 0) - rhos(1, 0)))/(x)))/(3.d0*rhos(1, 0))
         case (6) !twocomponent
            func = cdmodel(d2m, rhos(2, 0), x)
         case (7) !twolayer, twocomponent
            ftemp = qsimp(x, usqnrm)
            func = fx - ftemp
         case (8) !twolayer, twocomponent
            ftemp = qsimp(usqnrm, x)
            func = fx - ftemp
         case (9) !funcv
            func = densp*c0*((zlams/(ztot - zlams))*((ztot - x)/x))**pnstmp +&
           &dengas*(1.d0 - (c0*((zlams/(ztot - zlams))*((ztot - x)/x))**pnstmp))
         case (10) !funcv
            func = dengas + (densp - dengas)*c0*((z0/(ztot - z0))*((ztot - x)/x))**pns
			!write(*,*)x, dengas, pns, ztot, densp, func
         case (11) !profiles
            if (x .gt. ztavg) then
               cavg = 0.d0
               denavg = denatm
               uavg = 0.d0
            else
               cavg = c0*((z0avg/(ztavg - z0avg))*((ztavg - x)/x))**pnsavg
               denavg = cavg*densp + (1.d0 - cavg)*dengas
               uavg = ushavg*((1.d0/kvk)*log(x/ks) + 8.5d0)
            end if
            if (x .le. z0avg) cavg = c0
            func = 0.5d0*denavg*uavg**2
         case (12)  !profiles
            if (x .gt. ztmax) then
               cmin = 0.d0
               dfmin = denatm
               umax = 0.d0
            else
               cmin = c0*((z0min/(ztmax - z0min))*((ztmax - x)/x))**pnsmin
               dfmin = cmin*densp + (1.d0 - cmin)*dengas
               umax = ushmax*((1.d0/kvk)*log(x/ks) + 8.5d0)
            end if
            if (x .le. z0min) cmin = c0
            func = 0.5d0*dfmin*umax**2
         case (13)  !profiles
            if (x .gt. ztmin) then
               cmax = 0.d0
               dfmax = denatm
               umin = 0.d0
            else
               cmax = c0*((z0max/(ztmin - z0max))*((ztmin - x)/x))**pnsmax
               dfmax = cmax*densp + (1.d0 - cmax)*dengas
               umin = ushmin*((1.d0/kvk)*log(x/ks) + 8.5d0)
            end if
            if (x .le. z0max) cmax = c0
            func = 0.5d0*dfmax*umin**2
         case (14)  !profiles
            if (x .gt. ztot) then
               func = 0.d0
            else
               func = c0*((z0avg/(ztavg - z0avg))*((ztavg - x)/x))**pnsavg
            end if
            if (x .le. z0avg) func = c0
         case (15)  !profiles
            if (x .gt. ztmin) then
               func = 0.d0
            else
               func = c0*((z0max/(ztmin - z0max))*((ztmin - x)/x))**pnsmax
            end if
            if (x .le. z0max) func = c0
         case (16) !profiles
            if (x .gt. ztmax) then
               func = 0.d0
            else
               func = c0*((z0min/(ztmax - z0min))*((ztmax - x)/x))**pnsmin
            end if
            if (x .le. z0min) func = c0
         case (17)  !probfunction
            func = px - cum(x)
         case (18)  !probfunction
            func = mxdstr**x - mudstr**x - (mudstr**x - mndstr**x)
         case (19)
            func = fxdistr - cum(x)
         case (20)
            func = probchi - qsimp(x, xmaxchi)
         end select
      end function func

      function func1(x)
         USE inoutdata; USE nrtype
         implicit none
         INTERFACE

            FUNCTION znorm(x)
               USE inoutdata; USE nrtype
               IMPLICIT NONE
               REAL(dp) :: znorm
               REAL(dp), INTENT(IN) :: x
            END FUNCTION znorm

            FUNCTION chiqdr(x)
               USE inoutdata; USE nrtype
               IMPLICIT NONE
               REAL(dp) :: chiqdr
               REAL(dp), INTENT(IN) :: x
            END FUNCTION chiqdr

         END INTERFACE
         real(dp) :: x, func1
         select case (nfunc1)
         case (1) !twolayer
            func1 = (4.d0*g*(rhos(1, 0) - (theta*g*dm_ent*dens_ent)/(x + theta*g*dm_ent)))/&
           &(3.d0*x*((theta*g*dm_ent*dens_ent)/(x + theta*g*dm_ent)))
         case (2) !teocomponent
            func1 = (g*((3.d0*cd1*rhos(2, 0))/(g*d1m) + (4.d0*(rhos(2, 0) - rhos(1, 0)))/(x)))/(3.d0*rhos(1, 0))
         case (3)
            func1 = znorm(x)
         case (4)
            func1 = chiqdr(x)
         end select
      end function func1

      function cum(x)
         USE nrtype
         implicit none
         REAL(dp) :: cum
         REAL(dp), INTENT(IN) :: x
         cum = 0.5d0 + 0.5d0*erf(x/sqrt(2.d0))
      end function cum

      function chiqdr(x)
         USE inoutdata; USE nrtype
         implicit none
         REAL(dp) :: chiqdr
         REAL(dp), INTENT(IN) :: x
         chiqdr = (1.d0/((2.d0**(nfreedchi/2.d0))*dgamma(nfreedchi/2.d0)))*(x**(nfreedchi/2.d0 - 1.d0))*exp(-x/2.d0)
         return
      end

      function znorm(x)
         USE nrtype
         implicit none
         REAL(dp) :: znorm
         REAL(dp), INTENT(IN) :: x
         znorm = (1.d0/sqrt(2.d0*3.1416d0))*exp(-(x**2)/2.d0)
         return
      end

      function dlog2(xx)
         USE nrtype
         implicit none
         REAL(dp) :: dlog2
         REAL(dp), INTENT(IN) :: xx
         dlog2 = dlog(xx)/dlog(2.d0)
      end function dlog2

      function rad(xx)
         USE nrtype
         implicit none
         REAL(dp) :: rad
         REAL(dp), INTENT(IN) :: xx
         rad = ((2.d0*3.1416d0)/360.d0)*xx
      end function rad

      function grad(xx)
         USE nrtype
         implicit none
         REAL(dp) :: grad
         REAL(dp), INTENT(IN) :: xx
         grad = (360.d0/(2.d0*3.1416d0))*xx
      end function grad

