      subroutine twocomponent
         USE inoutdata; USE nrtype
         implicit none
         INTERFACE

            FUNCTION func(x)
               USE inoutdata; USE nrtype
               IMPLICIT NONE
               REAL(dp) :: func
               REAL(dp), INTENT(IN) :: x
            END FUNCTION func

            FUNCTION dlog2(xx)
               USE inoutdata; USE nrtype
               IMPLICIT NONE
               REAL(dp) :: dlog2
               REAL(dp), INTENT(IN) :: xx
            END FUNCTION dlog2

            FUNCTION zbrent(x1, x2)
               USE inoutdata; USE nrtype; USE nrutil, ONLY: nrerror
               IMPLICIT NONE
               REAL(dp), INTENT(IN) :: x1, x2
               REAL(dp) :: zbrent
            END FUNCTION zbrent

            FUNCTION qsimp(a, b)
               USE nrtype; USE nrutil, ONLY: nrerror
               REAL(dp), INTENT(IN) :: a, b
               REAL(dp) :: qsimp
            END FUNCTION qsimp

         END INTERFACE
         character(len=20) :: nmodel, cmd1
         REAL(dp) :: cdd2, cdd100, usq2, usq100, cddavg, denmod, cddnrm, aleft, aright, cddlft, cddrgt,&
        &usqmin, usqmax, cddmin, cddmax, cd2, d2mod, d2modm, s
         integer :: ntest, nrest, i
         nmodel = 'Two component method'

         nclasst = nclass(2)
         alfa = probt/2.d0
         d1m = d50mm(1)/1000.d0                                                  !First component diameter (m)
         d2m = d50mm(2)/1000.d0                                                   !Second component diameter (m)
         nfunc = 1
         ishape = 1
         jshape = 0
         s = qsimp(df2, df100)
         cd1 = (1.d0/(df100 - df2))*s                                          ! Cd avg first component
         cdd2 = (cd1*(rhos(2, 0) - df2))/(d1m*(rhos(1, 0) - df2))                          !Cd/d for flow density = 2 kg/m^3
         cdd100 = (cd1*(rhos(2, 0) - df100))/(d1m*(rhos(1, 0) - df100))                    !Cd/d for flow density = 100 kg/m^3
         write (52, 350) nmodel, cd1, cdd2, cdd100
         write (*, 350) nmodel, cd1, cdd2, cdd100
         usq2 = -(4.d0*d1m*g*(rhos(1, 0) - rhos(2, 0)))/(3.d0*(cdd2*d1m*rhos(1, 0) - cd1*rhos(2, 0)))       !ushear^2 for flow density = 2 kg/m^3
         usq100 = -(4.d0*d1m*g*(rhos(1, 0) - rhos(2, 0)))/(3.d0*(cdd100*d1m*rhos(1, 0) - cd1*rhos(2, 0)))   !ushear^2 for flow density = 100 kg/m^3
         write (52, 351) usq2, usq100
         write (*, 351) usq2, usq100
!     Cd/davg
         nfunc = 5
         s = qsimp(usq100, usq2)
         cddavg = (1.d0/(usq2 - usq100))*s
         nfunc = 6
         ishape = 2
         jshape = 0
         s = qsimp(df2, df100)
         cd2 = (1.d0/(df100 - df2))*s                                          ! Cd avg second component
         write (52, 361) cddavg, cd2
         write (*, 361) cddavg, cd2
!     Model diameter of the second component (m and phi units)
         d2modm = cd2/cddavg
         d2mod = -dlog2(d2modm*1000.d0)
         write (52, 362) d2modm, d2mod
         write (*, 362) d2modm, d2mod
!     t-Student test for comparing dmod and d of the second component
304      call testt(phi50(2), d2mod, sorting(2), nrest)
         write (*, 305) ttab, tcalc
305      format(/, 't tabulated =', f7.3, 2x, 't calculated =', f7.3,/)
         if (nrest .eq. 2) then
            write (*, 363) probt
            write (52, 363) probt
            probt = probt - 0.1d0*probt
            alfa = probt/2.d0
            write (*, 364) probt
            write (52, 364) probt
            goto 304
         else
            write (*, 365) probt
            write (52, 365) probt
         end if
!     Model flow density
         denmod = (cd1*rhos(2, 0) - cddavg*d1m*rhos(1, 0))/(cd1 - cddavg*d1m)
!     Normalized flow density
         dennrm = (-denmod*d2m*rhos(2, 0))/(d2modm*denmod - d2m*denmod - d2modm*rhos(2, 0))
!     Normalized Cd/d
         cddnrm = (cd1*(rhos(2, 0) - dennrm))/(d1m*(rhos(1, 0) - dennrm))
         write (52, 355) denmod, dennrm, cddnrm
         write (*, 355) denmod, dennrm, cddnrm
!     Normalized ushear^2
         fx = cddnrm
         usqnrm = (4.d0*g*(rhos(2, 0) - rhos(1, 0)))/(3.d0*(rhos(1, 0)*fx - (cd1*rhos(2, 0))/d1m))
         write (52, 356) usqnrm
         write (*, 356) usqnrm
         if (usqnrm .ge. usq2 .or. usqnrm .le. usq100) then
            write (*, *) 'Warning! ush^2 norm. is not in the range ush^2 (2 kg/m^3) - ush^2(100 kg/m^3)'
            write (*, *) 'This is due to the very low accuracy level of the T - Student test'
            write (52, *) 'Warning! ush^2 norm. is not in the range ush^2 (2 kg/m^3) - ush^2(100 kg/m^3)'
            write (52, *) 'This is due to the very low accuracy level of the T - Student test'
            write (*, *) 'Type "stop" if you want to stop the calculations'
            write (*, *) 'Any other input will continue the calculations'
            read (*, *) cmd1
            if (cmd1 .eq. 'stop') stop
         end if
!     Cd/d(34% left), Cd/d(34% right)
         nfunc = 5
         aleft = qsimp(usq100, usqnrm)
         cddlft = 0.68d0*aleft
         aright = qsimp(usqnrm, usq2)
         cddrgt = 0.68d0*aright
         write (52, 357) aleft, cddlft, aright, cddrgt
         write (*, 357) aleft, cddlft, aright, cddrgt
!     ushear^2 max and min
         nfunc1 = 2
         nfunc = 7
         fx = cddlft
         write (52, *) 'ush^2 min calculation residuals'
         write (*, *) 'ush^2 min calculation residuals'
         usqmin = zbrent(1.d-5, 600.d0)
         nfunc = 8
         fx = cddrgt
         write (52, *) ''
         write (52, *) 'ush^2 max calculation residuals'
         write (*, *) ''
         write (*, *) 'ush^2 max calculation residuals'
         usqmax = zbrent(1.d-5, 600.d0)
         write (52, 358) usqmin, usqmax
         write (*, 358) usqmin, usqmax
!     Cd/d(ushear^2 max) Cd/d(ushear^2 min)
         nfunc = 5
         cddmin = func(usqmin)
         cddmax = func(usqmax)
         write (52, 359) cddmin, cddmax
         write (*, 359) cddmin, cddmax
!     Flow densities as a function of Cd/d(ushear^2 max) Cd/d(ushear^2 min)
         denmax = (cd1*rhos(2, 0) - cddmin*d1m*rhos(1, 0))/(cd1 - cddmin*d1m)
         denmin = (cd1*rhos(2, 0) - cddmax*d1m*rhos(1, 0))/(cd1 - cddmax*d1m)
         write (52, 360) denmax, denmin
         write (*, 360) denmax, denmin
!     Shear stresses
         tauavg = dennrm*usqnrm
         taumin = denmax*usqmin
         taumax = denmin*usqmax
!     Shear velocities
         ushavg = sqrt(usqnrm)
         ushmin = sqrt(usqmin)
         ushmax = sqrt(usqmax)
350      format(a20, /, 'Cd1 = ', f8.3, /, 'Cd/d(2 kg/m^3) =', f12.4, 1x, 'Cd/d(100 kg/m^3) =', f12.4,/)
351      format('ush^2(2kg/m^3) =', f8.4, 1x, 'ush^2(100 kg/m^3) =', f8.4,/)
361      format('Cd/d avg =', f12.4, 1x, 'Cd2 =', f8.4,/)
362      format('d2mod (m) =', f8.4, 1x, 'd2mod (phi)', f8.4,/)
355      format('Denmod =', f8.4, 1x, 'Den norm =', f8.4, 1x, 'Cd/d norm =', f12.4,/)
356      format('ush^2 norm =', f8.4,/)
357      format('Aleft =', f10.4, 1x, 'Cd/d left =', f12.4, 1x, 'Aright =', f10.4, 1x, 'Cd/d right =', f12.4,/)
358      format('ush^2 min =', f8.4, 1x, 'ush^2 max =', f8.4,/)
359      format('Cd/d(ush^2 min) =', f12.4, 1x, 'Cd/d(ush^2 max) =', f12.4,/)
360      format('Den(Cd/d(ush^2 min)) =', f8.4, 1x, 'Den(Cd/d(ush^2 max)) =', f8.4,/)
363      format('t-test failed at the significance level = 'e8.2,/)
364      format('New significance level = 'e8.2,/)
365      format('Test t OK at a significance level = 'e8.2,/)
      end subroutine twocomponent

