      FUNCTION rtflsp(x1, x2)
         USE inoutdata; USE nrtype
         implicit none
         INTERFACE

            FUNCTION func(x)
               USE inoutdata; USE nrtype
               IMPLICIT NONE
               REAL(dp) :: func
               REAL(dp), INTENT(IN) :: x
            END FUNCTION func

         END INTERFACE
         integer, parameter :: MAXIT = 100000
         real(dp) :: rtflsp, x1, x2, del, dx, f, fh, fl, swap, xh, xl
         INTEGER :: j
         fl = func(x1)
         fh = func(x2)
!      if(fl*fh.gt.0.d0) pause 'root must be bracketed in rtflsp'
         if (fl .lt. 0.d0) then
            xl = x1
            xh = x2
         else
            xl = x2
            xh = x1
            swap = fl
            fl = fh
            fh = swap
         end if
         dx = xh - xl
         do j = 1, MAXIT
            rtflsp = xl + dx*fl/(fl - fh)
            f = func(rtflsp)
            if (f .lt. 0.d0) then
               del = xl - rtflsp
               xl = rtflsp
               fl = f
            else
               del = xh - rtflsp
               xh = rtflsp
               fh = f
            end if
            dx = xh - xl
            write (flog, 354) del, f
            write (*, 354) del, f
354         format(2(e10.3, 2x))
            if (abs(del) .lt. xacc .or. f .eq. 0.d0) return
         end do
!      pause 'rtflsp exceed maximum iterations'
      END

