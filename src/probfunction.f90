   subroutine probfunction
      USE inoutdata; USE nrtype
      implicit none
      INTERFACE
         FUNCTION zbrent(x1, x2)
            USE inoutdata; USE nrtype; USE nrutil, ONLY: nrerror
            IMPLICIT NONE
            REAL(dp) :: x1, x2
            REAL(dp) :: zbrent
         END FUNCTION zbrent
      END INTERFACE
      real(dp), dimension(20) :: mpz, mupz, sigpz, mcz, mdpcz, mucz, mudpcz, sigcz, sigdpcz, mtpz, mutpz, sigtpz
      real(dp) :: temp1, temp2, val, tempz1, tempz2, zstd, mp10, mc2, mdpc2, mup10, muc2, mudpc2, sigc2, sigdpc2, sigp10
      real(dp) :: mden, muden, sigden, mush, muush, sigush, mtp, mutp, sigtp
      real(dp) :: mrtot, mtdep, murtot, sigrtot, mutdep, sigtdep
      real(dp) :: x1tmp, x2tmp, x3tmp, x4tmp
      integer :: j, l, njc, njpr
      open(fout, file=trim(output_dir)//trim(path_sep)//'results.dat', position="append", status="old", action="write")
!      x1tmp=-50.d0;x2tmp=-1.d-10;x3tmp=1.d-10;x4tmp=50.d0
      x1tmp = -10.d0; x2tmp = -1.d-10; x3tmp = 1.d-10; x4tmp = 10.d0

      if (only_deprates .and. n_solutions .lt. 3) then
         write (*, *) 'WARNING! It is not possible to calculate PDFs of deposition rate and time'
         write (*, *) 'N_SOLUTIONS must be at least 3'
         write (flog, *) 'WARNING! It is not possible to calculate PDFs of deposition rate and time'
         write (flog, *) 'N_SOLUTIONS must be at least 3'
         return
      else
      end if
      if (only_deprates) goto 599
!     Symmetric probability distribution parameters for Pdyn at 10 m
      checkzbr = .false.
      mudstr = p10avg
      mxdstr = p10max
      mndstr = p10min
      nfunc = 18
      write (*, *) 'Pdyn 10 m simmetrization exponent calc. residuals'
      write (flog, *) 'Pdyn 10 m simmetrization exponent calc. residuals'
600   temp1 = zbrent(x3tmp, x4tmp)
      if (checkzbr) temp1 = 0.d0
      checkzbr = .false.
      write (flog, 369) temp1
      write (*, 369) temp1
      if (temp1 .ne. 0.d0 .and. temp1 .ge. 1.d-9 .and. temp1 .ne. 50.d0) then
         mp10 = temp1
      else
         temp2 = zbrent(x1tmp, x2tmp)
         if (checkzbr) then
            x1tmp = x1tmp - 50.d0
            x2tmp = x2tmp - 50.d0
            x3tmp = x3tmp + 50.d0
            x4tmp = x4tmp + 50.d0
            goto 600
         end if
         write (flog, 369) temp2
         write (*, 369) temp2
         if (temp2 .ne. 0.d0 .and. temp2 .le. -1.d-9 .and. temp2 .ne. -50.d0) then
            mp10 = temp2
         else
            write (*, *) 'Warning. Unable to find a simmetrization coefficient'
            write (*, *) 'for 10 m dynamic pressure prob. function'
         end if
      end if
      musim = mudstr**mp10
      sigsim = abs(mxdstr**mp10 - mudstr**mp10)
      mup10 = musim
      sigp10 = sigsim
      write (fout, *) '10 m dynamic pressure probability function'
      write (fout, 410) mp10, mup10, sigp10
      write (*, *) '10 m dynamic pressure probability function'
      write (*, 410) mp10, mup10, sigp10

!     Symmetric probability distribution parameters for C at 2 m
      checkzbr = .false.
      mudstr = c2avg
      mxdstr = c2max
      mndstr = c2min
      write (*, *) 'C 2 m simmetrization exponent calc. residuals'
      write (flog, *) 'C 2 m simmetrization exponent calc. residuals'
601   temp1 = zbrent(x3tmp, x4tmp)
      if (checkzbr) temp1 = 0.d0
      checkzbr = .false.
      write (flog, 369) temp1
      write (*, 369) temp1
      if (temp1 .ne. 0.d0 .and. temp1 .ge. 1.d-9 .and. temp1 .ne. 50.d0) then
         mc2 = temp1
      else
         temp2 = zbrent(x1tmp, x2tmp)
         if (checkzbr) then
            x1tmp = x1tmp - 50.d0
            x2tmp = x2tmp - 50.d0
            x3tmp = x3tmp + 50.d0
            x4tmp = x4tmp + 50.d0
            goto 601
         end if
         write (flog, 369) temp2
         write (*, 369) temp2
         if (temp2 .ne. 0.d0 .and. temp2 .le. -1.d-9 .and. temp2 .ne. -50.d0) then
            mc2 = temp2
         else
            write (*, *) 'Warning. Unable to find a simmetrization coefficient'
            write (*, *) 'for 2 m particle concentration prob. function'
         end if
      end if
      musim = mudstr**mc2
      sigsim = abs(mxdstr**mc2 - mudstr**mc2)
      muc2 = musim
      sigc2 = sigsim
      write (fout, *) '2 m particle concentration probability function'
      write (fout, 410) mc2, muc2, sigc2
      write (*, *) '2 m particle concentration probability function'
      write (*, 410) mc2, muc2, sigc2

!     Symmetric probability distribution parameters for 2 m depth-averaged C
      checkzbr = .false.
      mudstr = c2dpavg
      mxdstr = c2dpmax
      mndstr = c2dpmin
      write (*, *) '2 m depth-averaged C simmetrization exponent calc. residuals'
      write (flog, *) '2 m depth-averaged C simmetrization exponent calc. residuals'
602   temp1 = zbrent(x3tmp, x4tmp)
      if (checkzbr) temp1 = 0.d0
      checkzbr = .false.
      write (flog, 369) temp1
      write (*, 369) temp1
      if (temp1 .ne. 0.d0 .and. temp1 .ge. 1.d-9 .and. temp1 .ne. 50.d0) then
         mdpc2 = temp1
      else
         temp2 = zbrent(x1tmp, x2tmp)
         if (checkzbr) then
            x1tmp = x1tmp - 50.d0
            x2tmp = x2tmp - 50.d0
            x3tmp = x3tmp + 50.d0
            x4tmp = x4tmp + 50.d0
            goto 602
         end if
         write (flog, 369) temp2
         write (*, 369) temp2
         if (temp2 .ne. 0.d0 .and. temp2 .le. -1.d-9 .and. temp2 .ne. -50.d0) then
            mdpc2 = temp2
         else
            write (*, *) 'Warning. Unable to find a simmetrization coefficient'
            write (*, *) 'for 2 m depth-averaged particle concentration prob. function'
         end if
      end if
      musim = mudstr**mdpc2
      sigsim = abs(mxdstr**mdpc2 - mudstr**mdpc2)
      mudpc2 = musim
      sigdpc2 = sigsim
      write (fout, *) '2 m depth-averaged particle concentration probability function'
      write (fout, 410) mdpc2, mudpc2, sigdpc2
      write (*, *) '2 m depth-averaged particle concentration probability function'
      write (*, 410) mdpc2, mudpc2, sigdpc2


!     Symmetric probability distribution parameters for flow density
      checkzbr = .false.
      mudstr = dennrm
      mxdstr = denmax
      mndstr = denmin
      write (*, *) 'Density simmetrization exponent calc. residuals'
      write (flog, *) 'Density simmetrization exponent calc. residuals'
603   temp1 = zbrent(x3tmp, x4tmp)
      if (checkzbr) temp1 = 0.d0
      checkzbr = .false.
      write (flog, 369) temp1
      write (*, 369) temp1
      if (temp1 .ne. 0.d0 .and. temp1 .ge. 1.d-9 .and. temp1 .ne. 50.d0) then
         mden = temp1
      else
         temp2 = zbrent(x1tmp, x2tmp)
         if (checkzbr) then
            x1tmp = x1tmp - 50.d0
            x2tmp = x2tmp - 50.d0
            x3tmp = x3tmp + 50.d0
            x4tmp = x4tmp + 50.d0
            goto 603
         end if
         write (flog, 369) temp2
         write (*, 369) temp2
         if (temp2 .ne. 0.d0 .and. temp2 .le. -1.d-9 .and. temp2 .ne. -50.d0) then
            mden = temp2
         else
            write (*, *) 'Warning. Unable to find a simmetrization coefficient'
            write (*, *) 'for flow density prob. function'
         end if
      end if
      musim = mudstr**mden
      sigsim = abs(mxdstr**mden - mudstr**mden)
      muden = musim
      sigden = sigsim
      write (fout, *) 'Flow density probability function'
      write (fout, 410) mden, muden, sigden
      write (*, *) 'Flow density probability function'
      write (*, 410) mden, muden, sigden

!     Symmetric probability distribution parameters for flow shear velocity
      checkzbr = .false.
      mudstr = ushavg
      mxdstr = ushmax
      mndstr = ushmin
      write (*, *) 'Shear velocity simmetrization exponent calc. residuals'
      write (flog, *) 'Shear velocity simmetrization exponent calc. residuals'
604   temp1 = zbrent(x3tmp, x4tmp)
      if (checkzbr) temp1 = 0.d0
      checkzbr = .false.
      write (flog, 369) temp1
      write (*, 369) temp1
      if (temp1 .ne. 0.d0 .and. temp1 .ge. 1.d-9 .and. temp1 .ne. 50.d0) then
         mush = temp1
      else
         temp2 = zbrent(x1tmp, x2tmp)
         if (checkzbr) then
            x1tmp = x1tmp - 50.d0
            x2tmp = x2tmp - 50.d0
            x3tmp = x3tmp + 50.d0
            x4tmp = x4tmp + 50.d0
            goto 604
         end if
         write (flog, 369) temp2
         write (*, 369) temp2
         if (temp2 .ne. 0.d0 .and. temp2 .le. -1.d-9 .and. temp2 .ne. -50.d0) then
            mush = temp2
         else
            write (*, *) 'Warning. Unable to find a simmetrization coefficient'
            write (*, *) 'for flow shear velocity prob. function'
         end if
      end if
      musim = mudstr**mush
      sigsim = abs(mxdstr**mush - mudstr**mush)
      muush = musim
      sigush = sigsim
      write (fout, *) 'Flow shear velocity probability function'
      write (fout, 410) mush, muush, sigush
      write (*, *) 'Flow shear velocity probability function'
      write (*, 410) mush, muush, sigush

      if(calc_t_mix) then
!     Symmetric probability distribution parameters for flow temperature
         checkzbr = .false.
         mudstr = tavg
         mxdstr = tmax
         mndstr = tmin
         write (*, *) 'Temperature simmetrization exponent calc. residuals'
         write (flog, *) 'Temperature simmetrization exponent calc. residuals'
605      temp1 = zbrent(x3tmp, x4tmp)
         if (checkzbr) temp1 = 0.d0
         checkzbr = .false.
         write (flog, 369) temp1
         write (*, 369) temp1
         if (temp1 .ne. 0.d0 .and. temp1 .ge. 1.d-9 .and. temp1 .ne. 50.d0) then
            mtp = temp1
         else
            temp2 = zbrent(x1tmp, x2tmp)
            if (checkzbr) then
               x1tmp = x1tmp - 50.d0
               x2tmp = x2tmp - 50.d0
               x3tmp = x3tmp + 50.d0
               x4tmp = x4tmp + 50.d0
               goto 605
            end if
            write (flog, 369) temp2
            write (*, 369) temp2
            if (temp2 .ne. 0.d0 .and. temp2 .le. -1.d-9 .and. temp2 .ne. -50.d0) then
               mtp = temp2
            else
               write (*, *) 'Warning. Unable to find a simmetrization coefficient'
               write (*, *) 'for flow shear velocity prob. function'
            end if
         end if
         musim = mudstr**mtp
         sigsim = abs(mxdstr**mtp - mudstr**mtp)
         mutp = musim
         sigtp = sigsim
         write (fout, *) 'Flow temperature probability function'
         write (fout, 410) mtp, mutp, sigtp
         write (*, *) 'Flow temperature probability function'
         write (*, 410) mtp, mutp, sigtp
         if (usr_z_t) then
            write (*, *) '###PROBABILITY FUNCTIONS FOR DEPTH-AVERAGED FLOW TEMPERATURE OVER USER REQUESTED HEIGHTS###'
            write (flog, *) '###PROBABILITY FUNCTIONS FOR DEPTH-AVERAGED FLOW TEMPERATURE OVER USER REQUESTED HEIGHTS###'
            if (itemp .eq. 0) then
               write (*, *) 'WARNING! No valid ZT to calculate the depth-averaged flow temperature at user requested height probability function'
               write (*, *) ''
               write (flog, *) 'WARNING! No valid ZT to calculate the depth-averaged flow temperature at user requested height probability function'
               write (flog, *) ''
            else
            end if
            do j = 1, 20
               if(zt(j).eq.undefined) cycle
               checkzbr = .false.
               mudstr = tzavg(j)
               mxdstr = tzmax(j)
               mndstr = tzmin(j)
               write (*, 425) zt(j)
               write (flog, 425) zt(j)
606            tempz1 = zbrent(x3tmp, x4tmp)
               if (checkzbr) tempz1 = 0.d0
               checkzbr = .false.
               write (flog, 369) tempz1
               write (*, 369) tempz1
               if (tempz1 .ne. 0.d0 .and. tempz1 .ge. 1.d-9 .and. tempz1 .ne. 50.d0) then
                  mtpz(j) = tempz1
               else
                  tempz2 = zbrent(x1tmp, x2tmp)
                  if (checkzbr) then
                     x1tmp = x1tmp - 50.d0
                     x2tmp = x2tmp - 50.d0
                     x3tmp = x3tmp + 50.d0
                     x4tmp = x4tmp + 50.d0
                     goto 606
                  end if
                  write (flog, 369) tempz2
                  write (*, 369) tempz2
                  if (tempz2 .ne. 0.d0 .and. tempz2 .le. -1.d-9 .and. tempz2 .ne. -50.d0) then
                     mtpz(j) = tempz2
                  else
                     write (*, *) 'Warning. Unable to find a simmetrization coefficient'
                     write (*, *) 'for flow temperature prob. function at z=', zt(j)
                     write (flog, *) 'Warning. Unable to find a simmetrization coefficient'
                     write (flog, *) 'for flow temperature prob. function at z=', zt(j)
                  end if
               end if
               musim = mudstr**mtpz(j)
               sigsim = abs(mxdstr**mtpz(j) - mudstr**mtpz(j))
               mutpz(j) = musim
               sigtpz(j) = sigsim
               write (fout, 413) zt(j), mtpz(j), mutpz(j), sigtpz(j)
               write (*, 413) zt(j), mtpz(j), mutpz(j), sigtpz(j)
            end do
         end if
      endif

      if (usr_z_dynpr) then
!     Determination of probability function of Pdyn at user requested heights
         write (*, *) '###PROBABILITY FUNCTIONS FOR DEPTH-AVERAGED DYNAMIC PRESSURE OVER USER REQUESTED HEIGHTS###'
         write (flog, *) '###PROBABILITY FUNCTIONS FOR DEPTH-AVERAGED DYNAMIC PRESSURE OVER USER REQUESTED HEIGHTS###'
         if (ipr .eq. 0) then
            write (*, *) 'WARNING! No valid ZDYNPR to calculate the depth-averaged dynamic pressure at user requested height probability function'
            write (*, *) ''
            write (flog, *) 'WARNING! No valid ZDYNPR to calculate the depth-averaged dynamic pressure at user requested height probability function'
            write (flog, *) ''
         else
         end if
         do j = 1, 20
            if(zdynpr(j).eq.undefined) cycle
            checkzbr = .false.
            mudstr = pzavg(j)
            mxdstr = pzmax(j)
            mndstr = pzmin(j)
            write (*, 425) zdynpr(j)
            write (flog, 425) zdynpr(j)
607         tempz1 = zbrent(x3tmp, x4tmp)
            if (checkzbr) tempz1 = 0.d0
            checkzbr = .false.
            write (flog, 369) tempz1
            write (*, 369) tempz1
            if (tempz1 .ne. 0.d0 .and. tempz1 .ge. 1.d-9 .and. tempz1 .ne. 50.d0) then
               mpz(j) = tempz1
            else
               tempz2 = zbrent(x1tmp, x2tmp)
               if (checkzbr) then
                  x1tmp = x1tmp - 50.d0
                  x2tmp = x2tmp - 50.d0
                  x3tmp = x3tmp + 50.d0
                  x4tmp = x4tmp + 50.d0
                  goto 607
               end if
               write (flog, 369) tempz2
               write (*, 369) tempz2
               if (tempz2 .ne. 0.d0 .and. tempz2 .le. -1.d-9 .and. tempz2 .ne. -50.d0) then
                  mpz(j) = tempz2
               else
                  write (*, *) 'Warning. Unable to find a simmetrization coefficient'
                  write (*, *) 'for dynamic pressure prob. function at z=', zdynpr(j)
                  write (flog, *) 'Warning. Unable to find a simmetrization coefficient'
                  write (flog, *) 'for dynamic pressure prob. function at z=', zdynpr(j)
               end if
            end if
            musim = mudstr**mpz(j)
            sigsim = abs(mxdstr**mpz(j) - mudstr**mpz(j))
            mupz(j) = musim
            sigpz(j) = sigsim
            write (fout, 411) zdynpr(j), mpz(j), mupz(j), sigpz(j)
            write (*, 411) zdynpr(j), mpz(j), mupz(j), sigpz(j)
         end do
      end if

      if (usr_z_c) then
!     Determination of probability function of C at user requested heights
         write (*, *) '###PROBABILITY FUNCTIONS FOR PARTICLE CONCENTRATION AT USER REQUESTED HEIGHTS###'
         write (flog, *) '###PROBABILITY FUNCTIONS FOR PARTICLE CONCENTRATION AT USER REQUESTED HEIGHTS###'
         if (ic .eq. 0) then
            write (*, *) 'WARNING! No valid ZC to calculate the particle concentration at user requested height probability function'
            write (*, *) ''
            write (flog, *) 'WARNING! No valid ZC to calculate the particle concentration at user requested height probability function'
            write (flog, *) ''
         else
         end if
         do j = 1, 20
            if(zc(j).eq.undefined) cycle
            checkzbr = .false.
            mudstr = czavg(j)
            mxdstr = czmax(j)
            mndstr = czmin(j)
            write (*, 426) zc(j)
            write (flog, 426) zc(j)
608         tempz1 = zbrent(x3tmp, x4tmp)
            if (checkzbr) tempz1 = 0.d0
            checkzbr = .false.
            write (flog, 369) tempz1
            write (*, 369) tempz1
            if (tempz1 .ne. 0.d0 .and. tempz1 .ge. 1.d-9 .and. tempz1 .ne. 50.d0) then
               mcz(j) = tempz1
            else
               tempz2 = zbrent(x1tmp, x2tmp)
               if (checkzbr) then
                  x1tmp = x1tmp - 50.d0
                  x2tmp = x2tmp - 50.d0
                  x3tmp = x3tmp + 50.d0
                  x4tmp = x4tmp + 50.d0
                  goto 608
               end if
               write (flog, 369) tempz2
               write (*, 369) tempz2
               if (tempz2 .ne. 0.d0 .and. tempz2 .le. -1.d-9 .and. tempz2 .ne. -50.d0) then
                  mcz(j) = tempz2
               else
                  write (*, *) 'Warning. Unable to find a simmetrization coefficient'
                  write (*, *) 'for particle concentration prob. function at z=', zc(j)
                  write (flog, *) 'Warning. Unable to find a simmetrization coefficient'
                  write (flog, *) 'for particle concentration prob. function at z=', zc(j)
               end if
            end if
            musim = mudstr**mcz(j)
            sigsim = abs(mxdstr**mcz(j) - mudstr**mcz(j))
            mucz(j) = musim
            sigcz(j) = sigsim
            write (fout, 412) zc(j), mcz(j), mucz(j), sigcz(j)
            write (*, 412) zc(j), mcz(j), mucz(j), sigcz(j)
         end do
         do j = 1, 20
            if(zc(j).eq.undefined) cycle
            checkzbr = .false.
            mudstr = czdpavg(j)
            mxdstr = czdpmax(j)
            mndstr = czdpmin(j)
            write (*, 426) zc(j)
            write (flog, 426) zc(j)
609         tempz1 = zbrent(x3tmp, x4tmp)
            if (checkzbr) tempz1 = 0.d0
            checkzbr = .false.
            write (flog, 369) tempz1
            write (*, 369) tempz1
            if (tempz1 .ne. 0.d0 .and. tempz1 .ge. 1.d-9 .and. tempz1 .ne. 50.d0) then
               mdpcz(j) = tempz1
            else
               tempz2 = zbrent(x1tmp, x2tmp)
               if (checkzbr) then
                  x1tmp = x1tmp - 50.d0
                  x2tmp = x2tmp - 50.d0
                  x3tmp = x3tmp + 50.d0
                  x4tmp = x4tmp + 50.d0
                  goto 609
               end if
               write (flog, 369) tempz2
               write (*, 369) tempz2
               if (tempz2 .ne. 0.d0 .and. tempz2 .le. -1.d-9 .and. tempz2 .ne. -50.d0) then
                  mcz(j) = tempz2
               else
                  write (*, *) 'Warning. Unable to find a simmetrization coefficient'
                  write (*, *) 'for z = ', zc(j), ' depth-averaged particle concentration prob. function'
                  write (flog, *) 'Warning. Unable to find a simmetrization coefficient'
                  write (flog, *) 'for z = ', zc(j), ' depth-averaged particle concentration prob. function'
               end if
            end if
            musim = mudstr**mdpcz(j)
            sigsim = abs(mxdstr**mdpcz(j) - mudstr**mdpcz(j))
            mudpcz(j) = musim
            sigdpcz(j) = sigsim
            write (fout, 412) zc(j), mdpcz(j), mudpcz(j), sigdpcz(j)
            write (*, 412) zc(j), mdpcz(j), mudpcz(j), sigdpcz(j)
         end do
      end if

!     Symmetric probability distribution parameters for Rtot
599      if (deprates) then
         checkzbr = .false.
         mudstr = rtot_avg
         mxdstr = rtot_max
         mndstr = rtot_min
         nfunc = 18
         write (*, *) 'Rtot simmetrization exponent calc. residuals'
         write (flog, *) 'Rtot simmetrization exponent calc. residuals'
610      temp1 = zbrent(x3tmp, x4tmp)
         if (checkzbr) temp1 = 0.d0
         checkzbr = .false.
         write (flog, 369) temp1
         write (*, 369) temp1
         if (temp1 .ne. 0.d0 .and. temp1 .ge. 1.d-9 .and. temp1 .ne. 50.d0) then
            mrtot = temp1
         else
            temp2 = zbrent(x1tmp, x2tmp)
            if (checkzbr) then
               x1tmp = x1tmp - 50.d0
               x2tmp = x2tmp - 50.d0
               x3tmp = x3tmp + 50.d0
               x4tmp = x4tmp + 50.d0
               goto 610
            end if
            write (flog, 369) temp2
            write (*, 369) temp2
            if (temp2 .ne. 0.d0 .and. temp2 .le. -1.d-9 .and. temp2 .ne. -50.d0) then
               mrtot = temp2
            else
               write (*, *) 'Warning. Unable to find a simmetrization coefficient'
               write (*, *) 'for 10 m dynamic pressure prob. function'
            end if
         end if
         musim = mudstr**mrtot
         sigsim = abs(mxdstr**mrtot - mudstr**mrtot)
         murtot = musim
         sigrtot = sigsim
         write (fout, *) 'Rtot probability function'
         write (fout, 410) mrtot, murtot, sigrtot
         write (*, *) 'Rtot probability function'
         write (*, 410) mrtot, murtot, sigrtot

!     Symmetric probability distribution parameters for tdep
         checkzbr = .false.
         mudstr = tdep_avg
         mxdstr = tdep_max
         mndstr = tdep_min
         nfunc = 18
         write (*, *) 'tdep simmetrization exponent calc. residuals'
         write (flog, *) 'tdep simmetrization exponent calc. residuals'
611      temp1 = zbrent(x3tmp, x4tmp)
         if (checkzbr) temp1 = 0.d0
         checkzbr = .false.
         write (flog, 369) temp1
         write (*, 369) temp1
         if (temp1 .ne. 0.d0 .and. temp1 .ge. 1.d-9 .and. temp1 .ne. 50.d0) then
            mtdep = temp1
         else
            temp2 = zbrent(x1tmp, x2tmp)
            if (checkzbr) then
               x1tmp = x1tmp - 50.d0
               x2tmp = x2tmp - 50.d0
               x3tmp = x3tmp + 50.d0
               x4tmp = x4tmp + 50.d0
               goto 611
            end if
            write (flog, 369) temp2
            write (*, 369) temp2
            if (temp2 .ne. 0.d0 .and. temp2 .le. -1.d-9 .and. temp2 .ne. -50.d0) then
               mtdep = temp2
            else
               write (*, *) 'Warning. Unable to find a simmetrization coefficient'
               write (*, *) 'for 10 m dynamic pressure prob. function'
            end if
         end if
         musim = mudstr**mtdep
         sigsim = abs(mxdstr**mtdep - mudstr**mtdep)
         mutdep = musim
         sigtdep = sigsim
         write (fout, *) 'tdep probability function'
         write (fout, 410) mtdep, mutdep, sigtdep
         write (*, *) 'tdep probability function'
         write (*, 410) mtdep, mutdep, sigtdep
      else
      end if

!     Calculation of function values at a desired percentile
      if (.not. usr_pcx_sol) return
      write (*, *) '###FUNCTIONS VALUES AT USER REQUESTED PERCENTILES###'
      write (fout, *) '###FUNCTIONS VALUES AT USER REQUESTED PERCENTILES###'
      write (flog, *) '###FUNCTIONS VALUES AT USER REQUESTED PERCENTILES###'
      l = 1
      do while (pcx(l) .ne. undefined)
         if (pcx(l) .lt. 0.d0 .or. pcx(l) .gt. 1.d0) then
            write (*, 428) l
            write (flog, 428) l
            write (*, *) 'PYFLOW 2.5 is going to be stopped'
            write (flog, *) 'PYFLOW 2.5 is going to be stopped'
            stop
         else
            nfunc = 17
         end if
         write (*, 427) pcx(l)
         write (fout, 427) pcx(l)
         write (flog, 427) pcx(l)
         if (only_deprates) goto 404
         write (*, *) 'Dynamic pressure 10 m'
         write (flog, *) 'Dynamic pressure 10 m'
         musim = mup10
         sigsim = sigp10
         mm = mp10
         if (mm .lt. 0.d0) then
            px = 1.d0 - pcx(l)
         else
            px = pcx(l)
         end if
         zstd = zbrent(-4.d0, 4.d0)
         val = zstd*sigsim + musim
         if (val .le. 0.d0) then
            write (*, *) 'Warning!!'
            write (*, *) 'The percentile is outside the range of calculation'
            write (flog, *) 'Warning!!'
            write (flog, *) 'The percentile is outside the range of calculation'
         else
            val = val**(1.d0/mp10)
            write (*, *) 'Dynamic pressure 10 m'
            write (*, 420) pcx(l), val
            write (fout, *) 'Dynamic pressure 10 m'
            write (fout, 420) pcx(l), val
         end if
         write (*, *) 'Particle concentration 2 m'
         write (flog, *) 'Particle concentration 2 m'
         musim = muc2
         sigsim = sigc2
         mm = mc2
         if (mm .lt. 0.d0) then
            px = 1.d0 - pcx(l)
         else
            px = pcx(l)
         end if
         zstd = zbrent(-4.d0, 4.d0)
         val = zstd*sigsim + musim
         if (val .le. 0.d0) then
            write (*, *) 'Warning!!'
            write (*, *) 'The percentile is outside the range of calculation'
            write (flog, *) 'Warning!!'
            write (flog, *) 'The percentile is outside the range of calculation'
         else
            val = val**(1.d0/mc2)
            write (*, *) 'Particle concentration 2 m'
            write (*, 421) pcx(l), val
            write (fout, *) 'Particle concentration 2 m'
            write (fout, 421) pcx(l), val
         end if
         write (*, *) '2 m depth-averaged particle concentration'
         write (flog, *) '2 m depth-averaged particle concentration'
         musim = mudpc2
         sigsim = sigdpc2
         mm = mdpc2
         if (mm .lt. 0.d0) then
            px = 1.d0 - pcx(l)
         else
            px = pcx(l)
         end if
         zstd = zbrent(-4.d0, 4.d0)
         val = zstd*sigsim + musim
         if (val .le. 0.d0) then
            write (*, *) 'Warning!!'
            write (*, *) 'The percentile is outside the range of calculation'
            write (flog, *) 'Warning!!'
            write (flog, *) 'The percentile is outside the range of calculation'
         else
            val = val**(1.d0/mdpc2)
            write (*, *) '2 m depth-averaged particle concentration'
            write (*, 421) pcx(l), val
            write (fout, *) '2 m depth-averaged particle concentration'
            write (fout, 421) pcx(l), val
         end if
         write (*, *) 'Flow density'
         write (flog, *) 'Flow density'
         musim = muden
         sigsim = sigden
         mm = mden
         if (mm .lt. 0.d0) then
            px = 1.d0 - pcx(l)
         else
            px = pcx(l)
         end if
         zstd = zbrent(-4.d0, 4.d0)
         val = zstd*sigsim + musim
         if (val .le. 0.d0) then
            write (*, *) 'Warning!!'
            write (*, *) 'The percentile is outside the range of calculation'
            write (flog, *) 'Warning!!'
            write (flog, *) 'The percentile is outside the range of calculation'
         else
            val = val**(1.d0/mden)
            write (*, *) 'Flow density'
            write (*, 420) pcx(l), val
            write (fout, *) 'Flow density'
            write (fout, 420) pcx(l), val
         end if
         write (*, *) 'Flow shear velocity'
         write (flog, *) 'Flow shear velocity'
         musim = muush
         sigsim = sigush
         mm = mush
         if (mm .lt. 0.d0) then
            px = 1.d0 - pcx(l)
         else
            px = pcx(l)
         end if
         zstd = zbrent(-4.d0, 4.d0)
         val = zstd*sigsim + musim
         if (val .le. 0.d0) then
            write (*, *) 'Warning!!'
            write (*, *) 'The percentile is outside the range of calculation'
            write (flog, *) 'Warning!!'
            write (flog, *) 'The percentile is outside the range of calculation'
         else
            val = val**(1.d0/mush)
            write (*, *) 'Flow shear velocity'
            write (*, 420) pcx(l), val
            write (fout, *) 'Flow shear velocity'
            write (fout, 420) pcx(l), val
         end if
         if(calc_t_mix) then
            write (*, *) 'Flow temperature'
            write (flog, *) 'Flow temperature'
            musim = mutp
            sigsim = sigtp
            mm = mtp
            if (mm .lt. 0.d0) then
               px = 1.d0 - pcx(l)
            else
               px = pcx(l)
            end if
            zstd = zbrent(-4.d0, 4.d0)
            val = zstd*sigsim + musim
            if (val .le. 0.d0) then
               write (*, *) 'Warning!!'
               write (*, *) 'The percentile is outside the range of calculation'
               write (flog, *) 'Warning!!'
               write (flog, *) 'The percentile is outside the range of calculation'
            else
               val = val**(1.d0/mtp)
               write (*, *) 'Flow temperature'
               write (*, 420) pcx(l), val
               write (fout, *) 'Flow temperature'
               write (fout, 420) pcx(l), val
            end if
            if(usr_z_t) then
               write (*, *) 'Flow temperature at user requested heights'
               write (flog, *) 'Flow temperature at user requested heights'
               do j = 1, 20
                  if(zt(j).eq.undefined) cycle
                  write (*, 424) j, zt(j)
                  write (flog, 424) j, zt(j)
                  musim = mutpz(j)
                  sigsim = sigtpz(j)
                  mm = mtpz(j)
                  if (mm .lt. 0.d0) then
                     px = 1.d0 - pcx(l)
                  else
                     px = pcx(l)
                  end if
                  zstd = zbrent(-4.d0, 4.d0)
                  val = zstd*sigsim + musim
                  if (val .le. 0.d0) then
                     write (*, *) 'Warning!!'
                     write (*, *) 'The percentile is outside the range of calculation'
                     write (flog, *) 'Warning!!'
                     write (flog, *) 'The percentile is outside the range of calculation'
                  else
                     val = val**(1.d0/mtpz(j))
                     write (*, 429) zt(j), pcx(l), val
                     write (fout, 429) zt(j), pcx(l), val
                  end if
               end do
            endif
         endif
         if (usr_z_dynpr) then
            write (*, *) 'Dynamic pressure at user requested heights'
            write (flog, *) 'Dynamic pressure at user requested heights'
            do j = 1, 20
               if(zdynpr(j).eq.undefined) cycle
               write (*, 424) j, zdynpr(j)
               write (flog, 424) j, zdynpr(j)
               musim = mupz(j)
               sigsim = sigpz(j)
               mm = mpz(j)
               if (mm .lt. 0.d0) then
                  px = 1.d0 - pcx(l)
               else
                  px = pcx(l)
               end if
               zstd = zbrent(-4.d0, 4.d0)
               val = zstd*sigsim + musim
               if (val .le. 0.d0) then
                  write (*, *) 'Warning!!'
                  write (*, *) 'The percentile is outside the range of calculation'
                  write (flog, *) 'Warning!!'
                  write (flog, *) 'The percentile is outside the range of calculation'
               else
                  val = val**(1.d0/mpz(j))
                  write (*, 422) zdynpr(j), pcx(l), val
                  write (fout, 422) zdynpr(j), pcx(l), val
               end if
            end do
         endif
         if (usr_z_c) then
            write (*, *) 'Particle concentration at user requested heights'
            write (flog, *) 'Particle concentration at user requested heights'
            do j = 1, 20
               if(zc(j).eq.undefined) cycle
               write (*, 424) j, zc(j)
               write (flog, 424) j, zc(j)
               musim = mucz(j)
               sigsim = sigcz(j)
               mm = mcz(j)
               if (mm .lt. 0.d0) then
                  px = 1.d0 - pcx(l)
               else
                  px = pcx(l)
               end if
               zstd = zbrent(-4.d0, 4.d0)
               val = zstd*sigsim + musim
               if (val .le. 0.d0) then
                  write (*, *) 'Warning!!'
                  write (*, *) 'The percentile is outside the range of calculation'
                  write (flog, *) 'Warning!!'
                  write (flog, *) 'The percentile is outside the range of calculation'
               else
                  val = val**(1.d0/mcz(j))
                  write (*, 423) zc(j), pcx(l), val
                  write (fout, 423) zc(j), pcx(l), val
               end if
            end do
            write (*, *) 'User requested heights depth-averaged particle concentration'
            write (flog, *) 'User requested heights depth-averaged particle concentration'
            do j = 1, 20
               if(zc(j).eq.undefined) cycle
               write (*, 424) j, zc(j)
               write (flog, 424) j, zc(j)
               musim = mudpcz(j)
               sigsim = sigdpcz(j)
               mm = mdpcz(j)
               if (mm .lt. 0.d0) then
                  px = 1.d0 - pcx(l)
               else
                  px = pcx(l)
               end if
               zstd = zbrent(-4.d0, 4.d0)
               val = zstd*sigsim + musim
               if (val .le. 0.d0) then
                  write (*, *) 'Warning!!'
                  write (*, *) 'The percentile is outside the range of calculation'
                  write (flog, *) 'Warning!!'
                  write (flog, *) 'The percentile is outside the range of calculation'
               else
                  val = val**(1.d0/mdpcz(j))
                  write (*, 430) zc(j), pcx(l), val
                  write (fout, 430) zc(j), pcx(l), val
               end if
            end do
         endif
404      if (.not. deprates) goto 405
         write (*, *) 'Total deposition rate'
         write (flog, *) 'Total deposition rate'
         musim = murtot
         sigsim = sigrtot
         mm = mrtot
         if (mm .lt. 0.d0) then
            px = 1.d0 - pcx(l)
         else
            px = pcx(l)
         end if
         zstd = zbrent(-4.d0, 4.d0)
         val = zstd*sigsim + musim
         if (val .le. 0.d0) then
            write (*, *) 'Warning!!'
            write (*, *) 'The percentile is outside the range of calculation'
            write (flog, *) 'Warning!!'
            write (flog, *) 'The percentile is outside the range of calculation'
         else
            val = val**(1.d0/mrtot)
            write (*, *) 'Total deposition rate'
            write (*, 421) pcx(l), val
            write (fout, *) 'Total deposition rate'
            write (fout, 421) pcx(l), val
         end if
         write (*, *) 'Total deposition time'
         write (flog, *) 'Total deposition time'
         musim = mutdep
         sigsim = sigtdep
         mm = mtdep
         if (mm .lt. 0.d0) then
            px = 1.d0 - pcx(l)
         else
            px = pcx(l)
         end if
         zstd = zbrent(-4.d0, 4.d0)
         val = zstd*sigsim + musim
         if (val .le. 0.d0) then
            write (*, *) 'Warning!!'
            write (*, *) 'The percentile is outside the range of calculation'
            write (flog, *) 'Warning!!'
            write (flog, *) 'The percentile is outside the range of calculation'
            return
         else
            val = val**(1.d0/mtdep)
            write (*, *) 'Total deposition time'
            write (*, 421) pcx(l), val
            write (fout, *) 'Total deposition time'
            write (fout, 421) pcx(l), val
         end if
405      l = l + 1
      end do

369   format('Temp. simmetrization exponent', f8.3,/)
410   format('Symmetrization exponent', f8.3, /,&
       &'Median', e12.5, /,&
       &'Standard deviation', e12.5, //)
411   format(f6.2, 1x, 'm dynamic pressure probability function', /,&
       &'Symmetrization exponent', f8.3, /,&
       &'Median', e12.5, /,&
       &'Standard deviation', e12.5, //)
412   format(f6.2, 1x, 'm particle concentration probability function', /,&
       &'Symmetrization exponent', f8.3, /,&
       &'Median', e12.5, /,&
       &'Standard deviation', e12.5, //)
413   format(f6.2, 1x, 'm flow temperature probability function', /,&
       &'Symmetrization exponent', f8.3, /,&
       &'Median', e12.5, /,&
       &'Standard deviation', e12.5, //)          
420   format('Percentile', f12.3, /,&
       &'Function value', f14.4,/)
421   format('Percentile', f12.3, /,&
       &'Function value', e14.4,/)
422   format(f6.2, 1x, 'm', 2x, 'average dynamic pressure (Pa)', /,&
       &'Percentile', f12.3, /,&
       &'Function value', f14.4,/)
423   format(f6.2, 1x, 'm', 2x, 'particle concentration', /,&
       &'Percentile', f12.3, /,&
       &'Function value', e14.4,/)
430   format(f6.2, 1x, 'm', 2x, 'depth-averaged particle concentration', /,&
       &'Percentile', f12.3, /,&
       &'Function value', e14.4,/)
424   format(i2, 2x, 'z = ', f6.2)
425   format(f6.2, 1x, 'Pdyn simmetrization exponent calc. residuals')
426   format(f6.2, 1x, 'C simmetrization exponent calc. residuals')
427   format('Percentile = ', f6.3,/)
428   format(/, 'FATAL ERROR! PCX('i2') > 1 or < 0',/)
429   format(f6.2, 1x, 'm', 2x, 'average flow temperature (K)', /,&
       &'Percentile', f12.3, /,&
       &'Function value', f14.4,/)
   end subroutine probfunction
