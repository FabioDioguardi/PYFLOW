      subroutine testchi(nmax, p, w, p50, sr)
         USE inoutdata; USE nrtype
         implicit none
         INTERFACE

            FUNCTION qsimp(a, b)
               USE nrtype; USE nrutil, ONLY: nrerror
               REAL(dp), INTENT(IN) :: a, b
               REAL(dp) :: qsimp
            END FUNCTION qsimp

         END INTERFACE
         integer :: nmax, i, j, l, mmax
         real(dp), dimension(nmax) :: p, w, z
         real(dp), dimension(20) :: e, o, ztest
         real(dp) :: p50, sr, chical, chitab, x1, x2, rtflsp, sens
         sens = sensgrchi
!     Standardization of the grainsize classes
44       do i = 1, nmax
            z(i) = (p(i) - p50)/sr
            o(i) = w(i)
         end do
!     Put together classes with wt<5
         i = 1
         j = 1
!      mmax=i    !old
         mmax = nmax !new
         o(i) = w(i)
         ztest(i) = z(i)
43       if (o(i) .lt. sens) then
            j = j + 1
            mmax = mmax - 1 !new
            if (j .gt. nmax) goto 42
            o(i) = o(i) + w(j)
            ztest(i) = z(j)
            goto 43
         else
            i = i + 1
!      mmax=i  !old
            if (i .gt. nmax) goto 42
            j = j + 1
            o(i) = w(j)
            ztest(i) = z(j)
            goto 43
         end if
42       if (o(mmax) .lt. sens) then
            mmax = mmax - 1
            o(mmax) = o(mmax) + o(mmax + 1)
         else
         end if
         ztest(mmax) = 5.d0
         nfunc1 = 3
!     Calculation of the expected values
         e(1) = qsimp(-5.d0, ztest(1))
         do l = 2, mmax
            e(l) = qsimp(ztest(l - 1), ztest(l))
         end do
!     Test
         chical = 0.d0
         do l = 1, mmax
            o(l) = o(l)*100.d0
            e(l) = e(l)*100.d0
            chical = chical + ((o(l) - (e(l)))**2)/e(l)
         end do
         nfreedchi = mmax - 3
         if (nfreedchi .le. 0.d0) then
            sens = sens - 0.01d0
            goto 44
         else
!      xmaxchi=300.d0
!      dxint=1.d0
            x1 = 200.d0
            x2 = 1.d-5
            nfunc = 20
            nfunc1 = 4
            write (flog, *) 'Chi tabulated calculation residuals'
            write (*, *) 'Chi tabulated calculation residuals'
            chitab = rtflsp(x1, x2)
            write (57, *) ''
            write (57, *) '********CHI SQUARED TEST************'
            write (57, *) 'ztest   Oj    Ej   (Oj-Ej)^2/Ej'
            write (*, *) ''
            write (*, *) '********CHI SQUARED TEST************'
            write (*, *) 'ztest   Oj    Ej   (Oj-Ej)^2/Ej'
            do i = 1, mmax
               write (57, 102) ztest(i), o(i), e(i), ((o(i) - (e(i)))**2)/e(i)
               write (*, 102) ztest(i), o(i), e(i), ((o(i) - (e(i)))**2)/e(i)
            end do
102         format(4(f5.2, 2x))
            write (*, 107) sens
            write (*, 103) probchi
            write (*, 104) nfreedchi
            write (*, 105) chical
            write (*, 106) chitab
            write (57, 107) sens
            write (57, 103) probchi
            write (57, 104) nfreedchi
            write (57, 105) chical
            write (57, 106) chitab
            if (chical .le. chitab) then
               write (57, *) 'The distribution is not significantly different from a normal standard distribution'
               write (*, *) 'The distribution is not significantly different from a normal standard distribution'
            else
               write (57, *) 'The distribution is significantly different from a normal standard distribution'
               write (*, *) 'The distribution is significantly different from a normal standard distribution'
            end if
103         format('Significance =        ', f4.2)
104         format('Degrees of freedom =   ', f4.1)
105         format('Chi calculated =    ', f6.2)
106         format('Chi tabulated =     ', f6.2,/)
107         format(/, 'Sensitivity =       ', f6.2)
         end if
         nfunc1 = 0
         return
      end
