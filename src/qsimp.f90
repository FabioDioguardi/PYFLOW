        FUNCTION qsimp(a,b)
        USE inoutdata; USE nrtype; USE nrutil, ONLY : nrerror
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: a,b
        REAL(DP) :: qsimp
        INTERFACE
                 FUNCTION funcq(x)
                 USE nrtype
                 REAL(DP), DIMENSION(:), INTENT(IN) :: x
                 REAL(DP), DIMENSION(size(x)) :: funcq
                 END FUNCTION funcq
        END INTERFACE
        INTEGER(I4B), PARAMETER :: JMAX=40
        REAL(DP), PARAMETER :: EPS=1.0e-6_DP
        INTEGER(I4B) :: j
        REAL(DP) :: os,ost,st
        ost=0.0
        os= 0.0
        do j=1,JMAX
        call trapzd(funcq,a,b,st,j)
        qsimp=(4.0_DP*st-ost)/3.0_DP
!        write(52,*)a,b,st,j,qsimp,abs(qsimp-os),EPS*abs(os)
         if(isnan(qsimp).or.abs(qsimp).gt.1.d100) then
                 if(nnewt.eq.2) then
                 write(52,*)qsimp
                 write(52,*)'ERROR. Not able to find a solution for the reference level'
                 write(52,*)'Try changing Z0AVGGUESS and/or Z0MAXGUESS and/or Z0MINGUESS'
                 write(*,*)'ERROR. Not able to find a solution for the reference level'
                 write(*,*)'Try changing Z0AVGGUESS and/or Z0MAXGUESS and/or Z0MINGUESS'
                 stop
                 else
                 endif
         endif
         if (j > 5) then
         if (abs(qsimp-os) < EPS*abs(os) .or. &
         (qsimp == 0.0 .and. os == 0.0)) RETURN
          end if
         os=qsimp
         ost=st
         end do
!         write(52,*)qsimp
        call nrerror('qsimp: too many steps')
        END FUNCTION qsimp

	SUBROUTINE trapzd(funcq,a,b,s,n)
	USE inoutdata; USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION funcq(x)
		USE inoutdata; USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: funcq
		END FUNCTION funcq
	END INTERFACE
	REAL(dp) :: del,fsum
	INTEGER(I4B) :: it
	if (n == 1) then
		s=0.5_dp*(b-a)*sum(funcq( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(funcq(arth(a+0.5_dp*del,del,it)))
		s=0.5_dp*(s+del*fsum)
	end if
!	if(nfunc.eq.10) write(*,*)'trapzd',a,b,s
	END SUBROUTINE trapzd

	FUNCTION funcq(x)
	USE inoutdata; USE nrtype
	implicit none
	INTERFACE
		FUNCTION func(y)
		USE inoutdata; USE nrtype
		REAL(dp), INTENT(IN) :: y
		REAL(dp) :: func
		END FUNCTION func

		FUNCTION func1(y)
		USE inoutdata; USE nrtype
		REAL(dp), INTENT(IN) :: y
		REAL(dp) :: func1
		END FUNCTION func1
	END INTERFACE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(size(x)) :: funcq
	integer :: i,n
	n=size(x)
      if(nfunc.eq.7.or.nfunc.eq.8.or.nfunc1.eq.3.or.nfunc1.eq.4) then
	do i=1,n
	funcq(i)=func1(x(i))
	enddo
      else
	do i=1,n
	funcq(i)=func(x(i))
	enddo
      end if
	end function funcq
