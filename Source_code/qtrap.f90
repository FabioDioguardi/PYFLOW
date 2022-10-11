	FUNCTION qtrap(a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
!	USE nr, ONLY : trapzd
	IMPLICIT NONE
	INTERFACE
		FUNCTION funcq(x)
		USE nrtype
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: funcq
		END FUNCTION funcq
	END INTERFACE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp) :: qtrap
	INTEGER(I4B), PARAMETER :: JMAX=40
	REAL(dp), PARAMETER :: EPS=1.0e-6_dp
	REAL(dp) :: olds
	INTEGER(I4B) :: j
	olds=0.0
	do j=1,JMAX
		call trapzd(funcq,a,b,qtrap,j)
		      if (isnan(qtrap)) return
		if (j > 5) then
			if (abs(qtrap-olds) < EPS*abs(olds) .or. &
				(qtrap == 0.0 .and. olds == 0.0)) RETURN
		end if
		olds=qtrap
	end do
	call nrerror('qtrap: too many steps')
	END FUNCTION qtrap

	SUBROUTINE trapzd(funcq,a,b,s,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b
	REAL(dp), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION funcq(x)
		USE nrtype
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
	END SUBROUTINE trapzd

	FUNCTION funcq(x)
	USE nrtype
	implicit none
	INTERFACE
		FUNCTION func(y)
		USE nrtype
		REAL(dp), INTENT(IN) :: y
		REAL(dp) :: func
		END FUNCTION func

		FUNCTION func1(y)
		USE nrtype
		REAL(dp), INTENT(IN) :: y
		REAL(dp) :: func1
		END FUNCTION func1
	END INTERFACE
	REAL(dp), DIMENSION(:), INTENT(IN) :: x
	REAL(dp), DIMENSION(size(x)) :: funcq
	integer :: i,n,nfunc
        common/nf/ nfunc
	n=size(x)
      if(nfunc.eq.7.or.nfunc.eq.8) then
	do i=1,n
	funcq(i)=func1(x(i))
	enddo
      else
	do i=1,n
	funcq(i)=func(x(i))
	enddo
      end if
	end function funcq
