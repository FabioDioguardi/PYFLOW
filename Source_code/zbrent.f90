	FUNCTION zbrent(x1,x2)
        USE inputdata; USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x1,x2
	REAL(dp) :: zbrent
	INTERFACE
		FUNCTION func(x)
		USE inputdata; USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(dp), PARAMETER :: EPS=epsilon(x1)
	INTEGER(I4B) :: iter
	REAL(dp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol,xm
	logical :: checkzbr
      common/warn/ checkzbr
      checkzbr=.false.
	tol=xacc
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)
	if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
		call nrerror('root must be bracketed for zbrent')
	c=b
	fc=fb
	do iter=1,ITMAX
		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
			fc=fa
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol=2.0_dp*EPS*abs(b)+0.5_dp*tol
		xm=0.5_dp*(c-b)
               write(52,354)xm,fb
               write(*,354)xm,fb
  354 format(2(e10.3,2x))
		if (abs(xm) <= tol .or. fb == 0.0) then
			zbrent=b
			RETURN
		end if
		if (abs(e) >= tol .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_dp*xm*s
				q=1.0_dp-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
				q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_dp*p  <  min(3.0_dp*xm*q-abs(tol*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol,xm), abs(d) > tol )
		fb=func(b)
	end do
	call nrerror('zbrent: exceeded maximum iterations')
	zbrent=b
	END FUNCTION zbrent
