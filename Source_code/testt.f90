      subroutine testt(d,dmod,sigma,nrest)
      USE inputdata; USE nrtype; USE nrutil
      implicit none
      INTERFACE
		FUNCTION zbrent(x1,x2)
		USE inputdata; USE nrtype; USE nrutil, ONLY : nrerror
                IMPLICIT NONE
	        REAL(dp), INTENT(IN) :: x1,x2
         	REAL(dp) :: zbrent
         	END FUNCTION zbrent

		FUNCTION func(x)
		USE inputdata; USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp) :: func
		END FUNCTION func
	END INTERFACE
      real(dp) :: d,dmod,sigma,tcalc,ttab,x1,x2,xmax,nfreed,nclass
      integer :: nfunc,nrest
      common/nf/ nfunc
      common/tres/ tcalc,ttab
      common/tinp/ nfreed,nclass
      nfreed=nclass-1.d0
      tcalc=(d-dmod)/(sigma*sqrt(1.d0/nfreed))
      tcalc=abs(tcalc)
      write(52,*)''
      write(52,*)'T tabulated calculation residuals'
      write(*,*)''
      write(*,*)'T tabulated calculation residuals'
      xmax=800.d0
      x1=1.d-5
      x2=800.d0
      nfunc=4
      ttab=zbrent(x1,x2)
      nrest=1
      if(tcalc.gt.ttab) nrest=2
      end subroutine testt

      function tstud(x)
      USE nrtype; USE nrutil
      implicit none
      INTERFACE
		FUNCTION betai_s(a,b,x)
		USE nrtype; USE nrutil, ONLY : assert
!        	USE nr, ONLY : betacf,gammln
                IMPLICIT NONE
        	REAL(dp), INTENT(IN) :: a,b,x
        	REAL(dp) :: betai_s
         	END FUNCTION betai_s
	END INTERFACE
      real(dp) :: a,b,x,nfreed,nclass,tstud
      common/tinp/ nfreed,nclass
      nfreed=nclass-1.d0
      tstud=betai_s(0.5d0*nfreed,0.5d0,nfreed/(nfreed+x**2))
      end function tstud
      
	FUNCTION betai_s(a,b,x)
	USE nrtype; USE nrutil, ONLY : assert
!	USE nr, ONLY : betacf,gammln
	IMPLICIT NONE
          INTERFACE
		FUNCTION betacf(a,b,x)
        	USE nrtype; USE nrutil, ONLY : nrerror
        	IMPLICIT NONE
        	REAL(dp), INTENT(IN) :: a,b,x
         	REAL(dp) :: betacf
         	END FUNCTION betacf
         	
        	FUNCTION gammln(xx)
        	USE nrtype; USE nrutil, ONLY : arth,assert
        	IMPLICIT NONE
          	REAL(dp), INTENT(IN) :: xx
          	REAL(dp) :: gammln
          	END FUNCTION gammln
         	
	END INTERFACE
	REAL(dp), INTENT(IN) :: a,b,x
	REAL(dp) :: betai_s
	REAL(dp) :: bt
	call assert(x >= 0.0, x <= 1.0, 'betai_s arg')
	if (x == 0.0 .or. x == 1.0) then
		bt=0.0
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)&
			+a*log(x)+b*log(1.0_dp-x))
	end if
	if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
		betai_s=bt*betacf(a,b,x)/a
	else
		betai_s=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
	end if
	END FUNCTION betai_s

	FUNCTION betacf(a,b,x)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,b,x
	REAL(dp) :: betacf
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(dp), PARAMETER :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
	REAL(dp) :: aa,c,d,del,h,qab,qam,qap
	INTEGER(I4B) :: m,m2
	qab=a+b
	qap=a+1.0_dp
	qam=a-1.0_dp
	c=1.0
	d=1.0_dp-qab*x/qap
	if (abs(d) < FPMIN) d=FPMIN
	d=1.0_dp/d
	h=d
	do m=1,MAXIT
		m2=2*m
		aa=m*(b-m)*x/((qam+m2)*(a+m2))
		d=1.0_dp+aa*d
		if (abs(d) < FPMIN) d=FPMIN
		c=1.0_dp+aa/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		h=h*d*c
		aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
		d=1.0_dp+aa*d
		if (abs(d) < FPMIN) d=FPMIN
		c=1.0_dp+aa/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		del=d*c
		h=h*del
		if (abs(del-1.0_dp) <= EPS) exit
	end do
	if (m > MAXIT)&
		call nrerror('a or b too big, or MAXIT too small in betacf_s')
	betacf=h
	END FUNCTION betacf

	FUNCTION gammln(xx)
	USE nrtype; USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: xx
	REAL(dp) :: gammln
	REAL(DP) :: tmp,x
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	call assert(xx > 0.0, 'gammln_s arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	gammln=tmp+log(stp*(1.000000000190015_dp+&
		sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
	END FUNCTION gammln

