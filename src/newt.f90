	SUBROUTINE newt(x,check,xsize)
	USE nrtype; USE nrutil, ONLY : nrerror,vabs
!	USE nr, ONLY : fdjac,lnsrch,lubksb,ludcmp
	USE fminln
	IMPLICIT NONE
	
        INTERFACE
        
                 SUBROUTINE fdjac(x,fvec,df)
	         USE nrtype; USE nrutil, ONLY : assert_eq
	         IMPLICIT NONE
	         REAL(dp), DIMENSION(:), INTENT(IN) :: fvec
	         REAL(dp), DIMENSION(:), INTENT(INOUT) :: x
	         REAL(dp), DIMENSION(:,:), INTENT(OUT) :: df
	         END SUBROUTINE fdjac
	         
	         SUBROUTINE ludcmp(a,indx,d)
	         USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
	         IMPLICIT NONE
	         REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: a
	         INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	         REAL(dp), INTENT(OUT) :: d
	         END SUBROUTINE ludcmp

	         SUBROUTINE lubksb(a,indx,b)
	         USE nrtype; USE nrutil, ONLY : assert_eq
	         IMPLICIT NONE
	         REAL(dp), DIMENSION(:,:), INTENT(IN) :: a
	         INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
          	 REAL(dp), DIMENSION(:), INTENT(INOUT) :: b
           	 END SUBROUTINE lubksb
        
	         SUBROUTINE lnsrch(xold,fold,gg,p,x,f,stpmax,check,funcln)
	         USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,vabs
	         IMPLICIT NONE
	         REAL(dp), DIMENSION(:), INTENT(IN) :: xold,gg
	         REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
	         REAL(dp), INTENT(IN) :: fold,stpmax
	         REAL(dp), DIMENSION(:), INTENT(OUT) :: x
	         REAL(dp), INTENT(OUT) :: f
	         LOGICAL(LGT), INTENT(OUT) :: check
          	 INTERFACE
		 FUNCTION funcln(x)
                 USE nrtype
		 IMPLICIT NONE
		 REAL(dp) :: funcln
		 REAL(dp), DIMENSION(:), INTENT(IN) :: x
	 	 END FUNCTION funcln
  	         END INTERFACE
	         END SUBROUTINE lnsrch
	         
                function funcv(x)
                USE inoutdata; USE nrtype
                implicit none
                      INTERFACE
                      
		                FUNCTION qsimp(a,b)
        	                USE nrtype; USE nrutil, ONLY : nrerror
         	                REAL(dp), INTENT(IN) :: a,b
        	                REAL(dp) :: qsimp
		                END FUNCTION qsimp
                       		
		                FUNCTION cdmodel(dpart,rhop,rho)
                                USE inoutdata; USE nrtype
         	                REAL(dp), INTENT(IN) :: dpart,rhop,rho
                         	REAL(dp) :: cdmodel
	                 	END FUNCTION cdmodel
	                 	
                      END INTERFACE
                REAL(dp), DIMENSION(:), INTENT(IN) :: x
                REAL(dp), DIMENSION(size(x)) :: funcv
                end function funcv
        
        END INTERFACE
        
      	REAL(dp), DIMENSION(xsize), INTENT(INOUT) :: x
	LOGICAL(LGT), INTENT(OUT) :: check
	INTEGER(I4B), PARAMETER :: MAXITS=10000000
	REAL(dp), PARAMETER :: TOLF=1.0e-4_dp,TOLMIN=1.0e-6_dp,TOLX=epsilon(x),&
		STPMX=100.0
	INTEGER(I4B) :: its
	INTEGER(I4B), DIMENSION(size(x)) :: indx
	REAL(dp) :: d,f,fold,stpmax
	REAL(dp), DIMENSION(size(x)) :: g,p,xold
	REAL(dp), DIMENSION(size(x)), TARGET :: fvec
	REAL(dp), DIMENSION(size(x),size(x)) :: fjac
	INTEGER :: xsize
	REAL(dp) :: res1,res2
	fmin_fvecp=>fvec
	f=fmin(x)
	res1=maxval(abs(fvec(:)))
	write(*,*)'res_func_vals',res1
	write(52,*)'res_func_vals',res1
	if (maxval(abs(fvec(:))) < 0.01_dp*TOLF) then
		check=.false.
		RETURN
	end if
	stpmax=STPMX*max(vabs(x(:)),real(size(x),dp))
	do its=1,MAXITS
		call fdjac(x,fvec,fjac)
!		write(*,*)'fjac',fjac
		g(:)=matmul(fvec(:),fjac(:,:))
!		write(*,*),g
		xold(:)=x(:)
		fold=f
		p(:)=-fvec(:)
		call ludcmp(fjac,indx,d)
		call lubksb(fjac,indx,p)
		call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin)
		res1=maxval(abs(fvec(:)))
        	write(*,*)'res_func_vals',res1
        	write(52,*)'res_func_vals',res1
 		if (maxval(abs(fvec(:))) < TOLF) then
			check=.false.
			RETURN
		end if
		if (check) then
			check=(maxval(abs(g(:))*max(abs(x(:)),1.0_dp) / &
				max(f,0.5_dp*size(x))) < TOLMIN)
			RETURN
		end if
		res2=maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0_dp))
		write(*,*)'res_dx       ',res2
		write(52,*)'res_dx       ',res2
		if (maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0_dp)) < TOLX) &
			RETURN
	end do
	call nrerror('MAXITS exceeded in newt')
	END SUBROUTINE newt

	SUBROUTINE lnsrch(xold,fold,gg,p,x,f,stpmax,check,funcln)
	USE inoutdata; USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,vabs
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: xold,gg
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
	REAL(dp), INTENT(IN) :: fold,stpmax
	REAL(dp), DIMENSION(:), INTENT(OUT) :: x
	REAL(dp), INTENT(OUT) :: f
	LOGICAL(LGT), INTENT(OUT) :: check
	INTERFACE
		FUNCTION funcln(x)
                USE nrtype
		IMPLICIT NONE
		REAL(dp) :: funcln
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION funcln
	END INTERFACE
	REAL(dp), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x)
	INTEGER(I4B) :: ndum
	REAL(dp) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
	ndum=assert_eq(size(gg),size(p),size(x),size(xold),'lnsrch')
	check=.false.
	pabs=vabs(p(:))
	if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
	slope=dot_product(gg,p)
	if (slope >= 0.0) call nrerror('roundoff problem in lnsrch')
	alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp))
	alam=1.0
	do
                x(:)=xold(:)+alam*p(:)
               if(nnewt.ne.3.and.x(2).lt.0.d0) x(2)=abs(x(2))
               if(nnewt.ne.3.and.x(2).gt.ztot) x(2)=ztot-1.d-5
		f=funcln(x)
		if (alam < alamin) then
			x(:)=xold(:)
			check=.true.
			RETURN
		else if (f <= fold+ALF*alam*slope) then
			RETURN
		else
			if (alam == 1.0) then
				tmplam=-slope/(2.0_dp*(f-fold-slope))
			else
				rhs1=f-fold-alam*slope
				rhs2=f2-fold-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
					(alam-alam2)
				if (a == 0.0) then
					tmplam=-slope/(2.0_dp*b)
				else
					disc=b*b-3.0_dp*a*slope
					if (disc < 0.0) then
						tmplam=0.5_dp*alam
					else if (b <= 0.0) then
						tmplam=(-b+sqrt(disc))/(3.0_dp*a)
					else
						tmplam=-slope/(b+sqrt(disc))
					end if
				end if
				if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
			end if
		end if
		alam2=alam
		f2=f
		alam=max(tmplam,0.1_dp*alam)
	end do
	END SUBROUTINE lnsrch


	SUBROUTINE fdjac(x,fvec,df)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: fvec
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: x
	REAL(dp), DIMENSION(:,:), INTENT(OUT) :: df
	INTERFACE
		FUNCTION funcv(x)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: x
		REAL(dp), DIMENSION(size(x)) :: funcv
		END FUNCTION funcv
	END INTERFACE
	REAL(dp), PARAMETER :: EPS=1.0e-4_dp
	INTEGER(I4B) :: j,n
	REAL(dp), DIMENSION(size(x)) :: xsav,xph,h
	n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
	xsav=x
	h=EPS*abs(xsav)
	where (h == 0.0) h=EPS
	xph=xsav+h
	h=xph-xsav
	do j=1,n
		x(j)=xph(j)
		df(:,j)=(funcv(x)-fvec(:))/h(j)
		x(j)=xsav(j)
	end do
	END SUBROUTINE fdjac


	SUBROUTINE ludcmp(a,indx,d)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
	IMPLICIT NONE
	REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(dp), INTENT(OUT) :: d
	REAL(dp), DIMENSION(size(a,1)) :: vv
	REAL(dp), PARAMETER :: TINY=1.0e-20_dp
	INTEGER(I4B) :: j,n,imax
	n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
	d=1.0
	vv=maxval(abs(a),dim=2)
	if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
	vv=1.0_dp/vv
	do j=1,n
		imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
		if (j /= imax) then
			call swap(a(imax,:),a(j,:))
			d=-d
			vv(imax)=vv(j)
		end if
		indx(j)=imax
		if (a(j,j) == 0.0) a(j,j)=TINY
		a(j+1:n,j)=a(j+1:n,j)/a(j,j)
		a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
	end do
	END SUBROUTINE ludcmp

	SUBROUTINE lubksb(a,indx,b)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:,:), INTENT(IN) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(dp), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,n,ii,ll
	REAL(dp) :: summ
	n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
	ii=0
	do i=1,n
		ll=indx(i)
		summ=b(ll)
		b(ll)=b(i)
		if (ii /= 0) then
			summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
		else if (summ /= 0.0) then
			ii=i
		end if
		b(i)=summ
	end do
	do i=n,1,-1
		b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
	end do
	END SUBROUTINE lubksb
