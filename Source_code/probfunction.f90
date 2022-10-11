      subroutine probfunction
      USE inputdata; USE nrtype
      implicit none
	INTERFACE
		FUNCTION zbrent(x1,x2)
                USE inputdata; USE nrtype; USE nrutil, ONLY : nrerror
                IMPLICIT NONE
	        REAL(dp) :: x1,x2
         	REAL(dp) :: zbrent
         	END FUNCTION zbrent
	END INTERFACE
      real(dp), dimension(20) :: zdynpr,pzavg,pzmax,pzmin,czavg,czmax,czmin,zc,mpz,mupz,&
     &sigpz,mcz,mucz,sigcz
      real(dp) :: den,p10avg,p10max,p10min,c2avg,c2max,c2min,px,pcx,sigsim,mudstr,mxdstr,mndstr,&
     &mm,musim,temp1,temp2,val,tempz1,tempz2,zstd,x1,x2,mp10,mc2,mup10,muc2,sigc2,sigp10
      real(dp) :: x1tmp,x2tmp,x3tmp,x4tmp
      integer :: j,ipr,ic,ndynpr,nc,nnewt,njc,njpr,npfunc,nsrch,nfunc
      logical :: checkzbr
      common/warn/ checkzbr
      common/nf/ nfunc
      common/cnewt2/ den,nnewt
      common/prof3/ zdynpr,p10avg,p10max,p10min,c2avg,c2max,c2min,pzavg,pzmax,pzmin,&
     &czavg,czmax,czmin,zc,ipr,ic,ndynpr,nc
      common/prob/ px,sigsim,mudstr,mxdstr,mndstr,mm,musim
      x1tmp=-50.d0;x2tmp=-1.d-10;x3tmp=1.d-10;x4tmp=50.d0
!     Symmetric probability distribution parameters for Pdyn at 10 m
      mudstr=p10avg
      mxdstr=p10max
      mndstr=p10min
      nfunc=18
      write(*,*)'Pdyn 10 m simmetrization exponent calc. residuals'
      write(52,*)'Pdyn 10 m simmetrization exponent calc. residuals'
  600 temp1=zbrent(x3tmp,x4tmp)
      if(checkzbr) temp1=0.d0
      write(52,369)temp1
      write(*,369)temp1
  369 format('Temp. simmetrization exponent',f8.3,/)
      if(temp1.ne.0.d0.and.temp1.ge.1.d-9.and.temp1.ne.50.d0) then
      mp10=temp1
      else
      temp2=zbrent(x1tmp,x2tmp)
      if(checkzbr) then
      x1tmp=x1tmp-50.d0
      x2tmp=x2tmp-50.d0
      x3tmp=x3tmp+50.d0
      x4tmp=x4tmp+50.d0
      goto 600
      endif
      write(52,369)temp2
      write(*,369)temp2
      if(temp2.ne.0.d0.and.temp2.le.-1.d-9.and.temp2.ne.-50.d0) then
      mp10=temp2
      else
      write(*,*)'Warning. Unable to find a simmetrization coefficient'
      write(*,*)'for 10 m dynamic pressure prob. function'
      endif
      endif
      musim=mudstr**mp10
      sigsim=abs(mxdstr**mp10-mudstr**mp10)
      mup10=musim
      sigp10=sigsim
      write(50,*)'10 m dynamic pressure probability function'
      write(50,410)mp10,mup10,sigp10
      write(*,*)'10 m dynamic pressure probability function'
      write(*,410)mp10,mup10,sigp10
!     Symmetric probability distribution parameters for C at 2 m
      mudstr=c2avg
      mxdstr=c2max
      mndstr=c2min
      write(*,*)'C 2 m simmetrization exponent calc. residuals'
      write(52,*)'C 2 m simmetrization exponent calc. residuals'
  601 temp1=zbrent(x3tmp,x4tmp)
      if(checkzbr) temp1=0.d0
      write(52,369)temp1
      write(*,369)temp1
      if(temp1.ne.0.d0.and.temp1.ge.1.d-9.and.temp1.ne.50.d0) then
      mc2=temp1
      else
      temp2=zbrent(x1tmp,x2tmp)
      if(checkzbr) then
      x1tmp=x1tmp-50.d0
      x2tmp=x2tmp-50.d0
      x3tmp=x3tmp+50.d0
      x4tmp=x4tmp+50.d0
      goto 601
      endif
      write(52,369)temp2
      write(*,369)temp2
      if(temp2.ne.0.d0.and.temp2.le.-1.d-9.and.temp2.ne.-50.d0) then
      mc2=temp2
      else
      write(*,*)'Warning. Unable to find a simmetrization coefficient'
      write(*,*)'for 2 m particle concentration prob. function'
      endif
      endif
      musim=mudstr**mc2
      sigsim=abs(mxdstr**mc2-mudstr**mc2)
      muc2=musim
      sigc2=sigsim
      write(50,*)'2 m particle concentration probability function'
      write(50,410)mc2,muc2,sigc2
      write(*,*)'2 m particle concentration probability function'
      write(*,410)mc2,muc2,sigc2
  410 format('Symmetrization exponent',f8.3,/,&
     &'Median',e12.4,/,&
     &'Standard deviation',e12.4,//)
      if(ndynpr.eq.2.and.nc.eq.2) goto 405
      if(ndynpr.eq.1) then
!     Determination of probability function of Pdyn at user requested heights
      write(*,*)'###PROBABILITY FUNCTIONS FOR AVERAGE DYNAMIC PRESSURE OVER USER REQUESTED HEIGHTS###'
      do j=1,ipr
      mudstr=pzavg(j)
      mxdstr=pzmax(j)
      mndstr=pzmin(j)
      write(*,425)zdynpr(j)
      write(52,425)zdynpr(j)
  602 tempz1=zbrent(x3tmp,x4tmp)
      if(checkzbr) tempz1=0.d0
      write(52,369)tempz1
      write(*,369)tempz1
      if(tempz1.ne.0.d0.and.tempz1.ge.1.d-9.and.tempz1.ne.50.d0) then
      mpz(j)=tempz1
      else
      tempz2=zbrent(x1tmp,x2tmp)
      if(checkzbr) then
      x1tmp=x1tmp-50.d0
      x2tmp=x2tmp-50.d0
      x3tmp=x3tmp+50.d0
      x4tmp=x4tmp+50.d0
      goto 602
      endif
      write(52,369)tempz2
      write(*,369)tempz2
      if(tempz2.ne.0.d0.and.tempz2.le.-1.d-9.and.tempz2.ne.-50.d0) then
      mpz(j)=tempz2
      else
      write(*,*)'Warning. Unable to find a simmetrization coefficient'
      write(*,*)'for dynamic pressure prob. function at z=',zdynpr(j)
      write(52,*)'Warning. Unable to find a simmetrization coefficient'
      write(52,*)'for dynamic pressure prob. function at z=',zdynpr(j)
      endif
      endif
      musim=mudstr**mpz(j)
      sigsim=abs(mxdstr**mpz(j)-mudstr**mpz(j))
      mupz(j)=musim
      sigpz(j)=sigsim
      write(50,411)zdynpr(j),mpz(j),mupz(j),sigpz(j)
      write(*,411)zdynpr(j),mpz(j),mupz(j),sigpz(j)
  411 format(f6.2,1x,'dynamic pressure probability function',/,&
     &'Symmetrization exponent',f8.3,/,&
     &'Median',e12.4,/,&
     &'Standard deviation',e12.4,//)
      enddo
      endif
      if(nc.eq.1) then
!     Determination of probability function of C at user requested heights
      write(*,*)'###PROBABILITY FUNCTIONS FOR PARTICLE CONCENTRATION AT USER REQUESTED HEIGHTS###'
      do j=1,ic
      mudstr=czavg(j)
      mxdstr=czmax(j)
      mndstr=czmin(j)
      write(*,426)zdynpr(j)
      write(52,426)zdynpr(j)
  603 tempz1=zbrent(x3tmp,x4tmp)
      if(checkzbr) tempz1=0.d0
      write(52,369)tempz1
      write(*,369)tempz1
      if(tempz1.ne.0.d0.and.tempz1.ge.1.d-9.and.tempz1.ne.50.d0) then
      mcz(j)=tempz1
      else
      tempz2=zbrent(x1tmp,x2tmp)
      if(checkzbr) then
      x1tmp=x1tmp-50.d0
      x2tmp=x2tmp-50.d0
      x3tmp=x3tmp+50.d0
      x4tmp=x4tmp+50.d0
      goto 603
      endif
      write(52,369)tempz2
      write(*,369)tempz2
      if(tempz2.ne.0.d0.and.tempz2.le.-1.d-9.and.tempz2.ne.-50.d0) then
      mcz(j)=tempz2
      else
      write(*,*)'Warning. Unable to find a simmetrization coefficient'
      write(*,*)'for dynamic pressure prob. function at z=',zc(j)
      write(52,*)'Warning. Unable to find a simmetrization coefficient'
      write(52,*)'for dynamic pressure prob. function at z=',zc(j)
      endif
      endif
      musim=mudstr**mcz(j)
      sigsim=abs(mxdstr**mcz(j)-mudstr**mcz(j))
      mucz(j)=musim
      sigcz(j)=sigsim
      write(50,412)zc(j),mcz(j),mucz(j),sigcz(j)
      write(*,412)zc(j),mcz(j),mucz(j),sigcz(j)
  412 format(f6.2,1x,'particle concentration probability function',/,&
     &'Symmetrization exponent',f8.3,/,&
     &'Median',e12.4,/,&
     &'Standard deviation',e12.4,//)
      enddo
      endif
!     Calculation of function values at a desired percentile
  405 write(*,*)'Do you want to calculate the function values?'
      write(*,*)'1:yes 2:no'
      read(*,*)nsrch
      select case (nsrch)
      case (1)
      write(*,*)'###### CALCULATION OF FUNCTION VALUES ######'
      write(50,*)'###### CALCULATION OF FUNCTION VALUES ######'
  400 write(*,*)'Select the probability function'
      write(*,*)'1: Dynamic pressure 10 m'
      write(*,*)'2: Particle concentration 2 m'
      write(*,*)'3: Dynamic pressure at user requested heights'
      write(*,*)'4: Particle concentration at user requested heights'
      read(*,*)npfunc
      select case (npfunc)
      case (1)
      musim=mup10
      sigsim=sigp10
      mm=mp10
  401 write(*,*)'Write the percentile'
      read(*,*)pcx
      if(pcx.lt.0.d0.or.pcx.gt.1.d0) then
      write(*,*)'Wrong percentile (0 < p <1)'
      goto 401
      else
      if(mm.lt.0.d0) then
      px=1.d0-pcx
      else
      px=pcx
      endif
      nfunc=17
      zstd=zbrent(-3.d0,3.d0)
      val=zstd*sigsim+musim
      if(val.le.0.d0) then
      write(*,*)'Warning!!'
      write(*,*)'The percentile is outside the range of calculation'
      goto 401
      else
      val=val**(1.d0/mp10)
      write(*,*)'Dynamic pressure 10 m'
      write(*,420)pcx,val
      write(50,*)'Dynamic pressure 10 m'
      write(50,420)pcx,val
      endif
      endif
      case (2)
      musim=muc2
      sigsim=sigc2
      mm=mc2
  404 write(*,*)'Write the percentile'
      read(*,*)pcx
      if(pcx.lt.0.d0.or.pcx.gt.1.d0) then
      write(*,*)'Wrong percentile (0 < p <1)'
      goto 404
      else
      if(mm.lt.0.d0) then
      px=1.d0-pcx
      else
      px=pcx
      endif
      nfunc=17
      zstd=zbrent(-3.d0,3.d0)
      val=zstd*sigsim+musim
      if(val.le.0.d0) then
      write(*,*)'Warning!!'
      write(*,*)'The percentile is outside the range of calculation'
      goto 404
      else
      val=val**(1.d0/mc2)
      write(*,*)'Particle concentration 2 m'
      write(*,421)pcx,val
      write(50,*)'Particle concentration 2 m'
      write(50,421)pcx,val
      endif
      endif
      case (3)
      if(ndynpr.eq.2) then
      write(*,*)'Warning. No user requested heights were provided'
      write(*,*)''
      goto 400
      endif
  430 write(*,*)'Choose an height'
      do j=1,ipr
      write(*,424)j,zdynpr(j)
      enddo
      read(*,*)njpr
      if(njpr.gt.ipr.or.njpr.lt.1) then
      write(*,*)'Wrong choice'
      goto 430
      endif
      musim=mupz(njpr)
      sigsim=sigpz(njpr)
      mm=mpz(njpr)
  431 write(*,*)'Write the percentile'
      read(*,*)pcx
      if(pcx.lt.0.d0.or.pcx.gt.1.d0) then
      write(*,*)'Wrong percentile (0 < p <1)'
      goto 431
      else
      if(mm.lt.0.d0) then
      px=1.d0-pcx
      else
      px=pcx
      endif
      nfunc=17
      zstd=zbrent(-3.d0,3.d0)
      val=zstd*sigsim+musim
      if(val.le.0.d0) then
      write(*,*)'Warning!!'
      write(*,*)'The percentile is outside the range of calculation'
      goto 431
      else
      val=val**(1.d0/mpz(njpr))
      write(*,422)zdynpr(njpr),pcx,val
      write(50,422)zdynpr(njpr),pcx,val
      endif
      endif
      case (4)
      if(nc.eq.2) then
      write(*,*)'Warning. No user requested heights were provided'
      write(*,*)''
      goto 400
      endif
  432 write(*,*)'Choose an height'
      do j=1,ic
      write(*,424)j,zc(j)
      enddo
      read(*,*)njc
      if(njc.gt.ipr.or.njc.lt.1) then
      write(*,*)'Wrong choice'
      goto 432
      endif
      musim=mucz(njc)
      sigsim=sigcz(njc)
      mm=mcz(njc)
  433 write(*,*)'Write the percentile'
      read(*,*)pcx
      if(pcx.lt.0.d0.or.pcx.gt.1.d0) then
      write(*,*)'Wrong percentile (0 < p <1)'
      goto 433
      else
      if(mm.lt.0.d0) then
      px=1.d0-pcx
      else
      px=pcx
      endif
      nfunc=17
      zstd=zbrent(-3.d0,3.d0)
      val=zstd*sigsim+musim
      if(val.le.0.d0) then
      write(*,*)'Warning!!'
      write(*,*)'The percentile is outside the range of calculation'
      goto 433
      else
      val=val**(1.d0/mcz(njc))
      write(*,423)zc(njc),pcx,val
      write(50,423)zc(njc),pcx,val
      endif
      endif
      case default
      write(*,*)'Wrong choice'
      goto 400
      end select
      case (2)
      return
      case default
      write(*,*)'Wrong choice'
      goto 405
      end select
  420 format('Percentile',f12.3,/,&
     &'Function value',f14.4)
  421 format('Percentile',f12.3,/,&
     &'Function value',e14.4)
  422 format(f6.2,2x,'average dynamic pressure (Pa)',/,&
     &'Percentile',f12.3,/,&
     &'Function value',f14.4)
  423 format(f6.2,2x,'particle concentration',/,&
     &'Percentile',f12.3,/,&
     &'Function value',e14.4)
  424 format(i2,2x,'z = ',f6.2)
  425 format(f6.2,1x,'Pdyn simmetrization exponent calc. residuals')
  426 format(f6.2,1x,'C simmetrization exponent calc. residuals')
      goto 405
      end subroutine probfunction
