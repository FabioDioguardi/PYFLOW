      subroutine probfunction
      USE inoutdata; USE nrtype
      implicit none
	INTERFACE
		FUNCTION zbrent(x1,x2)
                USE inoutdata; USE nrtype; USE nrutil, ONLY : nrerror
                IMPLICIT NONE
	        REAL(dp) :: x1,x2
         	REAL(dp) :: zbrent
         	END FUNCTION zbrent
	END INTERFACE
      real(dp), dimension(20) :: mpz,mupz,sigpz,mcz,mucz,sigcz
      real(dp) :: temp1,temp2,val,tempz1,tempz2,zstd,mp10,mc2,mup10,muc2,sigc2,sigp10
	  real(dp) :: mden, muden, sigden, mush, muush, sigush
      real(dp) :: mrtot,mtdep,murtot,sigrtot,mutdep,sigtdep
      real(dp) :: x1tmp,x2tmp,x3tmp,x4tmp
      integer :: j,l,njc,njpr
!      x1tmp=-50.d0;x2tmp=-1.d-10;x3tmp=1.d-10;x4tmp=50.d0
      x1tmp=-10.d0;x2tmp=-1.d-10;x3tmp=1.d-10;x4tmp=10.d0

      if(only_deprates.and.n_solutions.lt.3) then
      write(*,*)'WARNING! It is not possible to calculate PDFs of deposition rate and time'
      write(*,*)'N_SOLUTIONS must be at least 3'
      write(52,*)'WARNING! It is not possible to calculate PDFs of deposition rate and time'
      write(52,*)'N_SOLUTIONS must be at least 3'
      return
      else
      endif
      if(only_deprates) goto 599
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
!      if(usr_z_dynpr.eqv..FALSE..and.usr_z_c.eqv..FALSE.) goto 405


!     Symmetric probability distribution parameters for flow density
      mudstr=dennrm
      mxdstr=denmax
      mndstr=denmin
      write(*,*)'Density simmetrization exponent calc. residuals'
      write(52,*)'Density simmetrization exponent calc. residuals'
  602 temp1=zbrent(x3tmp,x4tmp)
      if(checkzbr) temp1=0.d0
      write(52,369)temp1
      write(*,369)temp1
                if(temp1.ne.0.d0.and.temp1.ge.1.d-9.and.temp1.ne.50.d0) then
                mden=temp1
                else
                temp2=zbrent(x1tmp,x2tmp)
                     if(checkzbr) then
                     x1tmp=x1tmp-50.d0
                     x2tmp=x2tmp-50.d0
                     x3tmp=x3tmp+50.d0
                     x4tmp=x4tmp+50.d0
                     goto 602
                     endif
                write(52,369)temp2
                write(*,369)temp2
                     if(temp2.ne.0.d0.and.temp2.le.-1.d-9.and.temp2.ne.-50.d0) then
                     mden=temp2
                     else
                     write(*,*)'Warning. Unable to find a simmetrization coefficient'
                     write(*,*)'for flow density prob. function'
                     endif
                endif
      musim=mudstr**mden
      sigsim=abs(mxdstr**mden-mudstr**mden)
      muden=musim
      sigden=sigsim
      write(50,*)'Flow density probability function'
      write(50,410)mden,muden,sigden
      write(*,*)'Flow density probability function'
      write(*,410)mden,muden,sigden


!     Symmetric probability distribution parameters for flow shear velocity
      mudstr=ushavg
      mxdstr=ushmax
      mndstr=ushmin
      write(*,*)'Density simmetrization exponent calc. residuals'
      write(52,*)'Density simmetrization exponent calc. residuals'
  603 temp1=zbrent(x3tmp,x4tmp)
      if(checkzbr) temp1=0.d0
      write(52,369)temp1
      write(*,369)temp1
                if(temp1.ne.0.d0.and.temp1.ge.1.d-9.and.temp1.ne.50.d0) then
                mush=temp1
                else
                temp2=zbrent(x1tmp,x2tmp)
                     if(checkzbr) then
                     x1tmp=x1tmp-50.d0
                     x2tmp=x2tmp-50.d0
                     x3tmp=x3tmp+50.d0
                     x4tmp=x4tmp+50.d0
                     goto 603
                     endif
                write(52,369)temp2
                write(*,369)temp2
                     if(temp2.ne.0.d0.and.temp2.le.-1.d-9.and.temp2.ne.-50.d0) then
                     mush=temp2
                     else
                     write(*,*)'Warning. Unable to find a simmetrization coefficient'
                     write(*,*)'for flow shear velocity prob. function'
                     endif
                endif
      musim=mudstr**mush
      sigsim=abs(mxdstr**mush-mudstr**mush)
      muush=musim
      sigush=sigsim
      write(50,*)'Flow shear velocity probability function'
      write(50,410)mush,muush,sigush
      write(*,*)'Flow shear velocity probability function'
      write(*,410)mush,muush,sigush



      if(usr_z_dynpr) then
!     Determination of probability function of Pdyn at user requested heights
      write(*,*)'###PROBABILITY FUNCTIONS FOR AVERAGE DYNAMIC PRESSURE OVER USER REQUESTED HEIGHTS###'
      write(52,*)'###PROBABILITY FUNCTIONS FOR AVERAGE DYNAMIC PRESSURE OVER USER REQUESTED HEIGHTS###'
                if(ipr.eq.0) then
                write(*,*)'FATAL ERROR! Command ZDYNPR(i) missing in input.dat'
                write(*,*)''
                write(52,*)'FATAL ERROR! Command ZDYNPR(i) missing in input.dat'
                write(52,*)''
                write(*,*)''
                write(*,*)'PYFLOW 2.0 is going to be stopped'
                write(52,*)''
                write(52,*)'PYFLOW 2.0 is going to be stopped'
                stop
                else
                endif
      do j=1,ipr
      mudstr=pzavg(j)
      mxdstr=pzmax(j)
      mndstr=pzmin(j)
      write(*,425)zdynpr(j)
      write(52,425)zdynpr(j)
  604 tempz1=zbrent(x3tmp,x4tmp)
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
                      goto 604
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
      enddo
      endif


      if(usr_z_c) then
!     Determination of probability function of C at user requested heights
      write(*,*)'###PROBABILITY FUNCTIONS FOR PARTICLE CONCENTRATION AT USER REQUESTED HEIGHTS###'
      write(52,*)'###PROBABILITY FUNCTIONS FOR PARTICLE CONCENTRATION AT USER REQUESTED HEIGHTS###'
                if(ic.eq.0) then
                write(*,*)'FATAL ERROR! Command ZC(i) missing in input.dat'
                write(*,*)''
                write(52,*)'FATAL ERROR! Command ZC(i) missing in input.dat'
                write(52,*)''
                write(*,*)''
                write(*,*)'PYFLOW 2.0 is going to be stopped'
                write(52,*)''
                write(52,*)'PYFLOW 2.0 is going to be stopped'
                stop
                else
                endif
      do j=1,ic
      mudstr=czavg(j)
      mxdstr=czmax(j)
      mndstr=czmin(j)
      write(*,426)zdynpr(j)
      write(52,426)zdynpr(j)
  605 tempz1=zbrent(x3tmp,x4tmp)
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
                      goto 605
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
      enddo
      endif


!     Symmetric probability distribution parameters for Rtot
  599 if(deprates) then
      mudstr=rtot_avg
      mxdstr=rtot_max
      mndstr=rtot_min
      nfunc=18
      write(*,*)'Rtot simmetrization exponent calc. residuals'
      write(52,*)'Rtot simmetrization exponent calc. residuals'
  606 temp1=zbrent(x3tmp,x4tmp)
      if(checkzbr) temp1=0.d0
      write(52,369)temp1
      write(*,369)temp1
             if(temp1.ne.0.d0.and.temp1.ge.1.d-9.and.temp1.ne.50.d0) then
             mrtot=temp1
             else
             temp2=zbrent(x1tmp,x2tmp)
                if(checkzbr) then
                x1tmp=x1tmp-50.d0
                x2tmp=x2tmp-50.d0
                x3tmp=x3tmp+50.d0
                x4tmp=x4tmp+50.d0
                goto 606
                endif
             write(52,369)temp2
             write(*,369)temp2
                if(temp2.ne.0.d0.and.temp2.le.-1.d-9.and.temp2.ne.-50.d0) then
                mrtot=temp2
                else
                write(*,*)'Warning. Unable to find a simmetrization coefficient'
                write(*,*)'for 10 m dynamic pressure prob. function'
                endif
             endif
      musim=mudstr**mrtot
      sigsim=abs(mxdstr**mrtot-mudstr**mrtot)
      murtot=musim
      sigrtot=sigsim
      write(50,*)'Rtot probability function'
      write(50,410)mrtot,murtot,sigrtot
      write(*,*)'Rtot probability function'
      write(*,410)mrtot,murtot,sigrtot

!     Symmetric probability distribution parameters for tdep
      mudstr=tdep_avg
      mxdstr=tdep_max
      mndstr=tdep_min
      nfunc=18
      write(*,*)'tdep simmetrization exponent calc. residuals'
      write(52,*)'tdep simmetrization exponent calc. residuals'
  607 temp1=zbrent(x3tmp,x4tmp)
      if(checkzbr) temp1=0.d0
      write(52,369)temp1
      write(*,369)temp1
             if(temp1.ne.0.d0.and.temp1.ge.1.d-9.and.temp1.ne.50.d0) then
             mtdep=temp1
             else
             temp2=zbrent(x1tmp,x2tmp)
                if(checkzbr) then
                x1tmp=x1tmp-50.d0
                x2tmp=x2tmp-50.d0
                x3tmp=x3tmp+50.d0
                x4tmp=x4tmp+50.d0
                goto 607
                endif
             write(52,369)temp2
             write(*,369)temp2
                if(temp2.ne.0.d0.and.temp2.le.-1.d-9.and.temp2.ne.-50.d0) then
                mtdep=temp2
                else
                write(*,*)'Warning. Unable to find a simmetrization coefficient'
                write(*,*)'for 10 m dynamic pressure prob. function'
                endif
             endif
      musim=mudstr**mtdep
      sigsim=abs(mxdstr**mtdep-mudstr**mtdep)
      mutdep=musim
      sigtdep=sigsim
      write(50,*)'tdep probability function'
      write(50,410)mtdep,mutdep,sigtdep
      write(*,*)'tdep probability function'
      write(*,410)mtdep,mutdep,sigtdep
      else
      endif


!     Calculation of function values at a desired percentile
      if(.not.usr_pcx_sol) return
      write(*,*)'###FUNCTIONS VALUES AT USER REQUESTED PERCENTILES###'
      write(50,*)'###FUNCTIONS VALUES AT USER REQUESTED PERCENTILES###'
      write(52,*)'###FUNCTIONS VALUES AT USER REQUESTED PERCENTILES###'
      l=1
      do while(pcx(l).ne.undefined)
                if(pcx(l).lt.0.d0.or.pcx(l).gt.1.d0) then
                write(*,428)l
                write(52,428)l
                write(*,*)'PYFLOW 2.0 is going to be stopped'
                write(52,*)'PYFLOW 2.0 is going to be stopped'
                stop
                else
                endif
      write(*,427)pcx(l)
      write(50,427)pcx(l)
      write(52,427)pcx(l)
      if(only_deprates) goto 404
      write(*,*)'Dynamic pressure 10 m'
      write(52,*)'Dynamic pressure 10 m'
      musim=mup10
      sigsim=sigp10
      mm=mp10
                if(mm.lt.0.d0) then
                px=1.d0-pcx(l)
                else
                px=pcx(l)
                endif
      nfunc=17
      zstd=zbrent(-4.d0,4.d0)
      val=zstd*sigsim+musim
                if(val.le.0.d0) then
                write(*,*)'Warning!!'
                write(*,*)'The percentile is outside the range of calculation'
                write(52,*)'Warning!!'
                write(52,*)'The percentile is outside the range of calculation'
                else
                val=val**(1.d0/mp10)
                write(*,*)'Dynamic pressure 10 m'
                write(*,420)pcx(l),val
                write(50,*)'Dynamic pressure 10 m'
                write(50,420)pcx(l),val
                endif
      write(*,*)'Particle concentration 2 m'
      write(52,*)'Particle concentration 2 m'
      musim=muc2
      sigsim=sigc2
      mm=mc2
                if(mm.lt.0.d0) then
                px=1.d0-pcx(l)
                else
                px=pcx(l)
                endif
      nfunc=17
      zstd=zbrent(-4.d0,4.d0)
      val=zstd*sigsim+musim
                if(val.le.0.d0) then
                write(*,*)'Warning!!'
                write(*,*)'The percentile is outside the range of calculation'
                write(52,*)'Warning!!'
                write(52,*)'The percentile is outside the range of calculation'
                else
                val=val**(1.d0/mc2)
                write(*,*)'Particle concentration 2 m'
                write(*,421)pcx(l),val
                write(50,*)'Particle concentration 2 m'
                write(50,421)pcx(l),val
                endif
      write(*,*)'Flow density'
      write(52,*)'Flow density'
      musim=muden
      sigsim=sigden
      mm=mden
                if(mm.lt.0.d0) then
                px=1.d0-pcx(l)
                else
                px=pcx(l)
                endif
      nfunc=17
      zstd=zbrent(-4.d0,4.d0)
      val=zstd*sigsim+musim
                if(val.le.0.d0) then
                write(*,*)'Warning!!'
                write(*,*)'The percentile is outside the range of calculation'
                write(52,*)'Warning!!'
                write(52,*)'The percentile is outside the range of calculation'
                else
                val=val**(1.d0/mden)
                write(*,*)'Flow density'
                write(*,420)pcx(l),val
                write(50,*)'Flow density'
                write(50,420)pcx(l),val
                endif
      write(*,*)'Flow shear velocity'
      write(52,*)'Flow shear velocity'
      musim=muush
      sigsim=sigush
      mm=mush
                if(mm.lt.0.d0) then
                px=1.d0-pcx(l)
                else
                px=pcx(l)
                endif
      nfunc=17
      zstd=zbrent(-4.d0,4.d0)
      val=zstd*sigsim+musim
                if(val.le.0.d0) then
                write(*,*)'Warning!!'
                write(*,*)'The percentile is outside the range of calculation'
                write(52,*)'Warning!!'
                write(52,*)'The percentile is outside the range of calculation'
                else
                val=val**(1.d0/mush)
                write(*,*)'Flow shear velocity'
                write(*,420)pcx(l),val
                write(50,*)'Flow shear velocity'
                write(50,420)pcx(l),val
                endif
      if(.not.usr_z_dynpr) goto 608
      write(*,*)'Dynamic pressure at user requested heights'
      write(52,*)'Dynamic pressure at user requested heights'
      do j=1,ipr
      write(*,424)j,zdynpr(j)
      write(52,424)j,zdynpr(j)
      musim=mupz(j)
      sigsim=sigpz(j)
      mm=mpz(j)
                if(mm.lt.0.d0) then
                px=1.d0-pcx(l)
                else
                px=pcx(l)
                endif
      nfunc=17
      zstd=zbrent(-4.d0,4.d0)
      val=zstd*sigsim+musim
                if(val.le.0.d0) then
                write(*,*)'Warning!!'
                write(*,*)'The percentile is outside the range of calculation'
                write(52,*)'Warning!!'
                write(52,*)'The percentile is outside the range of calculation'
                else
                val=val**(1.d0/mpz(j))
                write(*,422)zdynpr(j),pcx(l),val
                write(50,422)zdynpr(j),pcx(l),val
                endif
      enddo
  608 if(.not.usr_z_c) goto 404
      write(*,*)'Particle concentration at user requested heights'
      write(52,*)'Particle concentration at user requested heights'
          do j=1,ic
          write(*,424)j,zc(j)
          write(52,424)j,zc(j)
          musim=mucz(j)
          sigsim=sigcz(j)
          mm=mcz(j)
                if(mm.lt.0.d0) then
                px=1.d0-pcx(l)
                else
                px=pcx(l)
                endif
          nfunc=17
          zstd=zbrent(-4.d0,4.d0)
          val=zstd*sigsim+musim
                if(val.le.0.d0) then
                write(*,*)'Warning!!'
                write(*,*)'The percentile is outside the range of calculation'
                write(52,*)'Warning!!'
                write(52,*)'The percentile is outside the range of calculation'
                else
                val=val**(1.d0/mcz(j))
                write(*,423)zc(j),pcx(l),val
                write(50,423)zc(j),pcx(l),val
                endif
          enddo
  404 if(.not.deprates) goto 405 
      write(*,*)'Total deposition rate'
      write(52,*)'Total deposition rate'
      musim=murtot
      sigsim=sigrtot
      mm=mrtot
                if(mm.lt.0.d0) then
                px=1.d0-pcx(l)
                else
                px=pcx(l)
                endif
      nfunc=17
      zstd=zbrent(-4.d0,4.d0)
      val=zstd*sigsim+musim
                if(val.le.0.d0) then
                write(*,*)'Warning!!'
                write(*,*)'The percentile is outside the range of calculation'
                write(52,*)'Warning!!'
                write(52,*)'The percentile is outside the range of calculation'
                else
                val=val**(1.d0/mrtot)
                write(*,*)'Total deposition rate'
                write(*,421)pcx(l),val
                write(50,*)'Total deposition rate'
                write(50,421)pcx(l),val
                endif
      write(*,*)'Total deposition time'
      write(52,*)'Total deposition time'
      musim=mutdep
      sigsim=sigtdep
      mm=mtdep
                if(mm.lt.0.d0) then
                px=1.d0-pcx(l)
                else
                px=pcx(l)
                endif
      nfunc=17
      zstd=zbrent(-4.d0,4.d0)
      val=zstd*sigsim+musim
                if(val.le.0.d0) then
                write(*,*)'Warning!!'
                write(*,*)'The percentile is outside the range of calculation'
                write(52,*)'Warning!!'
                write(52,*)'The percentile is outside the range of calculation'
                return
                else
                val=val**(1.d0/mtdep)
                write(*,*)'Total deposition time'
                write(*,421)pcx(l),val
                write(50,*)'Total deposition time'
                write(50,421)pcx(l),val
                endif
  405 l=l+1
      enddo
      
  369 format('Temp. simmetrization exponent',f8.3,/)
  410 format('Symmetrization exponent',f8.3,/,&
     &'Median',e12.4,/,&
     &'Standard deviation',e12.4,//)
  411 format(f6.2,1x,'m dynamic pressure probability function',/,&
     &'Symmetrization exponent',f8.3,/,&
     &'Median',e12.4,/,&
     &'Standard deviation',e12.4,//)
  412 format(f6.2,1x,'m particle concentration probability function',/,&
     &'Symmetrization exponent',f8.3,/,&
     &'Median',e12.4,/,&
     &'Standard deviation',e12.4,//)
  420 format('Percentile',f12.3,/,&
     &'Function value',f14.4,/)
  421 format('Percentile',f12.3,/,&
     &'Function value',e14.4,/)
  422 format(f6.2,1x,'m',2x,'average dynamic pressure (Pa)',/,&
     &'Percentile',f12.3,/,&
     &'Function value',f14.4,/)
  423 format(f6.2,1x,'m',2x,'particle concentration',/,&
     &'Percentile',f12.3,/,&
     &'Function value',e14.4,/)
  424 format(i2,2x,'z = ',f6.2)
  425 format(f6.2,1x,'Pdyn simmetrization exponent calc. residuals')
  426 format(f6.2,1x,'C simmetrization exponent calc. residuals')
  427 format('Percentile = ',f6.3,/)
  428 format(/,'FATAL ERROR! PCX('i2') > 1 or < 0',/)
      end subroutine probfunction
