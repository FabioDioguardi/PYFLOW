      subroutine grainsize
      USE inoutdata; USE nrtype
      implicit none
      INTERFACE
               function rhofunc(ii,d)
               USE inoutdata; USE nrtype
               implicit none
               integer :: ii
               real(dp) :: d,rhofunc
               end function rhofunc

!               function shapefunc(d)
!               USE inoutdata; USE nrtype
!               implicit none
!               real(dp) :: d,shapefunc
!               end function shapefunc
               
               character(len=20) function str(k)
               implicit none
               integer, intent(in) :: k
               end function str
               
               function dlog2(xx)
               USE nrtype
               implicit none
               REAL(dp) :: dlog2
               REAL(dp), INTENT(IN) :: xx
               end function dlog2
               
      END INTERFACE
      real(dp) :: normd
      real(dp), dimension(30) :: wcum,wcumpr,phi,diam,wgr
      real(dp),dimension(3) :: var,phivar,dvar
      real(dp), dimension(6) :: d50mm_or                                ! Real d50, without normalization
      real(dp) :: m,w50,w84,w16,phifx,phiinf,phisup,sigma,winf,wsup,wtot,z1,z2,rtflsp
      integer :: nmax,i,j,k,kk

      if(.not.distr1) then
        if(davgeqsph(1).eq.UNDEFINED) then
        write(*,*)'WARNING! Command DAVGEQSPH(1) missing in input.dat'
        write(*,*)'Particle median size will be not normalized'
            if(d50mm(1).eq.UNDEFINED) then
            d50mm(1)=2.d0**(-phi50(1))
            else
            phi50(1)=-dlog2(d50mm(1))
            endif
        else
            if(phi50(1).eq.UNDEFINED) then
            phi50(1)=-dlog2(d50mm(1))
            d50mm(1)=normd(phi50(1),sorting(1),davgeqsph(1))
            else
            d50mm(1)=normd(phi50(1),sorting(1),davgeqsph(1))
            endif
        endif
      else
      endif

      if(.not.distr2) then
        if(model.eq.'TWOCOMPONENTS') then
                if(davgeqsph(2).eq.UNDEFINED) then
                write(*,*)'WARNING! Command DAVGEQSPH(2) missing in input.dat'
                write(*,*)'Particle median size will be not normalized'
                    if(d50mm(2).eq.UNDEFINED) then
                    d50mm(2)=2.d0**(-phi50(2))
                    else
                    phi50(2)=-dlog2(d50mm(2))
                    endif
                else
                    if(phi50(2).eq.UNDEFINED) then
                    phi50(2)=-dlog2(d50mm(2))
                    d50mm(2)=normd(phi50(2),sorting(1),davgeqsph(2))
                    else
                    d50mm(2)=normd(phi50(2),sorting(1),davgeqsph(2))
                    endif
                endif
        endif
      else
      endif

      if(.not.deprates) then
      !     Set ncomp in case deprates=.false. and it was not set in the input file or to overwrite wrong input (e.g. NCOMP>2)
        if(model.eq.'TWOCOMPONENTS') then
        ncomp=2
        else
        ncomp=1
        endif
      else
      endif

      if((.not.distr1).and.(.not.distr2)) then
        if(.not.deprates) then
        goto 98
        else
        endif
      else
      endif
      
      if(input_weight.eq.'MASS') then
      psreadtot=0.d0
      if(wtot_sample.ne.undefined) goto 99
      wtot_sample=0.d0
      do i=1,ncomp
         do j=1,nclass(i)
         wtot_sample=wtot_sample+weight(i,j)
         enddo
      enddo
  99  do i=1,ncomp
      wtot=0.d0
      ! Number of classes calculation performed in check_data routine
      ! Store the total mass of grainsize classes for each component for calculating weight fraction of each component for grainsize distribution calculations
        do j=1,nclass(i)
        wgr(j)=weight(i,j)
        wtot=wtot+weight(i,j)
        enddo
        !     Weight fraction (%)
        do j=1,nclass(i)
        weight(i,j)=weight(i,j)/wtot
        enddo
      ! Calculation of weight fractions over the whole sample and storage of number of classes for deposition
        do j=1,nclass(i)
        psread(i,j)=wgr(j)/wtot_sample
        psreadtot=psreadtot+psread(i,j)
        enddo
        nclass_dep(i)=nclass(i)
      enddo
      else
      psreadtot=0.d0
      do i=1,ncomp
      wtot=0.d0
      ! Number of classes calculation performed in check_data routine
      ! Store the total weight fraction of grainsize classes for each component for renormalizing weight fraction of each component for grainsize distribution calculations
        do j=1,nclass(i)
        wtot=wtot+weight(i,j)
        enddo
      ! Calculation of weight fraction over the whole sample for deposition
        do j=1,nclass(i)
        psread(i,j)=weight(i,j)/100.d0
        psreadtot=psreadtot+psread(i,j)
        enddo
        !     Weight fraction (%)
        do j=1,nclass(i)
        weight(i,j)=weight(i,j)/wtot
        enddo
        nclass_dep(i)=nclass(i)
      enddo
      endif
      do i=1,ncomp
!      open(57,file=trim('component'//str(i))//'.dat')
      write(52,200)i
!      write(57,200)i
      write(*,200)i
        !     Cumulative distribution creation
      wcum(1)=weight(i,1)
      wcumpr(1)=wcum(1)*100.d0
      phi(1)=phimin(i)
      diam(1)=1.d0/(2.d0**(phi(1)))
      dread(i,1)=(2.d0**(-phi(1)))/1000.d0                              ! For deposition
      if(rhos(i,1).eq.undefined) rhos(i,1)=rhofunc(i,diam(1))
!      if(shapefact(i,1).eq.undefined) shapefact(i,1)=shapefunc(diam(1))
!      write(57,*)'Phi   d(mm)    Wt.%   Fr.Cum.Wt     Cum.Wt%'
!      write(57,100)phi(1),diam(1),weight(i,1)*100.d0,wcum(1),wcumpr(1)
        do j=2,nclass(i)
        wcum(j)=wcum(j-1)+weight(i,j)
        wcumpr(j)=wcum(j)*100.d0
        phi(j)=phi(j-1)+dphi(i)
        diam(j)=1.d0/(2.d0**(phi(j)))
        dread(i,j)=(2.d0**(-phi(j)))/1000.d0                            ! For deposition
        if(rhos(i,j).eq.undefined) rhos(i,j)=rhofunc(i,diam(j))
!        if(shapefact(i,j).eq.undefined) shapefact(i,j)=shapefunc(diam(j))
!        write(57,100)phi(j),diam(j),weight(i,j)*100.d0,wcum(j),wcumpr(j)
        enddo
!     Search for the interval d50
      write(52,*)'***GRAINSIZE ANALYSIS CALCULATIONS***'
      write(*,*)'***GRAINSIZE ANALYSIS CALCULATIONS***'
      w50=0.5d0
      var(1)=w50
      w84=0.84d0
      var(2)=w84
      w16=0.16d0
      var(3)=w16
      nfunc=19
            do kk=1,3
                  do j=1,nclass(i)-1
                  if(var(kk).ge.wcum(j).and.var(kk).lt.wcum(j+1)) k=j
                  enddo
            winf=wcum(k)
            phiinf=phi(k)
            wsup=wcum(k+1)
            phisup=phi(k+1)
            !     Search for phisnr value from inverted gaussian
            fxdistr=wsup
            write(52,120)var(kk)
            write(*,120)var(kk)
            z2=rtflsp(phi(k),phi(k+2))
            !     Search for phiinr value from inverted gaussian
            fxdistr=winf
            write(52,121)var(kk)
            write(*,121)var(kk)
            z1=rtflsp(phi(k-1),phi(k+1))
            !     Corresponding phi value
            fxdistr=var(kk)                                                         !per distinguere le funzioni chiamate in rtflsp
            write(52,122)var(kk)
            write(*,122)var(kk)
            phifx=rtflsp(phiinf,phisup)
            m=(z2-z1)/(phisup-phiinf)
            phivar(kk)=(phifx-z1+m*phiinf)/m
            dvar(kk)=1.d0/(2.d0**(phivar(kk)))
            enddo
        sigma=(phivar(2)-phivar(3))/2.d0
!        write(57,*)'PARTICLE 1 GRAINSIZE ANALYSIS RESULTS'
!        write(57,101)phivar(1),dvar(1),phivar(2),dvar(2),phivar(3),dvar(3),sigma
        write(*,102)i
        write(*,101)phivar(1),dvar(1),phivar(2),dvar(2),phivar(3),dvar(3),sigma
        write(52,102)i
        write(52,101)phivar(1),dvar(1),phivar(2),dvar(2),phivar(3),dvar(3),sigma
        phi50(i)=phivar(1)
        d50mm_or(i)=2.d0**(-phi50(i))                                   ! Store d50mm for printing it into the component grainsize analysis files, where the real d50 has to be printed and not the normalized one
        sorting(i)=sigma
        phi84(i)=phivar(2)
        d84mm(i)=dvar(2)
        phi16(i)=phivar(3)
        d16mm(i)=dvar(3)
        select case (model)
            case ('TWOLAYERS')
                if(i.eq.1) then
!                    if(davgeqsph(1).eq.UNDEFINED) then
                    if(davgeqsph(i).eq.UNDEFINED) then
                    write(*,103)i
                    write(52,103)i
                    !d50mm(1)=2.d0**(-phi50(1))
                    d50mm(i)=2.d0**(-phi50(i))
                    else
                    !d50mm(1)=normd(phi50(1),sorting(1),davgeqsph(1))
                    d50mm(i)=normd(phi50(i),sorting(i),davgeqsph(i))
                    endif
                else
                endif

            case ('TWOCOMPONENTS')
                 if(i.le.2) then
                    if(davgeqsph(i).eq.UNDEFINED) then
                    write(*,103)i
                    write(52,103)i
                    d50mm(i)=2.d0**(-phi50(i))
                    else
                    d50mm(i)=normd(phi50(i),sorting(i),davgeqsph(i))
                    endif
                 else
                 endif
!                 if(i.eq.1) then
!                    if(davgeqsph(1).eq.UNDEFINED) then
!                    write(*,*)'WARNING! Command DAVGEQSPH(1) missing in input.dat'
!                    write(*,*)'Particle median size will be not normalized'
!                    d50mm(1)=2.d0**(-phi50(1))
!                    else
!                    d50mm(1)=normd(phi50(1),sorting(1),davgeqsph(1))
!                    endif
!                else
!                endif
!                if(i.eq.2) then
!                    if(davgeqsph(2).eq.UNDEFINED) then
!                    write(*,*)'WARNING! Command DAVGEQSPH(2) missing in input.dat'
!                    write(*,*)'Particle median size will be not normalized'
!                    d50mm(2)=2.d0**(-phi50(2))
!                    else
!                    d50mm(2)=normd(phi50(2),sorting(2),davgeqsph(2))
!                    endif
!                else
!                endif
        end select

      close(57)
      ! Calculate density and shape factor of the median grainsizes if they are not given at input
      enddo
   98 if(rhos(1,0).eq.undefined) rhos(1,0)=rhofunc(1,d50mm(1))
 !     if(shapefact(1,0).eq.undefined) shapefact(1,0)=shapefunc(d50mm(1))
      if(model.eq.'TWOCOMPONENTS') then
      if(rhos(2,0).eq.undefined) rhos(2,0)=rhofunc(2,d50mm(2))
 !     if(shapefact(2,0).eq.undefined) shapefact(2,0)=shapefunc(d50mm(2))
      else
      endif
      
      if(model.eq.'TWOCOMPONENTS'.and.rhos(2,0).gt.rhos(1,0)) call component_swap
!     Just to write component summary files
      if(.not.deprates) then
         if(distr1) then
          i=1
          wcum(1)=weight(i,1)
          wcumpr(1)=wcum(1)*100.d0
          phi(1)=phimin(i)
          diam(1)=1.d0/(2.d0**(phi(1)))
          open(57,file=trim('component'//str(i))//'.dat')
          write(57,200)i
          write(57,*)'Phi   d(mm)    Wt.%   Fr.Cum.Wt     Cum.Wt%'
          write(57,100)phi(1),diam(1),weight(i,1)*100.d0,wcum(1),wcumpr(1)
                do j=2,nclass(i)
                wcum(j)=wcum(j-1)+weight(i,j)
                wcumpr(j)=wcum(j)*100.d0
                phi(j)=phi(j-1)+dphi(i)
                diam(j)=1.d0/(2.d0**(phi(j)))
                write(57,100)phi(j),diam(j),weight(i,j)*100.d0,wcum(j),wcumpr(j)
                enddo
                if(dotestchi(i)) then
                probchi=siglevchi(i)
                sensgrchi=senschi(i)
                call testchi(nclass(i),phi,weight(i,:),phi50(i),sorting(i))
                else
                write(*,201)i
                write(52,201)i
                write(57,201)i
                endif
          write(57,102)i
          write(57,101)phi50(i),d50mm_or(i),phi84(i),d84mm(i),phi16(i),d16mm(i),sorting(i)
          close(57)
          else
          endif
          if(distr2) then
          i=2
          wcum(1)=weight(i,1)
          wcumpr(1)=wcum(1)*100.d0
          phi(1)=phimin(i)
          diam(1)=1.d0/(2.d0**(phi(1)))
          open(57,file=trim('component'//str(i))//'.dat')
          write(57,200)i
          write(57,*)'Phi   d(mm)    Wt.%   Fr.Cum.Wt     Cum.Wt%'
          write(57,100)phi(1),diam(1),weight(i,1)*100.d0,wcum(1),wcumpr(1)
                do j=2,nclass(i)
                wcum(j)=wcum(j-1)+weight(i,j)
                wcumpr(j)=wcum(j)*100.d0
                phi(j)=phi(j-1)+dphi(i)
                diam(j)=1.d0/(2.d0**(phi(j)))
                write(57,100)phi(j),diam(j),weight(i,j)*100.d0,wcum(j),wcumpr(j)
                enddo
                if(dotestchi(i)) then
                probchi=siglevchi(i)
                sensgrchi=senschi(i)
                call testchi(nclass(i),phi,weight(i,:),phi50(i),sorting(i))
                else
                write(*,201)i
                write(52,201)i
                write(57,201)i
                endif
          write(57,102)i
          write(57,101)phi50(i),d50mm_or(i),phi84(i),d84mm(i),phi16(i),d16mm(i),sorting(i)
          close(57)
          else
          endif
      return
      else
      endif
      do i=1,ncomp
      wcum(1)=weight(i,1)
      wcumpr(1)=wcum(1)*100.d0
      phi(1)=phimin(i)
      diam(1)=1.d0/(2.d0**(phi(1)))
      open(57,file=trim('component'//str(i))//'.dat')
      write(57,200)i
      write(57,*)'Phi   d(mm)    Wt.%   Fr.Cum.Wt     Cum.Wt%'
      write(57,100)phi(1),diam(1),weight(i,1)*100.d0,wcum(1),wcumpr(1)
        do j=2,nclass(i)
        wcum(j)=wcum(j-1)+weight(i,j)
        wcumpr(j)=wcum(j)*100.d0
        phi(j)=phi(j-1)+dphi(i)
        diam(j)=1.d0/(2.d0**(phi(j)))
        write(57,100)phi(j),diam(j),weight(i,j)*100.d0,wcum(j),wcumpr(j)
        enddo
            if(dotestchi(i)) then
            probchi=siglevchi(i)
            sensgrchi=senschi(i)
            call testchi(nclass(i),phi,weight(i,:),phi50(i),sorting(i))
            else
            write(*,201)i
            write(52,201)i
            write(57,201)i
            endif
        write(57,102)i
        write(57,101)phi50(i),d50mm_or(i),phi84(i),d84mm(i),phi16(i),d16mm(i),sorting(i)
      close(57)
      enddo
      
  100 format(f4.1,2x,f6.3,3x,f5.2,3x,f8.4,3x,f9.4)
  120 format(/,f4.1,1x,'Phisup calculation residuals')
  121 format(/,f4.1,1x,'Phiinf calculation residuals')
  122 format(/,f4.1,1x,'Phi calculation residuals')
  101 format(/,'phi50 = ',f8.5,3x,'d50 (mm) = ',f8.5,/,'phi84 = ',f8.5,3x,'d84 (mm) = ',f8.5,/,&
      'phi16 = ',f8.5,3x,'d16 (mm) = ',f8.5,/,'sigma = ',f8.5,/)
  102 format('COMPONENT ',i2,' GRAINSIZE ANALYSIS RESULTS')
  103 format('WARNING! Command DAVGEQSPH(',i2,') missing in input.dat',/,'Particle median size will be not normalized',/)
  200 format('****** GRAINSIZE DISTRIBUTION PARAMETERS CALCULATION FOR COMPONENT ',i2,' *******',/)
  201 format(/,'WARNING! Chi squared test on the grainsize distribution not done for component ',i2)
      end subroutine grainsize

      function normd(dphi,sigma,dsph)
      USE nrtype
      implicit none
	INTERFACE
		FUNCTION cum(x)
                USE inoutdata; USE nrtype
		IMPLICIT NONE
		REAL(dp) :: cum
		REAL(dp), INTENT(IN) :: x
		END FUNCTION cum
        END INTERFACE
      real(dp) :: dphi,sigma,dsph
      real(dp) :: phiup,phidown,z1,z2,cumz1,cumz2,zdsph,cumzdsph,zdsphnorm,dsphnorm,normd
      call phisearch(dphi,phiup,phidown)
      zdsph=(dsph-dphi)/sigma
      cumzdsph=cum(zdsph)
      zdsphnorm=(zdsph*0.5d0)/cumzdsph
      dsphnorm=zdsphnorm*sigma+dphi
      normd=2.d0**(-dsphnorm)
      end function normd

      subroutine phisearch(dphi,phiup,phidown)
      USE nrtype
      implicit none
      real(dp) :: dphi,phiup,phidown,incr
      phiup=-6.d0
      phidown=-5.5d0
      incr=0.5d0
  123 if(dphi.ge.phiup.and.dphi.lt.phidown)return
      phiup=phiup+incr
      phidown=phidown+incr
      goto 123
      end subroutine phisearch
      
      function rhofunc(ii,d)
      USE inoutdata; USE nrtype
      implicit none
      integer :: ii
      real(dp) :: d,rhofunc
      select case (rholaw(ii))
             case ('POLLENA')
             rhofunc=(-0.31d0*dlog(d)+1.83d0)*1000.d0
             if(rhofunc.gt.2760.d0) rhofunc=2760.d0
             case ('AVERNO2')
             rhofunc=(-0.30d0*dlog(d)+1.42d0)*1000.d0
             if(rhofunc.gt.2600.d0) rhofunc=2600.d0
             case('AMS')
             rhofunc=(1.0599d0*d**(-0.332d0))*1000.d0
             if(rhofunc.gt.2600.d0) rhofunc=2560.d0
             case ('POMPEI')
             rhofunc=(1.8223d0*d**(-0.167d0))*1000.d0
             if(rhofunc.gt.2700.d0) rhofunc=2700.d0
             case ('MERCATO')
             rhofunc=(0.9556d0*d**(-0.191d0))*1000.d0
             if(rhofunc.gt.2400.d0) rhofunc=2400.d0
             case ('ASTRONI')
             rhofunc=(0.8128d0*d**(-0.213d0))*1000.d0
             if(rhofunc.gt.2400.d0) rhofunc=2510.d0
             case ('SIAL_XX')
             rhofunc=2400.d0
             case ('FEM_XX')
             rhofunc=3280.d0
             case ('LITHIC')
             rhofunc=2570.d0
             case ('CUSTOM')
             rhofunc=rhos(ii,0)
      end select
      end function rhofunc

!      function shapefunc(d)
!      USE inoutdata; USE nrtype
!      implicit none
!      real(dp) :: d,shapefunc
!      shapefunc=0.458d0*(d/10.d0)+0.227d0
!      end function shapefunc
      
      character(len=20) function str(k)
      implicit none
      integer, intent(in) :: k
      write(str,*)k
      str=adjustl(str)
      end function str
      
      subroutine component_swap
      use inoutdata; use nrtype
      implicit none
      real(dp) :: dcphitemp,d50mmtemp,sigmatemp,denstemp,siglevchitemp,dphitemp,phimintemp,phimaxtemp
      logical :: dotestchitemp
      character(len=10):: cdlawtemp
      real(dp), dimension(30) :: weighttemp
      real(dp),dimension(0:30):: shpar1temp,shpar2temp,shpar3temp
!      logical, dimension(0:30) :: shpar4temp
      logical :: shpar4temp
      integer :: i,nclasstemp,nclasstemp_bis
      dcphitemp=phi50(1)
      d50mmtemp=d50mm(1)
      sigmatemp=sorting(1)
      nclasstemp=nclass(1)
      nclasstemp_bis=nclass(1)
      denstemp=rhos(1,0)
      dotestchitemp=dotestchi(1)
      siglevchitemp=siglevchi(1)
      dphitemp=dphi(1)
      phimintemp=phimin(1)
      phimaxtemp=phimax(1)
      shpar4temp=shpar4(1)
      do i=1,30
      weighttemp(i)=weight(1,i)
      enddo
      do i=0,30
      shpar1temp(i)=shpar1(1,i)
      shpar2temp(i)=shpar2(1,i)
      shpar3temp(i)=shpar3(1,i)
!      shpar4temp(i)=shpar4(1,i)
      enddo
      cdlawtemp=cdlaw(1)

      phi50(1)=phi50(2)
      d50mm(1)=d50mm(2)
      sorting(1)=sorting(2)
      nclass(1)=nclass(2)
      rhos(1,0)=rhos(2,0)
      dotestchi(1)=dotestchi(2)
      siglevchi(1)=siglevchi(2)
      dphi(1)=dphi(2)
      phimin(1)=phimin(2)
      phimax(1)=phimax(2)
      shpar4(1)=shpar4(2)
      do i=1,30
      weight(1,i)=weight(2,i)
      enddo
      do i=0,30
      shpar1(1,i)=shpar1(2,i)
      shpar2(1,i)=shpar2(2,i)
      shpar3(1,i)=shpar3(2,i)
!      shpar4(1,i)=shpar4(2,i)
      enddo
      cdlaw(1)=cdlaw(2)

      phi50(2)=dcphitemp
      d50mm(2)=d50mmtemp
      sorting(2)=sigmatemp
      nclass(2)=nclasstemp
      nclass(2)=nclasstemp_bis
      rhos(2,0)=denstemp
      dotestchi(2)=dotestchitemp
      siglevchi(2)=siglevchitemp
      dphi(2)=dphitemp
      phimin(2)=phimintemp
      phimax(2)=phimaxtemp
      shpar4(2)=shpar4temp
      do i=1,30
      weight(2,i)=weighttemp(i)
      enddo
      do i=0,30
      shpar1(2,i)=shpar1temp(i)
      shpar2(2,i)=shpar2temp(i)
      shpar3(2,i)=shpar3temp(i)
!      shpar4(2,i)=shpar4temp(i)
      enddo
      cdlaw(2)=cdlawtemp
      end subroutine component_swap

