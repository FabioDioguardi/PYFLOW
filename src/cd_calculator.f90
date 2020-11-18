      function cdmodel(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      REAL(dp), INTENT(INOUT) :: dpart,rhop,rho
      REAL(dp) :: sphere,haidlev,swamoj,ganser,chien,trancong,dellino,holzsomm,diogmele,diog2017,diog2018
      REAL(dp) :: cdmodel
          select case (cdlaw(ishape))
          case ('SPHERE')
          cdmodel=sphere(dpart,rhop,rho)
          case ('HAIDLEV')
          cdmodel=haidlev(dpart,rhop,rho)
          case ('SWAMOJ')
          cdmodel=swamoj(dpart,rhop,rho)
          case ('GANSER')
          cdmodel=ganser(dpart,rhop,rho)
          case ('CHIEN')
          cdmodel=chien(dpart,rhop,rho)
          case ('TRANCONG')
          cdmodel=trancong(dpart,rhop,rho)
          case ('DELLINO')
          cdmodel=dellino(dpart,rhop,rho)
          case ('HOLZSOMM')
          cdmodel=holzsomm(dpart,rhop,rho)
          case ('DIOGMELE')
          cdmodel=diogmele(dpart,rhop,rho)
          case ('DIOG2017')
          cdmodel=diog2017(dpart,rhop,rho)
          case ('DIOG2018')
          cdmodel=diog2018(dpart,rhop,rho)
          end select
      end function cdmodel
      
      function sphere(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,cdsphere,sphere,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
!     START OF ITERATIONS
      reold=restart
      cdold=(24.d0/reold)*(1.d0+0.15d0*(reold**0.687d0))+0.42d0/(1.d0+42500.d0/(reold**1.16))
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      cdnew=(24.d0/renew)*(1.d0+0.15d0*(renew**0.687d0))+0.42d0/(1.d0+42500.d0/(renew**1.16))
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      sphere=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function sphere

      function haidlev(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,haidlev,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: par                                                 ! Shape parameter considered in the law
      par=shpar1(ishape,jshape)
!     START OF ITERATIONS
      reold=restart
      cdold=(24.d0/reold)*(1.d0+(8.1716*exp(-4.0655d0*par))*(reold**(0.0964d0+0.5565d0*par))) &
     &+(73.69d0*reold*exp(-5.0748d0*par))/(reold+5.378*exp(6.2122*par))
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      cdnew=(24.d0/renew)*(1.d0+(8.1716*exp(-4.0655d0*par))*(renew**(0.0964d0+0.5565d0*par))) &
     &+(73.69d0*renew*exp(-5.0748d0*par))/(renew+5.378*exp(6.2122*par))
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      haidlev=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function haidlev

      function swamoj(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,swamoj,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: par                                                 ! Shape parameter considered in the law
      par=shpar1(ishape,jshape)
!      if(shpar.le.small) shpar=dshort(npart)/sqrt(dlong(npart)*dmed(npart))
!     START OF ITERATIONS
      reold=restart
      cdold=48.5d0/(((1.d0+4.5d0*par**0.35d0)**0.8d0)*reold**0.64d0)+ &
     &((reold/(reold+100.d0+100.d0*par))**0.32d0)*(1.d0/(par**18.d0+1.05d0*(par**0.8d0)))
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      cdnew=48.5d0/(((1.d0+4.5d0*par**0.35d0)**0.8d0)*renew**0.64d0)+ &
     &((renew/(renew+100.d0+100.d0*par))**0.32d0)*(1.d0/(par**18.d0+1.05d0*(par**0.8d0)))
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      swamoj=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function swamoj

      function ganser(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,ganser,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: k1,k2,par1,par2,par3                                  ! Shape parameter considered in the law
      logical :: par4
      par1=shpar1(ishape,jshape)
      par2=shpar2(ishape,jshape)
      par3=shpar3(ishape,jshape)
!      par4=shpar4(ishape,jshape)
      par4=shpar4(ishape)
      if(par4.eqv..TRUE.) then
      k1=(1.d0/3.d0+(2.d0/3.d0)*par1**(-0.5d0))**(-1.d0)
      else
      k1=((par3/(3.d0*par2))+(2.d0/3.d0)*par1**(-0.5d0))**(-1.d0)
      endif
      k2=10.d0**(1.8148d0*(-dlog10(par1))**0.5743d0)
!     START OF ITERATIONS
      reold=restart
      cdold=((24.d0/(reold*k1*k2))*(1.d0+0.1118d0*(reold*k1*k2)**0.6567d0)+0.4305d0/(1.d0+3305.d0/(reold*k1*k2)))*k2
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      cdnew=((24.d0/(renew*k1*k2))*(1.d0+0.1118d0*(renew*k1*k2)**0.6567d0)+0.4305d0/(1.d0+3305.d0/(renew*k1*k2)))*k2
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      ganser=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function ganser

      function chien(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,chien,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: par                                                 ! Shape parameter considered in the law
      par=shpar1(ishape,jshape)
!     START OF ITERATIONS
      reold=restart
      cdold=30.d0/reold+67.289d0*exp(-5.03d0*par)
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      cdnew=30.d0/renew+67.289d0*exp(-5.03d0*par)
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      chien=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function chien

      function trancong(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,trancong,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: par1,par2                                         ! Shape parameter considered in the law
      par1=shpar1(ishape,jshape)
      par2=shpar2(ishape,jshape)
!      if(shpar2.le.small) shpar2=sureqsphd(npart)/voleqsphd(npart)
!     START OF ITERATIONS
      reold=restart
      cdold=(24.d0/reold)*par2*(1.d0+(0.15d0/sqrt(par1))*(par2*reold)**0.687d0)+ &
     &(0.42d0*par2**2.d0)/(sqrt(par1)*(1.d0+4.25d4*(par2*reold)**(-1.16d0)))
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      cdnew=(24.d0/renew)*par2*(1.d0+(0.15d0/sqrt(par1))*(par2*renew)**0.687d0)+ &
     &(0.42d0*par2**2.d0)/(sqrt(par1)*(1.d0+4.25d4*(par2*renew)**(-1.16d0)))
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      trancong=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function trancong

      function dellino(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho                                             ! User's input
      real(dp) :: cdnew,dellino
      real(dp) :: par                                                 ! Shape parameter considered in the law
      par=shpar1(ishape,jshape)
      dellino=(0.69d0*g*(dpart**3)*rho*(1.33d0*rhop-1.33d0*rho))/((mu**2)*(((g*(par**1.6d0)*&
     &(dpart**3)*rho*(rhop-rho))/(mu**2))**1.0412))
      cdnew=dellino
      wt=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      end function dellino
      
      function holzsomm(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,holzsomm,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: par1,par2,par3                                        ! Shape parameter considered in the law
      par1=shpar1(ishape,jshape)
      par2=shpar2(ishape,jshape)
      par3=shpar3(ishape,jshape)
!     START OF ITERATIONS
      reold=restart
      cdold=(8.d0/reold)*(1.d0/sqrt(par2))+(16.d0/reold)*(1.d0/sqrt(par1))+&
     &(3.d0/reold)*(1.d0/(par1**(3.d0/4.d0)))+(0.4210d0**(0.4d0*((-log(par1))**0.2d0)))*(1.d0/par3)
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      cdnew=(8.d0/renew)*(1.d0/sqrt(par2))+(16.d0/renew)*(1.d0/sqrt(par1))+&
     &(3.d0/renew)*(1.d0/(par1**(3.d0/4.d0)))+(0.4210d0**(0.4d0*((-log(par1))**0.2d0)))*(1.d0/par3)
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      holzsomm=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function holzsomm
      
      function diogmele(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: exp
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,cdsphere,diogmele,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: par                                                 ! Shape parameter considered in the law
      integer :: it,it2
      par=shpar1(ishape,jshape)
!     START OF ITERATIONS
      reold=restart
      if(restart.le.50.d0) then
      exp=exp1
      else
      exp=exp2
      endif
      exp=reold**exp
      cdsphere=(24.d0/reold)*(1.d0+0.15d0*(reold**0.687d0))+0.42d0/(1.d0+42500.d0/(reold**1.16))
      cdold=(cdsphere*((reold/acd)**(1.d0/bcd)))/((reold**2)*(par**exp))
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
      do it=1,100
      if(reold.le.50.d0) then
      exp=exp1
      else
      exp=exp2
      endif
      exp=reold**exp
      cdsphere=(24.d0/reold)*(1.d0+0.15d0*(reold**0.687d0))+0.42d0/(1.d0+42500.d0/(reold**1.16))
      cdnew=(cdsphere*((reold/acd)**(1.d0/bcd)))/((reold**2)*(par**exp))
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      renew=(rho*wtnew*dpart)/mu
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      diogmele=cdnew
      wt=wtnew
      exit
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      endif
      enddo
      if(it.ge.100) then      ! Apply the relationship that is valid in Re=10-100
      re=reold
      diogmele=cdold
      wt=wtold
      reold=50.d0
      else
      endif
      if(isnan(rescd)) then                                             ! Apply Stokes
      wt=(2*(rhop-rho)*g*((dpart/2.d0)**2))/(9.d0*mu)
      re=(rho*wt*dpart)/mu
      diogmele=(24.d0/re)*(1.d0+0.15d0*(re**0.687d0))+0.42d0/(1.d0+42500.d0/(re**1.16))
!      write(100,*)'   stokes',ishape,jshape,rho,g,mu,rhop,dpart,re,diogmele,wt
      else
!      write(100,*)'no stokes',ishape,jshape,rho,g,mu,rhop,dpart,re,diogmele,wt
      endif
      end function diogmele
      
      function diog2017(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,cdsphere,diog2017,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: par1                                                  ! Shape parameter considered in the law
      real(dp) :: a,b,exp_1,exp_2
      par1=shpar1(ishape,jshape)
      if(fractal(ishape)) then
      a=a_fd
      b=b_fd
      exp_1=exp1_fd
!      exp_2=exp2_fd
      else
      a=a_sph
      b=b_sph
      exp_1=exp1_sph
!      exp_2=exp2_sph
      endif
!     START OF ITERATIONS
      reold=restart
      if(fractal(ishape)) then
      exp_2=reold**exp2_fd
      else
      exp_2=-(reold**exp2_sph)
      endif
      cdsphere=(24.d0/reold)*(1.d0+0.15d0*(reold**0.687d0))+0.42d0/(1.d0+42500.d0/(reold**1.16))
      cdold=(4.d0/3.d0)*(a*cdsphere*(((reold**exp_1)*(par1**exp_2))**b))/(reold**2.d0)
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      if(fractal(ishape)) then
      exp_2=renew**exp2_fd
      else
      exp_2=-(renew**exp2_sph)
      endif
      cdsphere=(24.d0/renew)*(1.d0+0.15d0*(renew**0.687d0))+0.42d0/(1.d0+42500.d0/(renew**1.16))
      cdnew=(4.d0/3.d0)*(a*cdsphere*(((renew**exp_1)*(par1**exp_2))**b))/(renew**2.d0)
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.reswt.le.tolcd) then
      re=renew
      diog2017=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function diog2017
      
      function diog2018(dpart,rhop,rho)
      USE inoutdata; USE nrtype
      implicit none
      real(dp) :: dpart,rhop,rho
      real(dp) :: renew,reold,re,cdnew,cdold,diog2018,wtnew,wtold
      real(dp) :: resre,rescd,reswt                                     ! Residuals
      real(dp) :: par                                                   ! Shape parameter considered in the law
      par=shpar1(ishape,jshape)
!     START OF ITERATIONS
      reold=restart
      cdold=(24.d0/reold)*(((1.d0-par)/reold+1.d0)**exp_a)+(24.d0/reold)*0.1806d0*(reold**0.6459d0)*(par**(-(reold**exp_b)))&
      +0.4251d0/(1.d0+6880.95d0/(reold*(par**exp_c)))
      wtold=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdold*rho))
   50 renew=(rho*wtold*dpart)/mu
      cdnew=(24.d0/renew)*(((1.d0-par)/renew+1.d0)**exp_a)+(24.d0/renew)*0.1806d0*(renew**0.6459d0)*(par**(-(renew**exp_b)))&
      +0.4251d0/(1.d0+6880.95d0/(renew*(par**exp_c)))
      wtnew=sqrt((4.d0*g*dpart*(rhop-rho))/(3.d0*cdnew*rho))
      resre=abs(renew-reold)
      rescd=abs(cdnew-cdold)
      reswt=abs(wtnew-wtold)
      if(resre.lt.tolcd.and.rescd.le.tolcd.and.reswt.le.tolcd) then
      re=renew
      diog2018=cdnew
      wt=wtnew
      else
      cdold=cdnew
      wtold=wtnew
      reold=renew
      goto 50
      endif
      end function diog2018

