      subroutine deposition
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
		
		FUNCTION zbrent(x1,x2)
                USE inoutdata; USE nrtype; USE nrutil, ONLY : nrerror
                IMPLICIT NONE
	        REAL(dp), INTENT(IN) :: x1,x2
         	REAL(dp) :: zbrent
         	END FUNCTION zbrent
	END INTERFACE
      integer :: i,j,k
      integer, dimension (6) :: ncltemp1                                ! To store temporary classes number
      integer, dimension(6) :: nclasses                                 ! To store final classes number (after merging and cut)
      integer, dimension(6) :: jpn,jend                                 ! To store the class index for which Pn<5 and for which Pn>0.8
      real(dp), dimension(6,30) :: d,ps                                 ! Particle size and weight fraction
      real(dp), dimension(6,30) :: wtp,cdp,pn,c,r               ! terminal velocity, drag coefficient, Rouse number, deposition probability, concentration, deposition rates
      real(dp), dimension(6,30) :: c_massive,ps_massive,wtp_massive,r_massive     ! Particle concentration, weight fraction, terminal velocity and deposition rate of the wash load part that contributes to the fine massive layer
      real(dp), dimension(6,30) :: b, csi, wstar, qb, taur              ! Wilcock and Crowe parameters
      real(dp) :: tau_rsm, tau_rsms, ref_a, ref_w                       ! Wilcock and Crowe parameters
      real(dp), dimension(4,5) :: store_sol                             ! To store PYFLOW solutions
      real(dp), dimension(6,30,5) :: d_dep,ps_dep,rhos_dep,wtp_dep,pn_dep,ps_dep_norm,shpar1_dep,shpar2_dep,shpar3_dep ! Input data stored for each solutions after reworking (e.g. cutting classes)
      real(dp) ::  dweight,rhosweight,shpar1weight,shpar2weight,shpar3weight
      logical :: noskip2,noskip3,switch
      real(dp) ::  rhospst,ctottemp,pstot                               ! Temporary macroscopic density, temporary total concentration, total weight fraction
      real(dp) :: rhoflow_wash                                          ! Flow density when wash load was freely settling
      real(dp), dimension(5):: ps_massive_tot,ps_lam_tot,ps_wash_tot,rhospst_susp,rhospst_wash,rhospst_massive
      real(dp) :: weight_ave_rhos                                           ! Weighted average particle density
      real(dp), dimension(6,30) :: r_scaled,c_max,c_ratio,ps_scaled               ! Deposition rate scaled with ps to sum up to R_meas, particle volumetric concentration at sedimentation
      real(dp) :: r_corrector                      ! Correction parameter of Dellino 2018 model

      open(60,file='detailed_dep_results.dat')
      open(61,file='deposition_summary.dat')
!      open(70,file='C_check.dat')
!      open(71,file='C_check_final.dat')
      write(*,*)'*******DEPOSITION RATES AND TIMES CALCULATION********'
      write(52,*)'*******DEPOSITION RATES AND TIMES CALCULATION********'
      write(60,*)'*******DEPOSITION RATES AND TIMES CALCULATION********'
      write(61,*)'Solution Component Class     dp (m)     ps(-)  wt(m s^-1)     Pn(-)       C(-) rhos(kg m^-3) R(kg m^-1 s^-2)'
      noskip2=.TRUE.                                                    !Flag to skip the calculation for merging classes for the shape parameters 2 and 3
      noskip3=.TRUE.

      r_corrector = 0.d0
      if(dep_model.eq.'DELLINO_2018') r_corrector = -0.01d0
! Normalize layer thickness to the considered components total weight fractions
      if(psreadtot.le.1.d0) zlam=zlam*psreadtot                         ! Avoid normalization if wt_tot>1 (100%), although this should not be the case
! Closure to 100% of the considered components and definitions of shape factors matrices
      do i=1,ncomp
!      nclass_dep(i)=nclass(i)
!         do j=1,nclass(i)      !old
         do j=1,nclass_dep(i)  !new
         psread(i,j)=psread(i,j)/psreadtot
         rhosread(i,j)=rhos(i,j)
         select case (cdlaw(i))
                case ('HAIDLEV')
                shpar1read(i,j)=sphericity(i,j)
                case ('SWAMOJ')
                shpar1read(i,j)=corey(i,j)
                case ('GANSER')
                shpar1read(i,j)=sphericity(i,j)
                shpar2read(i,j)=voleqsphd(i,j)
                shpar3read(i,j)=circeqard(i,j)
                case ('CHIEN')
                shpar1read(i,j)=sphericity(i,j)
                case ('TRANCONG')
                shpar1read(i,j)=circularity(i,j)
                shpar2read(i,j)=flatratio(i,j)
                case ('DELLINO')
                shpar1read(i,j)=shapefact(i,j)
                case ('HOLZSOMM')
                shpar1read(i,j)=sphericity(i,j)
                shpar2read(i,j)=longspher(i,j)
                shpar3read(i,j)=crossspher(i,j)
                case ('DIOGMELE')
                shpar1read(i,j)=shapefact(i,j)
                case ('DIOG2017')
                     if(fractal(i)) then
                     shpar1read(i,j)=fractdim(i,j)
                     else
                     shpar1read(i,j)=sphericity(i,j)
                     endif
                case ('DIOG2018')
                shpar1read(i,j)=shapefact(i,j)
         end select
         enddo
      enddo

! Delete classes with wt%=0
      do i=1,ncomp
            j=1
  192       if(psread(i,j).lt.1.d-8) then
            nclass_dep(i)=nclass_dep(i)-1
               do k=j,nclass_dep(i)
               psread(i,k)=psread(i,k+1)
               rhosread(i,k)=rhosread(i,k+1)
               dread(i,k)=dread(i,k+1)
               shpar1read(i,k)=shpar1read(i,k+1)
               if(noskip2) shpar2read(i,k)=shpar2read(i,k+1)
               if(noskip3) shpar3read(i,k)=shpar3read(i,k+1)
               enddo
            j=1
            goto 192
            else
            j=j+1
               if(j.gt.nclass_dep(i)) then
               cycle
               else
               goto 192
               endif
            endif
      enddo

      do i=1,ncomp
      if(merge_classes(i)) then
      
!     Merge grainsizes
      ncltemp1(i)=1
      dweight=0.d0
      rhosweight=0.d0
      shpar1weight=0.d0
      shpar2weight=0.d0
      shpar3weight=0.d0
         do j=1,nclass_dep(i)-1
         ps(i,ncltemp1(i))=ps(i,ncltemp1(i))+psread(i,j)
         dweight=dweight+psread(i,j)*dread(i,j)
         rhosweight=rhosweight+psread(i,j)*rhosread(i,j)
         shpar1weight=shpar1weight+psread(i,j)*shpar1read(i,j)
                if(shpar2(i,j).ne.undefined) then
                shpar2weight=shpar2weight+psread(i,j)*shpar2read(i,j)
                else
                noskip2=.FALSE.
                endif
                if(shpar3(i,j).ne.undefined) then
                shpar3weight=shpar3weight+psread(i,j)*shpar3read(i,j)
                else
                noskip3=.FALSE.
                endif
         d(i,ncltemp1(i))=dweight/ps(i,ncltemp1(i))
         rhos(i,ncltemp1(i))=rhosweight/ps(i,ncltemp1(i))
         shpar1(i,ncltemp1(i))=shpar1weight/ps(i,ncltemp1(i))
         if(noskip2) shpar2(i,ncltemp1(i))=shpar2weight/ps(i,ncltemp1(i))
         if(noskip3) shpar3(i,ncltemp1(i))=shpar3weight/ps(i,ncltemp1(i))
         if(ps(i,ncltemp1(i)).lt.sensmerge(i)) then
         cycle
         else
         ncltemp1(i)=ncltemp1(i)+1
         dweight=0.d0
         rhosweight=0.d0
         shpar1weight=0.d0
         shpar2weight=0.d0
         shpar3weight=0.d0
         endif
         enddo
      if(psread(i,nclass_dep(i)).lt.sensmerge(i)) then
      ps(i,ncltemp1(i))=ps(i,ncltemp1(i))+psread(i,nclass_dep(i))
      dweight=dweight+psread(i,nclass_dep(i))*dread(i,nclass_dep(i))
      rhosweight=rhosweight+psread(i,nclass_dep(i))*rhosread(i,nclass_dep(i))
      shpar1weight=shpar1weight+psread(i,nclass_dep(i))*shpar1read(i,nclass_dep(i))
      if(noskip2) shpar2weight=shpar2weight+psread(i,nclass_dep(i))*shpar2read(i,nclass_dep(i))
      if(noskip3) shpar3weight=shpar3weight+psread(i,nclass_dep(i))*shpar3read(i,nclass_dep(i))
      d(i,ncltemp1(i))=dweight/ps(i,ncltemp1(i))
      rhos(i,ncltemp1(i))=rhosweight/ps(i,ncltemp1(i))
      shpar1(i,ncltemp1(i))=shpar1weight/ps(i,ncltemp1(i))
      if(noskip2) shpar2(i,ncltemp1(i))=shpar2weight/ps(i,ncltemp1(i))
      if(noskip3) shpar3(i,ncltemp1(i))=shpar3weight/ps(i,ncltemp1(i))
   78   if(ps(i,ncltemp1(i)).lt.sensmerge(i)) then
        ncltemp1(i)=ncltemp1(i)-1
        dweight=dweight+ps(i,ncltemp1(i))*d(i,ncltemp1(i))
        rhosweight=rhosweight+ps(i,ncltemp1(i))*rhos(i,ncltemp1(i))
        shpar1weight=shpar1weight+ps(i,ncltemp1(i))*shpar1(i,ncltemp1(i))
        if(noskip2) shpar2weight=shpar2weight+ps(i,nclass_dep(i))*shpar2(i,nclass_dep(i))
        if(noskip3) shpar3weight=shpar3weight+ps(i,nclass_dep(i))*shpar3(i,nclass_dep(i))
        ps(i,ncltemp1(i))=ps(i,ncltemp1(i))+ps(i,ncltemp1(i)+1)
        d(i,ncltemp1(i))=dweight/ps(i,ncltemp1(i))
        rhos(i,ncltemp1(i))=rhosweight/ps(i,ncltemp1(i))
        shpar1(i,ncltemp1(i))=shpar1weight/ps(i,ncltemp1(i))
        if(noskip2) shpar2(i,ncltemp1(i))=shpar2weight/ps(i,ncltemp1(i))
        if(noskip3) shpar3(i,ncltemp1(i))=shpar3weight/ps(i,ncltemp1(i))
        goto 78
        else
        endif
      else
      ncltemp1(i)=ncltemp1(i)+1
      ps(i,ncltemp1(i))=ps(i,ncltemp1(i))+psread(i,nclass_dep(i))
      d(i,ncltemp1(i))=dread(i,nclass_dep(i))
      rhos(i,ncltemp1(i))=rhos(i,nclass_dep(i))
      shpar1(i,ncltemp1(i))=shpar1(i,nclass_dep(i))
      if(noskip2) shpar2(i,ncltemp1(i))=shpar2(i,nclass_dep(i))
      if(noskip3) shpar3(i,ncltemp1(i))=shpar3(i,nclass_dep(i))
      endif

      else
!     Not merging classes
      ncltemp1(i)=nclass_dep(i)
      do j=1,ncltemp1(i)
      ps(i,j)=psread(i,j)
      d(i,j)=dread(i,j)
      rhos(i,j)=rhosread(i,j)
      shpar1(i,j)=shpar1read(i,j)
      if(noskip2) shpar2(i,j)=shpar2read(i,j)
      if(noskip3) shpar3(i,j)=shpar3read(i,j)
      enddo
      endif
      enddo

!     Store data from the three flow solutions
      if(only_deprates) then
         do i=1,n_solutions
         store_sol(1,i)=rho_flow(i)
         store_sol(2,i)=ush_flow(i)
         store_sol(3,i)=ztot_flow(i)
         store_sol(4,i)=pns_flow(i)
         enddo
      kmax=n_solutions
      else
      store_sol(1,1)=dennrm
      store_sol(1,2)=denmin
      store_sol(1,3)=denmax
      store_sol(2,1)=ushavg
      store_sol(2,2)=ushmax
      store_sol(2,3)=ushmin
      store_sol(3,1)=ztavg
      store_sol(3,2)=ztmax
      store_sol(3,3)=ztmin
      store_sol(4,1)=pnsavg
      store_sol(4,2)=pnsmin
      store_sol(4,3)=pnsmax
      kmax=3
      n_solutions=kmax
      endif

!     Start cycle for the three deposition solutions
      do k=1,kmax
      rhoflow=store_sol(1,k)
      ushearflow=store_sol(2,k)
      ztotflow=store_sol(3,k)
      pnsflow=store_sol(4,k)
      zlam_final(k)=zlam
      if(.not.only_deprates) then
        if(k.eq.1) then
        write(*,*)'### AVERAGE SOLUTION ###'
        write(52,*)'### AVERAGE SOLUTION ###'
        write(60,*)'### AVERAGE SOLUTION ###'
        else
            if(k.eq.2) then
            write(*,*)'### MAXIMUM SOLUTION ###'
            write(52,*)'### MAXIMUM SOLUTION ###'
            write(60,*)'### MAXIMUM SOLUTION ###'
            else
            write(*,*)'### MINIMUM SOLUTION ###'
            write(52,*)'### MINIMUM SOLUTION ###'
            write(60,*)'### MINIMUM SOLUTION ###'
            endif
        endif
      else
      write(*,92)k
      write(52,92)k
      write(60,92)k
      endif

!     Store classes for which Pn>5
      pstot=0.d0
      do i=1,ncomp
      nclasses(i)=1
      jpn(i)=0
      ishape=i                                                        ! Identifies the component for Cd calcualation
      switch=.FALSE.
         do j=1,ncltemp1(i)-1
         jshape=j
         cdp(i,j)=cdmodel(d(i,j),rhos(i,j),rhoflow)
         wtp(i,j)=wt
         pn(i,j)=wtp(i,j)/(kvk*ushearflow)
                if(.not.switch.and.pn(i,j).ge.5.d0) then
                jpn(i)=j-1
                cycle
                else
                switch=.TRUE.
                nclasses(i)=nclasses(i)+1
                endif
         enddo
      j=ncltemp1(i)
      jshape=j
      cdp(i,j)=cdmodel(d(i,j),rhos(i,j),rhoflow)
      wtp(i,j)=wt
      pn(i,j)=wtp(i,j)/(kvk*ushearflow)
      nclasses(i)=nclasses(i)+1
      if(pn(i,1).lt.5.d0.or..not.pn_cut) then
          nclasses(i)=ncltemp1(i)                                           ! Correction in case there is no need to cut classes
               do j=1,nclasses(i)
               d_dep(i,j,k)=d(i,j)
               ps_dep(i,j,k)=ps(i,j)
               pstot=pstot+ps_dep(i,j,k)
               rhos_dep(i,j,k)=rhos(i,j)
               wtp_dep(i,j,k)=wtp(i,j)
               pn_dep(i,j,k)=pn(i,j)
               shpar1_dep(i,j,k)=shpar1(i,j)
               if(noskip2) shpar2_dep(i,j,k)=shpar2(i,j)
               if(noskip3) shpar3_dep(i,j,k)=shpar3(i,j)
               enddo
          else
               do j=1,nclasses(i)
               d_dep(i,j,k)=d(i,j+jpn(i))
               ps_dep(i,j,k)=ps(i,j+jpn(i))
               pstot=pstot+ps_dep(i,j,k)
               rhos_dep(i,j,k)=rhos(i,j+jpn(i))
               wtp_dep(i,j,k)=wtp(i,j+jpn(i))
               pn_dep(i,j,k)=pn(i,j+jpn(i))
               shpar1_dep(i,j,k)=shpar1(i,j+jpn(i))
               if(noskip2) shpar2_dep(i,j,k)=shpar2(i,j+jpn(i))
               if(noskip3) shpar3_dep(i,j,k)=shpar3(i,j+jpn(i))
               enddo
          endif
      enddo

      zlam_final(k)=zlam_final(k)*pstot
!     Re-normalize classes weights to close to 100%                     ! TO BE CHECKED!
      do i=1,ncomp
          do j=1,nclasses(i)
          ps_dep(i,j,k)=ps_dep(i,j,k)/pstot
          enddo
      enddo
!     Start calculation

      rtot_susp(k)=0.d0
      rtot_massive(k)=0.d0
      ctot_susp(k)=0.d0
      ctot_massive(k)=0.d0
      ps_massive_tot(k)=0.d0
      ps_lam_tot(k)=0.d0
      rhospst_susp(k)=0.d0
      rhospst_massive(k)=0.d0
      ctottemp=0.d0
      rhoflow_wash=0.d0

!     Calculate average concentration in the flow
      weight_ave_rhos=0.d0
      do i=1,ncomp
         do j=1,nclasses(i)
         weight_ave_rhos=weight_ave_rhos+rhos_dep(i,j,k)*ps_dep(i,j,k)
         enddo
      enddo
      ctotflow=(rhoflow-dengas)/(weight_ave_rhos-dengas)
      ctot_flow(k)=ctotflow

      write(52,*)'Imput data summary'
      write(60,*)'Imput data summary'
      write(*,*)'Imput data summary'
      write(52,89)pstot,zlam_final(k),rhoflow,ushearflow,ztotflow
      write(60,89)pstot,zlam_final(k),rhoflow,ushearflow,ztotflow
      write(*,89)pstot,zlam_final(k),rhoflow,ushearflow,ztotflow

!     Just to calculate jend in case I decide to keep the calculation of the fine massive layer
      do i=1,ncomp
      jend(i)=0
         do j=1,nclasses(i)
                if(pn_dep(i,j,k).ge.0.8d0) then
                jend(i)=j
                else
                exit
                endif
         enddo
      enddo
      
      do i=1,ncomp
         do j=1,nclasses(i)
         ps_lam_tot(k)=ps_lam_tot(k)+ps_dep(i,j,k)
         enddo
      enddo

!     Calculate C from C_max
      call c_calc(ps_dep,rhos_dep,wtp_dep,pn_dep,nclasses,k,c)
      
!     Layer thicknesses
      zlam_susp(k)=zlam_final(k)*ps_lam_tot(k)

!     Particle volumetric concentration of each class in the wash load forming the fine massive layer, flow density after laminae deposition
!     Assuming that the fines in the massive layer have the same componentry and grainsize distribution of the fines trapped in the laminated
!     layer, I calculate the equivalent grainsize distribution in the fine massive layer. In a future ideal case we will have direct measurements
!     for this layer
      if(zlam_massive.eq.undefined) goto 111
      do i=1,ncomp
         do j=jend(i)+1,nclasses(i)
         ps_wash_tot(k)=ps_wash_tot(k)+ps_dep(i,j,k)
         enddo
      enddo
      do i=1,ncomp
         do j=jend(i)+1,nclasses(i)
         ps_massive(i,j)=ps_dep(i,j,k)/ps_wash_tot(k)
         ps_massive_tot(k)=ps_massive_tot(k)+ps_massive(i,j)
         c_massive(i,j)=(ps_massive(i,j)*zlam_massive)/ztotflow                     ! Assuming the flow is still that thick...
         ctot_massive(k)=ctot_massive(k)+c_massive(i,j)
         rhoflow_wash=rhoflow_wash+c_massive(i,j)*rhos_dep(i,j,k)
         enddo
      enddo
      rhoflow_wash=rhoflow_wash+(1.d0-ctot_massive(k))*1.22d0
!     Terminal velocity calculation of particles in the washload forming the fine massive layer
      do i=1,ncomp
      ishape=i
         do j=jend(i)+1,nclasses(i)
         jshape=j
         cdp(i,j)=cdmodel(d_dep(i,j,k),rhos_dep(i,j,k),rhoflow_wash)
         wtp_massive(i,j)=wt
         enddo
      enddo

!     Macroscopic density in the massive layer
      do i=1,ncomp
         do j=jend(i)+1,nclasses(i)
         rhospst_massive(k)=rhospst_massive(k)+rhos_dep(i,j,k)*ps_massive(i,j)
         enddo
      enddo

!     Deposition rates and times in the massive layer
      do i=1,ncomp
         do j=jend(i)+1,nclasses(i)
         r_massive(i,j)=rhos_dep(i,j,k)*wtp_massive(i,j)*c_massive(i,j)
         rtot_massive(k)=rtot_massive(k)+r_massive(i,j)
         enddo
      enddo

!     Macroscopic density in the laminated part: suspension component amd wash load
  111 do i=1,ncomp
         do j=1,nclasses(i)
         rhospst_susp(k)=rhospst_susp(k)+rhos_dep(i,j,k)*ps_dep(i,j,k)
         enddo
      enddo
      
!     Deposition rates and times in the laminated part from turbulent suspension
      do i=1,ncomp
         do j=1,nclasses(i)
         r(i,j)=rhos_dep(i,j,k)*wtp_dep(i,j,k)*c(i,j)+r_corrector
         rtot_susp(k)=rtot_susp(k)+r(i,j)
         ctot_dep(k)=ctot_dep(k)+c(i,j)
         enddo
      ctot_susp(k)=ctot_flow(k)-ctot_dep(k)
      enddo

!     Deposition times
      if(rtot_susp(k).gt.0.d0) tdep_susp(k)=(zlam_susp(k)*c0*rhospst_susp(k))/rtot_susp(k)                   !Deposition time of the laminated part from turbulent suspension
      if(rtot_massive(k).gt.0.d0) tdep_massive(k)=(zlam_massive*c0*rhospst_massive(k))/rtot_massive(k)       !Deposition time of the massive part
      tdep(k)=tdep_susp(k)+tdep_massive(k)
      rtot(k)=rtot_susp(k)+rtot_massive(k)
      ctot(k)=ctot_susp(k)+ctot_massive(k)
      !ar(k) = rtot_susp(k)/rhospst_susp(k) !Accretion rate of the laminated part
      ar(k) = rtot(k)/rhospst_susp(k) !Total accretion rate
      ref_a = 1.d0
      ref_w = 1.d0
      srw(k) = ar(k) * (ref_a/ref_w)
!     Bedload transportation rate (Dellino et al. 2020 after Wilcock & Crowe)
      tau_rsms = 0.75d0*0.015d0 !To be verified
      tau_rsm = tau_rsms * (dep_median*g*(rhos_median - rhoflow))
      qtot(k) = 0.d0
      do i=1,ncomp
         do j=1,nclasses(i)
         b(i,j) = 0.67d0 / (1.d0 + exp(1.5d0 - d_dep(i,j,k)/dep_median))
         taur(i,j) = tau_rsm * ((d_dep(i,j,k)/dep_median)**b(i,j))
         csi(i,j) = (rhoflow * ushearflow ** 2)/taur(i,j)
         if(csi(i,j).lt.1.35d0) then
            wstar(i,j) = 0.002d0 * csi(i,j) ** 7.5d0
         else
            wstar(i,j) = 14.d0*(1.d0-0.894d0/(csi(i,j)**0.5d0))
         endif
         qb(i,j) = (wstar(i,j)*ps_dep(i,j,k)*ushearflow**3.d0)/((rhos_median/rhoflow - 1.d0)*g)
         qtot(k) = qtot(k) + qb(i,j)
         enddo
      enddo
      sqratio(k) = srw(k)/qtot(k)

      do i=1,ncomp
      write(60,94)i
      write(*,94)i
           write(60,*)'   SUSPENSION LOAD'
           write(*,*)'   SUSPENSION LOAD'
           do j=1,nclasses(i)
           write(60,93)j,r(i,j),c(i,j),wtp_dep(i,j,k),pn_dep(i,j,k)
           write(*,93)j,r(i,j),c(i,j),wtp_dep(i,j,k),pn_dep(i,j,k)
           enddo
           if(zlam_massive.eq.undefined) cycle
           write(60,*)'   WASH LOAD CONTRIBUTING TO THE FINE MASSIVE LAYER'
           write(*,*)'   WASH LOAD CONTRIBUTING TO THE FINE MASSIVE LAYER'
           do j=jend(i)+1,nclasses(i)
           write(60,96)j,r_massive(i,j),c_massive(i,j),wtp_massive(i,j)
           write(*,96)j,r_massive(i,j),c_massive(i,j),wtp_massive(i,j)
           enddo
      enddo
      do i=1,ncomp
           do j=1,nclasses(i)
           write(61,97)k,i,j,d_dep(i,j,k),ps_dep(i,j,k),wtp_dep(i,j,k),pn_dep(i,j,k),c(i,j),rhos_dep(i,j,k),r(i,j)
           enddo
      enddo
      write(61,*)''
      enddo
!     Store average, maximum and minimum solution for building probability functions of deposition rate and time
      if(n_solutions.ge.3) call sort_rates_times

   89 format('Total weight fraction = ',f5.3,/,'Normalized layer thickness = ',f5.3,/,&
     &'Flow density (kg/m^3) = ',f7.3,/,'Flow shear velocity (m/s) = ',f6.3,/,&
     &'Flow total thickness (m) = 'f8.3,/)
   90 format(/,'Component ',i2,/)
   91 format(i2,2x,e8.3,10x,f6.4,10x,f7.2,17x,f6.3,7x,f6.3,/)
   92 format(/,'### SOLUTION N ',i1,'###',/)
   93 format('Class n. = ',i3,/,'Rate of deposition (kg m^-2 s^-1) = ',e10.3,/, &
      'Particle concentration = ',e10.3,/, &
      'Particle terminal velocity (m s^-1) = ',f10.3,/, &
      'Particle Rouse number = ',f7.3,/)
   95 format('Class n. = ',i3,/,'Rate of deposition (kg m^-2 s^-1) = ',e10.3,/, &
      'Particle concentration = ',e10.3,/, &
      'Particle terminal velocity (m s^-1) = ',f10.3,/, &
      'Particle Rouse number = ',f7.3,/)
   96 format('Class n. = ',i3,/,'Rate of deposition (kg m^-2 s^-1) = ',e10.3,/, &
      'Particle concentration = ',e10.3,/, &
      'Particle terminal velocity (m s^-1) = ',f10.3,/)
   94 format('*** Component ',i2,' ***')
   97 format(6x,i2,9x,i2,2x,i4,2x,e9.3,2x,e9.3,2x,f8.4,3x,f8.4,2x,e9.3,6x,f8.3,7x,e9.3)
      end subroutine deposition
      

      subroutine sort_rates_times
      USE inoutdata; USE nrtype
      implicit none
	INTERFACE
                 SUBROUTINE sort_pick(arr)
                 USE nrtype
                 IMPLICIT NONE
                 REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
                 END SUBROUTINE sort_pick
	END INTERFACE
      real(dp), dimension(n_solutions):: rtot_new,tdep_new
      integer :: i
!     Resize arrays
      do i=1,n_solutions
      rtot_new(i)=rtot(i)
      tdep_new(i)=tdep(i)
      enddo
!     Sort Rtot
      call sort_pick(rtot_new)
      rtot_max=rtot_new(n_solutions)
      rtot_min=rtot_new(1)
      if(n_solutions.eq.3) then
      rtot_avg=rtot_new(2)
      else
      rtot_avg=sum(rtot_new)/n_solutions
      endif
!     Sort tdep
      call sort_pick(tdep_new)
      tdep_max=tdep_new(n_solutions)
      tdep_min=tdep_new(1)
      if(n_solutions.eq.3) then
      tdep_avg=tdep_new(2)
      else
      tdep_avg=sum(tdep_new)/n_solutions
      endif
      end subroutine sort_rates_times
      
      SUBROUTINE sort_pick(arr)
      USE nrtype
      IMPLICIT NONE
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER(I4B) :: i,j,n
      REAL(DP) :: a
      n=size(arr)
      do j=2,n
      a=arr(j)
              do i=j-1,1,-1
              if (arr(i) <= a) exit
              arr(i+1)=arr(i)
              end do
              arr(i+1)=a
      end do
      END SUBROUTINE sort_pick

      subroutine c_calc(ps_dep,rhos_dep,wtp_dep,pn_dep,nclasses,k,c)
      USE inoutdata; USE nrtype
      implicit none
      integer, dimension(6) :: nclasses
      real(dp), dimension(6,30) :: c
      real(dp), dimension(6,30,5) :: ps_dep,rhos_dep,wtp_dep,pn_dep
      real(dp), dimension(:), allocatable :: ps_vec,rhos_vec,wtp_vec,c_vec,pn_vec,numerator,denominator,sums
      real(dp), dimension(:), allocatable :: c_trans,gamma_proxy
      real(dp) :: rhos_prod,wtp_prod,pn_avg,pn_star,den2,ctottemp
      integer :: nmax,i,j,k,l

!     Calculate the vector length
      nmax=0
      do i=1,ncomp
      nmax=nmax+nclasses(i)
      enddo
      allocate(ps_vec(nmax),rhos_vec(nmax),wtp_vec(nmax),c_vec(nmax),pn_vec(nmax),numerator(nmax),sums(nmax),denominator(nmax))
      allocate(c_trans(nmax),gamma_proxy(nmax))

!     Create vectors of variables
      l=1
      do i=1,ncomp
         do j=1,nclasses(i)
         ps_vec(l)=ps_dep(i,j,k)
         rhos_vec(l)=rhos_dep(i,j,k)
         wtp_vec(l)=wtp_dep(i,j,k)
         pn_vec(l)=pn_dep(i,j,k)
         c_vec(l)=0.d0
         l=l+1
         enddo
      enddo


      !The model presented in Dioguardi et al. 2018 is solved
      do l=1,nmax
      rhos_prod=1.d0
      wtp_prod=1.d0
         do j=1,l-1
         rhos_prod=rhos_prod*rhos_vec(j)
         wtp_prod=wtp_prod*wtp_vec(j)
         enddo
         do j=l+1,nmax
         rhos_prod=rhos_prod*rhos_vec(j)
         wtp_prod=wtp_prod*wtp_vec(j)
         enddo
      numerator(l)=ctotflow*ps_vec(l)*rhos_prod*wtp_prod
      enddo
      do l=1,nmax
      rhos_prod=1.d0
      wtp_prod=1.d0
         do j=1,l-1
         rhos_prod=rhos_prod*rhos_vec(j)
         wtp_prod=wtp_prod*wtp_vec(j)
         enddo
         do j=l+1,nmax
         rhos_prod=rhos_prod*rhos_vec(j)
         wtp_prod=wtp_prod*wtp_vec(j)
         enddo
      sums(l)=ps_vec(l)*rhos_prod*wtp_prod
      enddo
      do l=1,nmax
      denominator(l)=0.d0
         do j=1,nmax
         denominator(l)=denominator(l)+sums(j)
         enddo
      enddo
      do l=1,nmax
      c_vec(l)=numerator(l)/denominator(l)
      enddo

      if(dep_model.eq.'DELLINO_2018') then
!     The model of Dellino et al. (2018) is solved
      pn_avg=0.d0
      den2=0.d0
      do i=1,nmax
      den2=den2+ps_vec(i)/rhos_vec(i)
      enddo
      do i=1,nmax
      c_trans(i)=((ps_vec(i)/rhos_vec(i))/den2)*ctotflow
      pn_avg=pn_avg+pn_vec(i)*c_trans(i)
      enddo
      pn_avg=pn_avg/ctotflow
      pn_star=pn_avg/pnsflow
      ctottemp=0
      do i=1,nmax
      gamma_proxy(i)=c_trans(i)/(10.065d0*(pn_vec(i)/pn_star)+0.1579d0)
      enddo
      do i=1,nmax-1
      c_vec(i)=0.7d0*gamma_proxy(i)+0.3d0*gamma_proxy(i+1)
      ctottemp=ctottemp+c_vec(i)
      enddo
      c_vec(nmax)=0.7d0*gamma_proxy(nmax)
      endif

!     Redistribute c_vec among components and classes
      l=1
      do i=1,ncomp
         do j=1,nclasses(i)
         c(i,j)=c_vec(l)
         l=l+1
         enddo
      enddo
      deallocate(ps_vec,rhos_vec,wtp_vec,c_vec,pn_vec,numerator,sums,denominator)

      end subroutine c_calc
      
      

