      subroutine check_data
      use inoutdata; use nrtype
      implicit none
	INTERFACE
		character(len=20) function str(k)
		implicit none
		integer, intent(in) :: k
		end function str
         END INTERFACE
      integer :: i,j,ii,jj,jjstart,kk,ncomp_checks
      logical :: error_d1,error_d2,error_d3
      logical :: error,error_drag_final,error_distr,error_deprates            ! Flags to stop PYFLOW in case of errors in the input data
      logical, dimension(6) :: error_drag,error_grainsize
!      error_d1=.FALSE.
!      error_d2=.FALSE.
!      error_d3=.FALSE.
      if(deprates) then
      ncomp_checks=ncomp
      else
             if(model.eq.'TWOLAYERS') then
             ncomp_checks=1
             else
             ncomp_checks=2
             endif
      endif

      do i=1,ncomp_checks
      error_drag(i)=.FALSE.
      error_grainsize(i)=.FALSE.
      enddo
      error=.FALSE.
      error_drag_final=.FALSE.
      error_distr=.FALSE.
      error_deprates=.FALSE.
      write(*,*)''
      write(52,*)''
      if(mu.eq.UNDEFINED) then
      write(*,*)'FATAL ERROR! Command MU missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command MU missing in input.dat'
      write(52,*)''
      error=.TRUE.
      else
      endif

      if(dengas.eq.UNDEFINED) then
      write(*,*)'FATAL ERROR! Command DENGAS missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command DENGAS missing in input.dat'
      write(52,*)''
      error=.TRUE.
      else
      endif

      if(c0.eq.UNDEFINED) then
      write(*,*)'WARNING! Command C0 missing in input.dat'
      write(*,*)'Setting C0 = 0.75'
      c0=0.75d0
      else
      endif

      if(zlam.eq.UNDEFINED) then
      write(*,*)'FATAL ERROR! Command ZLAM missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command ZLAM missing in input.dat'
      write(52,*)''
      error=.TRUE.
      else
      endif
      if(only_deprates) then
          if(.not.deprates) then
          write(*,*)'FATAL ERROR! Command DEPRATES missing or set to .FALSE.'
          write(*,*)''
          write(52,*)'FATAL ERROR! Command DEPRATES missing or set to .FALSE.'
          write(52,*)''
          error=.TRUE.
          else
          endif
          if(n_solutions.eq.undefined_i) then
          write(*,*)'FATAL ERROR! Command N_SOLUTIONS missing in input.dat'
          write(*,*)''
          write(52,*)'FATAL ERROR! Command N_SOLUTIONS missing in input.dat'
          write(52,*)''
          error=.TRUE.
          else
              if(n_solutions.gt.5) then
              write(*,*)'FATAL ERROR! N_SOLUTIONS cannot be larger than 5'
              write(*,*)''
              write(52,*)'FATAL ERROR! N_SOLUTIONS cannot be larger than 5'
              write(52,*)''
              error=.TRUE.
              else
                  do kk=1,n_solutions
                     if(rho_flow(kk).eq.undefined) then
                     write(*,152)kk
                     write(52,152)kk
                     error=.TRUE.
                     else
                     endif
                     if(ztot_flow(kk).eq.undefined.and.zlam_massive.ne.undefined) then
                     write(*,153)kk
                     write(52,153)kk
                     error=.TRUE.
                     else
                     endif
                     if(ush_flow(kk).eq.undefined) then
                     write(*,155)kk
                     write(52,155)kk
                     error=.TRUE.
                     else
                     endif
                     if(pns_flow(kk).eq.undefined.and.dep_model.eq.'DELLINO_2018') then
                     write(*,157)kk
                     write(52,157)kk
                     error=.TRUE.
                     else
                     endif
                  enddo
              endif
          endif
          if(ncomp.eq.UNDEFINED_I) then
          write(*,*)'FATAL ERROR! Command NCOMP missing in input.dat'
          write(*,*)''
          write(52,*)'FATAL ERROR! Command NCOMP missing in input.dat'
          write(52,*)''
          error=.TRUE.
          else
              do jj=1,ncomp
              call check_distribution(jj,error_grainsize(i))
              call check_drag_input(jj,1,error_drag(i))
              enddo
          call check_deprates(error_deprates)
          endif
          goto 567
      else
      endif


      if(model.eq.UNDEFINED_C) then
      write(*,*)'FATAL ERROR! Solution method not defined in input.dat'
      write(*,*)'Please write solution method in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Solution method not defined in input.dat'
      write(52,*)'Please write solution method in input.dat'
      write(52,*)''
      error=.TRUE.

      else
        if(model.ne.'TWOLAYERS'.and.model.ne.'TWOCOMPONENTS') then
        write(*,*)'Wrong method choice'
        write(*,*)'Please write the correct method: TWOLAYERS or TWOCOMPONENTS'
        write(*,*)''
        write(52,*)'Wrong method choice'
        write(52,*)'Please write the correct method: TWOLAYERS or TWOCOMPONENTS'
        write(52,*)''
        error=.TRUE.
        else
        endif
      endif

      if(probt.eq.UNDEFINED) then
      write(*,*)'FATAL ERROR! Command PROBT missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command PROBT missing in input.dat'
      write(52,*)''
      error=.TRUE.
      else
      alfa=probt/2.d0
      endif

      if(ks.eq.UNDEFINED) then
      write(*,*)'FATAL ERROR! Command KS missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command KS missing in input.dat'
      write(52,*)''
      error=.TRUE.
      else
      endif
!     At least these data are needed. If model is two layers, they will change from 1 to 2
      ii=0                                                              ! Class identifier. 0 means median class in case grainsize distributions are not used
      jj=1                                                       ! Component identifier
      if(rholaw(jj).eq.undefined_c) then
      write(*,*)'FATAL ERROR! Command RHOLAW(1) missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command RHOLAW(1) missing in input.dat'
      write(52,*)''
      error=.TRUE.
      else
          if(rholaw(jj).ne.'POLLENA'.and.rholaw(jj).ne.'AVERNO2'.and.rholaw(jj).ne.'AMS'&
     &    .and.rholaw(jj).ne.'POMPEI'.and.rholaw(jj).ne.'SIAL_XX'.and.rholaw(jj).ne.'FEM_XX'&
     &    .and.rholaw(jj).ne.'LITHIC'.and.rholaw(jj).ne.'MERCATO'.and.rholaw(jj).ne.'ASTRONI'&
     &    .and.rholaw(jj).ne.'CUSTOM') then
          write(*,139)jj
          write(52,139)jj
          error=.TRUE.
          else
              if(rhos(jj,0).eq.UNDEFINED.and.rholaw(jj).eq.'CUSTOM') then
              write(*,*)'FATAL ERROR! Command RHOS(1,0) missing in input.dat'
              write(*,*)''
              write(52,*)'FATAL ERROR! Command RHOS(1,0) missing in input.dat'
              write(52,*)''
              error=.TRUE.
              else
              endif
          endif
      endif
      call check_drag_input(jj,ii,error_drag(jj))
      if(distr1) then                                                             ! Component identifier
        if(dotestchi(jj).and.siglevchi(jj).eq.undefined) then
        write(*,140)jj
        write(52,140)jj
        error=.TRUE.
        else
        endif
        call check_distribution(jj,error_grainsize(jj))
      else
        if(phi50(1).eq.UNDEFINED) then
                if(d50mm(1).eq.UNDEFINED) then
                write(*,*)'FATAL ERROR! Command PHI50(1) or "D50MM(1)" missing in input.dat'
                write(*,*)''
                write(52,*)'FATAL ERROR! Command PHI50(1) or "D50MM(1)" missing in input.dat'
                write(52,*)''
                error=.TRUE.
                else
                endif
        else
        endif
        if(sorting(1).eq.UNDEFINED) then
        write(*,*)'FATAL ERROR! Command SORTING(1) missing in input.dat'
        write(*,*)''
        write(52,*)'FATAL ERROR! Command SORTING(1) missing in input.dat'
        write(52,*)''
        error=.TRUE.
        else
        endif
        if(nclass(1).eq.UNDEFINED_I) then
        write(*,*)'FATAL ERROR! Command NCLASS(1) missing in input.dat'
        write(*,*)''
        write(52,*)'FATAL ERROR! Command NCLASS(1) missing in input.dat'
        write(52,*)''
        error=.TRUE.
        else
        endif
      endif
      select case (model)
             case ('TWOLAYERS')
                if(dens_ent.eq.UNDEFINED) then
                write(*,*)'FATAL ERROR! Command DENS_ENT missing in input.dat'
                write(*,*)''
                write(52,*)'FATAL ERROR! Command DENS_ENT missing in input.dat'
                write(52,*)''
                error=.TRUE.
                else
                endif
                if(dm_ent.eq.UNDEFINED) then
                write(*,*)'FATAL ERROR! Command DM_ENT missing in input.dat'
                write(*,*)''
                write(52,*)'FATAL ERROR! Command DM_ENT missing in input.dat'
                write(52,*)''
                error=.TRUE.
                else
                endif
             case ('TWOCOMPONENTS')
                jj=jj+1
                if(rholaw(jj).eq.undefined_c) then
                write(*,*)'FATAL ERROR! Command RHOLAW(2) missing in input.dat'
                write(*,*)''
                write(52,*)'FATAL ERROR! Command RHOLAW(2) missing in input.dat'
                write(52,*)''
                else
                    if(rholaw(jj).ne.'POLLENA'.and.rholaw(jj).ne.'AVERNO2'.and.rholaw(jj).ne.'AMS'&
     &              .and.rholaw(jj).ne.'POMPEI'.and.rholaw(jj).ne.'SIAL_XX'.and.rholaw(jj).ne.'FEM_XX'&
     &              .and.rholaw(jj).ne.'LITHIC'.and.rholaw(jj).ne.'MERCATO'.and.rholaw(jj).ne.'ASTRONI'&
     &              .and.rholaw(jj).ne.'CUSTOM') then
                    write(*,139)jj
                    write(52,139)jj
                    error=.TRUE.
                    else
                        if(rhos(jj,0).eq.UNDEFINED.and.rholaw(jj).eq.'CUSTOM') then
                        write(*,*)'FATAL ERROR! Command RHOS(2,0) missing in input.dat'
                        write(*,*)''
                        write(52,*)'FATAL ERROR! Command RHOS(2,0) missing in input.dat'
                        write(52,*)''
                        error=.TRUE.
                        else
                        endif
                    endif
                endif
                call check_drag_input(jj,ii,error_drag(jj))
                if(distr2) then
                        if(dotestchi(jj).and.siglevchi(jj).eq.undefined) then
                        write(*,140)jj
                        write(52,140)jj
                        error=.TRUE.
                        else
                        endif
                        call check_distribution(jj,error_grainsize(jj))
                else
                       if(phi50(2).eq.UNDEFINED) then
                                 if(d50mm(2).eq.UNDEFINED) then
                                 write(*,*)'FATAL ERROR! Command PHI50(2) or D50MM(2) missing in input.dat'
                                 write(*,*)''
                                 write(52,*)'FATAL ERROR! Command PHI50(2) or D50MM(2) missing in input.dat'
                                 write(52,*)''
                                 error=.TRUE.
                                 else
                                 endif
                       else
                       endif
                       if(sorting(2).eq.UNDEFINED) then
                       write(*,*)'FATAL ERROR! Command SORTING(2) missing in input.dat'
                       write(*,*)''
                       write(52,*)'FATAL ERROR! Command SORTING(2) missing in input.dat'
                       write(52,*)''
                       error=.TRUE.
                       else
                       endif
                       if(nclass(2).eq.UNDEFINED_I) then
                       write(*,*)'FATAL ERROR! Command NCLASS(2) missing in input.dat'
                       write(*,*)''
                       write(52,*)'FATAL ERROR! Command NCLASS(2) missing in input.dat'
                       write(52,*)''
                       error=.TRUE.
                       else
                       endif
                endif
      end select
      ii=1
      jjstart=1
      if(deprates) then
      
      if(dep_model.eq.undefined_c) then
      dep_model='DEFAULT'
      else
          if(dep_model.ne.'DELLINO_2018') then
          write(*,*)'FATAL ERROR! DEP_MODEL must be either DEFAULT or DELLINO_2018'
          write(*,*)''
          write(52,*)'FATAL ERROR! DEP_MODEL must be either DEFAULT or DELLINO_2018'
          write(52,*)''
          error=.TRUE.
          endif
      endif
      
      if(ncomp.eq.UNDEFINED_I) then
      write(*,*)'FATAL ERROR! Command NCOMP missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command NCOMP missing in input.dat'
      write(52,*)''
      error=.TRUE.
      else
           if(model.eq.'TWOLAYERS'.and.ncomp.lt.1) then
           write(*,*)'FATAL ERROR! NCOMP MUST be greater than or equal to 1'
           write(*,*)''
           write(52,*)'FATAL ERROR! NCOMP MUST be greater than or equal to 1'
           write(52,*)''
           error=.TRUE.
           else
             if(model.eq.'TWOCOMPONENTS'.and.ncomp.lt.2) then
             write(*,*)'FATAL ERROR! NCOMP MUST be greater than or equal to 2'
             write(*,*)''
             write(52,*)'FATAL ERROR! NCOMP MUST be greater than or equal to 2'
             write(52,*)''
             error=.TRUE.
             else
             endif
           endif
           if(model.eq.'TWOLAYERS'.and.ncomp.ne.1) jjstart=2
           if(model.eq.'TWOCOMPONENTS'.and.ncomp.ne.2) jjstart=3
           do jj=jjstart,ncomp
           call check_distribution(jj,error_grainsize(jj))
           enddo
           do jj=1,ncomp
           call check_drag_input(jj,ii,error_drag(jj))
           enddo
      call check_deprates(error_deprates)
      endif
      else
      endif

!      if(deprates.or.distr1.or.distr2) then

!      endif
      do i=1,ncomp_checks
      error_drag_final=error_drag_final.or.error_drag(i)
      error_distr=error_distr.or.error_grainsize(i)
      enddo
  567 if(error.or.error_drag_final.or.error_distr.or.error_deprates) call exit_pyflow


  139 format('FATAL ERROR! RHOLAW(',i1,') selection wrong. Please choose one of the following possibilities:',/,&
     &'POLLENA, AVERNO2, AMS, POMPEI, MERCATO, ASTRONI, SIAL_XX, FEM_XX, LITHIC, CUSTOM',//)
  140 format('FATAL ERROR! Command SIGLEVCHI(',i1,') missing in input.dat',//)
  152 format('FATAL ERROR! Command RHO_FLOW(',i1,') missing in input.dat',//)
  153 format('FATAL ERROR! Command ZTOT_FLOW(',i1,') missing in input.dat',//)
  154 format('FATAL ERROR! Command PN_FLOW(',i1,') missing in input.dat',//)
  155 format('FATAL ERROR! Command USH_FLOW(',i1,') missing in input.dat',//)
  156 format('FATAL ERROR! Command Z0_FLOW(',i1,') missing in input.dat',/,&
     &'Z0_FLOW(',i1,') mandatory if Z0_REC(',i1,') = .FALSE.',//)
  157 format('FATAL ERROR! Command PNS_FLOW(',i1,') missing in input.dat',//)
      end subroutine check_data

      subroutine check_distribution(jj,error_distr)
      use inoutdata; use nrtype
      implicit none
      integer :: i,jj,kk
      logical :: error_distr
      real(dp):: wt_input_check                         ! Sum of input weight fraction to check if it is 1
      error_distr=.FALSE.
      if(dphi(jj).eq.undefined) then
      write(*,143)jj
      write(52,143)jj
      error_distr=.TRUE.
      else
      endif
      if(phimin(jj).eq.undefined) then
      write(*,144)jj
      write(52,144)jj
      error_distr=.TRUE.
      else
      endif
      if(phimax(jj).eq.undefined) then
      write(*,145)jj
      write(52,145)jj
      error_distr=.TRUE.
      else
      endif
      xnmax=(phimax(jj)-phimin(jj))/dphi(jj)
      nclass(jj)=idint(xnmax)+1
      do i=1,nclass(jj)
         if(weight(jj,i).eq.undefined) then              ! Check only weights here, density and shape parameterd are needed for deposition rate calculations
         write(*,141)jj,i
         write(52,141)jj,i
         error_distr=.TRUE.
         else
         endif
      enddo
      
  141 format('FATAL ERROR! Command WEIGHT(',i1,','i2,') missing in input.dat',//)
  143 format('FATAL ERROR! Command DPHI(',i1,') missing in input.dat',//)
  144 format('FATAL ERROR! Command PHIMIN(',i1,') missing in input.dat',//)
  145 format('FATAL ERROR! Command PHIMAX(',i1,') missing in input.dat',//)
      end subroutine check_distribution

      subroutine check_deprates(error_deprates)
      use inoutdata; use nrtype
      implicit none
      integer :: i,jj
      logical :: error_deprates
      real(dp) :: wt_input_check
      error_deprates=.FALSE.
      
          if(zlam_massive.eq.undefined) then
          write(*,*)'WARNING! Command ZLAM_MASSIVE missing in input.dat'
          write(*,*)'Deposition rate calculations will not take this into account'
          write(*,*)''
          write(52,*)'WARNING! Command ZLAM_MASSIVE missing in input.dat'
          write(52,*)'Deposition rate calculations will not take this into account'
          write(52,*)''
          else
          endif
          
          do jj=1,ncomp
              if(merge_classes(jj).and.sensmerge(jj).eq.undefined) then
              write(*,158)jj
              write(52,158)jj
              error_deprates=.TRUE.
              else
              endif
              if(rholaw(jj).eq.undefined_c) then
              write(*,146)jj
              write(52,146)jj
              error_deprates=.TRUE.
              else
                  if(rholaw(jj).ne.'POLLENA'.and.rholaw(jj).ne.'AVERNO2'.and.rholaw(jj).ne.'AMS'&
     &            .and.rholaw(jj).ne.'POMPEI'.and.rholaw(jj).ne.'SIAL_XX'.and.rholaw(jj).ne.'FEM_XX'&
     &            .and.rholaw(jj).ne.'LITHIC'.and.rholaw(jj).ne.'MERCATO'.and.rholaw(jj).ne.'CUSTOM'&
     &            .and.rholaw(jj).ne.'ASTRONI') then
                  write(*,148)jj
                  write(52,148)jj
                  error_deprates=.TRUE.
                  else
                      if(rholaw(jj).eq.'CUSTOM') then
                                 if(rho_custom(jj).eq.undefined_c) then
                                 write(*,149)jj
                                 write(52,149)jj
                                 error_deprates=.TRUE.
                                 else
                                     select case (rho_custom(jj))
                                     case('CONSTANT')
                                     if(rhos(jj,0).eq.undefined) then
                                     write(*,150)jj
                                     write(52,150)jj
                                     error_deprates=.TRUE.
                                     else
                                     endif
                                     case('VARIABLE')
                                     do i=1,nclass(jj)
                                        if(rhos(jj,i).eq.undefined) then
                                        write(*,147)jj,i
                                        write(52,147)jj,i
                                        error_deprates=.TRUE.
                                        else
                                        endif
                                     enddo
                                     case default
                                     write(*,151)jj
                                     write(52,151)jj
                                     end select
                                 endif
                      endif
                  endif
              endif
          enddo
          
      if(input_weight.eq.undefined_c) then
      write(*,*)'FATAL ERROR! Command INPUT_WEIGHT missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command INPUT_WEIGHT missing in input.dat'
      write(52,*)''
      error_deprates=.TRUE.
      else
             if(input_weight.eq.'WT') then
             wt_input_check=0.d0
                do jj=1,ncomp
                   do i=1,nclass(jj)
                   wt_input_check=wt_input_check+weight(jj,i)
                   enddo
                enddo
                write(52,*)wt_input_check
                write(*,*)'wt_tot',wt_input_check
                if(wt_input_check.lt.100.d0) then
                write(*,174)wt_input_check
                write(52,174)wt_input_check
                else
                      if(wt_input_check-100.d0.gt.0.1d0) then
                        write(*,175)wt_input_check
                        write(52,175)wt_input_check
                        error_deprates=.TRUE.
                        else
                        write(*,176)wt_input_check
                        write(52,176)wt_input_check
                      endif
                endif
             else
!             wt_input_check=0.d0
!                do jj=1,ncomp
!                   do i=1,nclass(jj)
!                   wt_input_check=wt_input_check+weight(jj,i)
!                   enddo
!                enddo
!                if(wt_input_check.lt.WTOT_SAMPLE) then
!                write(*,177)wt_input_check,WTOT_SAMPLE
!                write(52,177)wt_input_check,WTOT_SAMPLE
!                else
!                      if(wt_input_check.gt.WTOT_SAMPLE) then
!                      write(*,178)wt_input_check,WTOT_SAMPLE
!                      write(52,178)wt_input_check,WTOT_SAMPLE
!                      error_deprates=.TRUE.
!                      else
!                      endif
!                endif
             endif
      endif
      if(dep_median.eq.undefined) then
      write(*,*)'FATAL ERROR! Command DEP_MEDIAN missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command DEP_MEDIAN missing in input.dat'
      write(52,*)''
      error_deprates=.TRUE.
      endif
      if(rhos_median.eq.undefined) then
      write(*,*)'FATAL ERROR! Command RHOS_MEDIAN missing in input.dat'
      write(*,*)''
      write(52,*)'FATAL ERROR! Command RHOS_MEDIAN missing in input.dat'
      write(52,*)''
      error_deprates=.TRUE.
      endif
  146 format('FATAL ERROR! Command RHOLAW(',i1,') missing in input.dat',//)
  147 format('FATAL ERROR! Command RHOS(',i1,','i2,') missing in input.dat',//)
  148 format('FATAL ERROR! RHOLAW selection wrong. Please choose one of the following possibilities:',/,&
     &'POLLENA, AVERNO2, AMS, POMPEI, MERCATO, ASTRONI, SIAL_XX, FEM_XX, LITHIC, CUSTOM',//)
  149 format('FATAL ERROR! Command RHO_CUSTOM(',i1,') missing in input.dat',//)
  150 format('FATAL ERROR! Command RHOS(',i1,',0) missing in input.dat',//)
  151 format('FATAL ERROR! Wrong assignment for command RHO_CUSTOM(',i1,')',/,&
     &'Please choose one of the following possibilities',/,&
     &'WARNING! IT IS CASE SENSITIVE!',/,'CONSTANT  or  VARIABLE',//)
  158 format('FATAL ERROR! Command SENSMERGE(',i1,') missing in input.dat',//)
  174 format('WARNING! Sum of weight fractions of the sample is ',f7.3,' < 100%',/,&
     &'Weight fraction will be normalized to sum up to 100%',//)
  175 format('FATAL ERROR! Sum of weight fractions of the sample ',f7.3,' is significantly > 100%',//)
  176 format('WARNING! Sum of weight fractions of the sample ',f7.3,' is slightly > 100%',//)
  177 format('WARNING! Sum of masses of grainsize classes is ',f7.3,' < ',f7.3,//)
  178 format('FATAL ERROR! Sum of masses of grainsize classes is ',f7.3,' > ',f7.3,//)
      end subroutine check_deprates

      subroutine check_drag_input(jj,ii,error_d)
      use inoutdata; use nrtype
      implicit none
      integer :: i,ii,jj
      logical :: error_d
      error_d=.FALSE.
      if(cdlaw(jj).eq.UNDEFINED_C) then
      write(*,142)jj
      write(52,142)jj
      error_d=.TRUE.
      else
              select case (cdlaw(jj))
              case ('SPHERE')
                   write(*,173)jj
                   write(52,173)jj
              case ('HAIDLEV')
                   if(ii.eq.0) then
                        if(sphericity(jj,0).eq.UNDEFINED) then
                        write(*,120)jj
                        write(52,120)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,0)=sphericity(jj,0)
                        endif
                   else
                     do i=1,nclass(jj)
                        if(sphericity(jj,i).eq.UNDEFINED) then
                                if(sphericity(jj,0).eq.UNDEFINED) then
                                write(*,161)jj,i,jj
                                write(52,161)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                sphericity(jj,i)=sphericity(jj,0)
                                shpar1(jj,i)=sphericity(jj,i)
                                endif
                        else
                        shpar1(jj,i)=sphericity(jj,i)
                        endif
                     enddo
                   endif
              case ('SWAMOJ')
                   if(ii.eq.0) then
                        if(corey(jj,0).eq.UNDEFINED) then
                        write(*,121)jj
                        write(52,121)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,0)=corey(jj,0)
                        endif
                   else
                     do i=1,nclass(jj)
                        if(corey(jj,i).eq.UNDEFINED) then
                                if(corey(jj,0).eq.UNDEFINED) then
                                write(*,162)jj,i,jj
                                write(52,162)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                corey(jj,i)=corey(jj,0)
                                shpar1(jj,i)=corey(jj,i)
                                endif
                        else
                        shpar1(jj,i)=corey(jj,i)
                        endif
                     enddo
                   endif
              case ('GANSER')
                   if(ii.eq.0) then
                        !Sphericity
!                        if(isometric(jj,0) then
                        if(isometric(jj)) then
                               if(sphericity(jj,0).eq.UNDEFINED) then
                               write(*,120)jj
                               write(52,120)jj
                               error_d=.TRUE.
                               else
                               shpar1(jj,0)=sphericity(jj,0)
                               endif
                        !Volume equivalent sphere diameter
                        else
                                if(voleqsphd(jj,0).eq.UNDEFINED) then
                                write(*,122)jj
                                write(52,122)jj
                                error_d=.TRUE.
                                else
                                shpar2(jj,0)=voleqsphd(jj,0)
                                endif
                                !Equal projected area circle diameter
                                if(circeqard(jj,0).eq.UNDEFINED) then
                                write(*,123)jj
                                write(52,123)jj
                                error_d=.TRUE.
                                else
                                shpar3(jj,0)=circeqard(jj,0)
                                endif
                        endif
                        !Isometric particle
!                        if(.not.isometric(jj)) shpar4(jj)=isometric(jj)
                        shpar4(jj)=isometric(jj)
                   else
                        if(isometric(jj)) then
                        !Sphericity
                                do i=1,nclass(jj)
                                   if(sphericity(jj,i).eq.UNDEFINED) then
                                        if(sphericity(jj,0).eq.UNDEFINED) then
                                        write(*,161)jj,i,jj
                                        write(52,161)jj,i,jj
                                        error_d=.TRUE.
                                        exit
                                        else
                                        sphericity(jj,i)=sphericity(jj,0)
                                        shpar1(jj,i)=sphericity(jj,i)
                                        endif
                                   else
                                   shpar1(jj,i)=sphericity(jj,i)
                                   endif
                                enddo
                         else
                        !Volume equivalent sphere diameter
                                do i=1,nclass(jj)
                                   if(voleqsphd(jj,i).eq.UNDEFINED) then
                                        if(voleqsphd(jj,0).eq.UNDEFINED) then
                                        write(*,163)jj,i,jj
                                        write(52,163)jj,i,jj
                                        error_d=.TRUE.
                                        exit
                                        else
                                        voleqsphd(jj,i)=voleqsphd(jj,0)
                                        shpar2(jj,i)=voleqsphd(jj,i)
                                        endif
                                   else
                                   shpar2(jj,i)=voleqsphd(jj,i)
                                   endif
                                enddo
                       !Equal projected area circle diameter
                                do i=1,nclass(jj)
                                   if(circeqard(jj,i).eq.UNDEFINED) then
                                        if(circeqard(jj,0).eq.UNDEFINED) then
                                        write(*,164)jj,i,jj
                                        write(52,164)jj,i,jj
                                        error_d=.TRUE.
                                        exit
                                        else
                                        circeqard(jj,i)=circeqard(jj,0)
                                        shpar3(jj,i)=circeqard(jj,i)
                                        endif
                                   else
                                   shpar3(jj,i)=circeqard(jj,i)
                                   endif
                                enddo
                      !Isometric particle
!                     do i=1,nclass(jj)
!                        if(.not.isometric(jj,i)) shpar4(jj,i)=isometric(jj,i)
!                     enddo
                        endif
!                        if(.not.isometric(jj)) shpar4(jj)=isometric(jj)
                        shpar4(jj)=isometric(jj)
                   endif
!                   write(52,*)'CHECK_DRAG',isometric(jj)
              case ('CHIEN')
                   if(ii.eq.0) then
                        if(sphericity(jj,0).eq.UNDEFINED) then
                        write(*,120)jj
                        write(52,120)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,0)=sphericity(jj,0)
                        endif
                   else
                     do i=1,nclass(jj)
                        if(sphericity(jj,i).eq.UNDEFINED) then
                                if(sphericity(jj,0).eq.UNDEFINED) then
                                write(*,161)jj,i,jj
                                write(52,161)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                sphericity(jj,i)=sphericity(jj,0)
                                shpar1(jj,i)=sphericity(jj,i)
                                endif
                        else
                        shpar1(jj,i)=sphericity(jj,i)
                        endif
                     enddo
                   endif
              case ('TRANCONG')
                   if(ii.eq.0) then
                        !Circularity
                        if(circularity(jj,0).eq.UNDEFINED) then
                        write(*,124)jj
                        write(52,124)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,i)=circularity(jj,0)
                        endif
                        !Flatness ratio
                        if(flatratio(jj,0).eq.UNDEFINED) then
                        write(*,125)jj
                        write(52,125)jj
                        error_d=.TRUE.
                        else
                        shpar2(jj,0)=flatratio(jj,0)
                        endif
                   else
                        !Circularity
                     do i=1,nclass(jj)
                        if(circularity(jj,i).eq.UNDEFINED) then
                                if(circularity(jj,0).eq.UNDEFINED) then
                                write(*,165)jj,i,jj
                                write(52,165)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                circularity(jj,i)=circularity(jj,0)
                                shpar1(jj,i)=circularity(jj,i)
                                endif
                        else
                        shpar1(jj,i)=circularity(jj,i)
                        endif
                     enddo
                        !Flatness ratio
                     do i=1,nclass(jj)
                        if(flatratio(jj,i).eq.UNDEFINED) then
                                if(flatratio(jj,0).eq.UNDEFINED) then
                                write(*,166)jj,i,jj
                                write(52,166)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                flatratio(jj,i)=flatratio(jj,0)
                                shpar2(jj,i)=flatratio(jj,i)
                                endif
                        else
                        shpar2(jj,i)=flatratio(jj,i)
                        endif
                     enddo
                   endif
              case ('DELLINO')
                        !Shape factor
                   if(ii.eq.0) then
                        if(shapefact(jj,0).eq.UNDEFINED) then
                        write(*,126)jj
                        write(52,126)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,0)=shapefact(jj,0)
                        endif
                   else
                     do i=1,nclass(jj)
                        if(shapefact(jj,i).eq.UNDEFINED) then
                                if(shapefact(jj,0).eq.UNDEFINED) then
                                write(*,168)jj,i,jj
                                write(52,168)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                shapefact(jj,i)=shapefact(jj,0)
                                shpar1(jj,i)=shapefact(jj,i)
                                endif
                        else
                        shpar1(jj,i)=shapefact(jj,i)
                        endif
                     enddo
                   endif
              case ('HOLZSOMM')
                   if(ii.eq.0) then
                        !Sphericity
                        if(sphericity(jj,0).eq.UNDEFINED) then
                        write(*,120)jj
                        write(52,120)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,i)=sphericity(jj,0)
                        endif
                        !Lengthwise sphericity
                        if(longspher(jj,0).eq.UNDEFINED) then
                        write(*,127)jj
                        write(52,127)jj
                        error_d=.TRUE.
                        else
                        shpar2(jj,0)=longspher(jj,0)
                        endif
                        !Crosswise sphericity
                        if(crossspher(jj,0).eq.UNDEFINED) then
                        write(*,128)jj
                        write(52,128)jj
                        error_d=.TRUE.
                        else
                        shpar3(jj,0)=crossspher(jj,0)
                        endif
                   else
                        !Sphericity
                     do i=1,nclass(jj)
                        if(sphericity(jj,i).eq.UNDEFINED) then
                                if(sphericity(jj,0).eq.UNDEFINED) then
                                write(*,161)jj,i,jj
                                write(52,161)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                sphericity(jj,i)=sphericity(jj,0)
                                shpar1(jj,i)=sphericity(jj,i)
                                endif
                        else
                        shpar1(jj,i)=sphericity(jj,i)
                        endif
                     enddo
                        !Lengthwise sphericity
                     do i=1,nclass(jj)
                        if(longspher(jj,i).eq.UNDEFINED) then
                                if(longspher(jj,0).eq.UNDEFINED) then
                                write(*,169)jj,i,jj
                                write(52,169)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                longspher(jj,i)=longspher(jj,0)
                                shpar2(jj,i)=longspher(jj,i)
                                endif
                        else
                        shpar2(jj,i)=longspher(jj,i)
                        endif
                     enddo
                       !Crosswise sphericity
                     do i=1,nclass(jj)
                        if(crossspher(jj,i).eq.UNDEFINED) then
                                if(crossspher(jj,0).eq.UNDEFINED) then
                                write(*,170)jj,i,jj
                                write(52,170)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                crossspher(jj,i)=crossspher(jj,0)
                                shpar3(jj,i)=crossspher(jj,i)
                                endif
                        else
                        shpar3(jj,i)=crossspher(jj,i)
                        endif
                     enddo
                   endif
              case ('DIOGMELE')
                        !Shape factor
                   if(ii.eq.0) then
                        if(shapefact(jj,0).eq.UNDEFINED) then
                        write(*,126)jj
                        write(52,126)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,0)=shapefact(jj,0)
                        endif
                   else
                     do i=1,nclass(jj)
                        if(shapefact(jj,i).eq.UNDEFINED) then
                                if(shapefact(jj,0).eq.UNDEFINED) then
                                write(*,168)jj,i,jj
                                write(52,168)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                shapefact(jj,i)=shapefact(jj,0)
                                shpar1(jj,i)=shapefact(jj,i)
                                endif
                        else
                        shpar1(jj,i)=shapefact(jj,i)
                        endif
                     enddo
                   endif
              case ('DIOG2017')
                if(fractal(jj)) then
                   ! 3D Fractal Dimension
                   if(ii.eq.0) then
                        if(fractdim(jj,0).eq.UNDEFINED) then
                        write(*,124)jj
                        write(52,124)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,i)=fractdim(jj,0)
                        endif
                   else
                     do i=1,nclass(jj)
                        if(fractdim(jj,i).eq.UNDEFINED) then
                                if(fractdim(jj,0).eq.UNDEFINED) then
                                write(*,171)jj,i,jj
                                write(52,171)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                fractdim(jj,i)=fractdim(jj,0)
                                shpar1(jj,i)=fractdim(jj,i)
                                endif
                        else
                        shpar1(jj,i)=fractdim(jj,i)
                        endif
                     enddo
                   endif
                else
                   ! 3D Sphericity
                   if(ii.eq.0) then
                        if(sphericity(jj,0).eq.UNDEFINED) then
                        write(*,120)jj
                        write(52,120)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,0)=sphericity(jj,0)
                        endif
                   else
                     do i=1,nclass(jj)
                        if(sphericity(jj,i).eq.UNDEFINED) then
                                if(sphericity(jj,0).eq.UNDEFINED) then
                                write(*,161)jj,i,jj
                                write(52,161)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                sphericity(jj,i)=sphericity(jj,0)
                                shpar1(jj,i)=sphericity(jj,i)
                                endif
                        else
                        shpar1(jj,i)=sphericity(jj,i)
                        endif
                     enddo
                   endif
                endif
              case ('DIOG2018')
                        !Shape factor
                   if(ii.eq.0) then
                        if(shapefact(jj,0).eq.UNDEFINED) then
                        write(*,126)jj
                        write(52,126)jj
                        error_d=.TRUE.
                        else
                        shpar1(jj,0)=shapefact(jj,0)
                        endif
                   else
                     do i=1,nclass(jj)
                        if(shapefact(jj,i).eq.UNDEFINED) then
                                if(shapefact(jj,0).eq.UNDEFINED) then
                                write(*,168)jj,i,jj
                                write(52,168)jj,i,jj
                                error_d=.TRUE.
                                exit
                                else
                                shapefact(jj,i)=shapefact(jj,0)
                                shpar1(jj,i)=shapefact(jj,i)
                                endif
                        else
                        shpar1(jj,i)=shapefact(jj,i)
                        endif
                     enddo
                   endif
              case default
                   write(52,172)jj
                   write(*,172)jj
                   error_d=.TRUE.
              end select
      endif
  142 format('FATAL ERROR! Command CDLAW(',i2,') missing in input.dat'//)
  120 format('FATAL ERROR! Command SPHERICITY(',i1,',0) missing in input.dat',//)
  121 format('FATAL ERROR! Command COREY(',i1,',0) missing in input.dat',//)
  122 format('FATAL ERROR! Command VOLEQSPHD(',i1,',0) missing in input.dat',//)
  123 format('FATAL ERROR! Command CIRCEQARD(',i1,',0) missing in input.dat',//)
  124 format('FATAL ERROR! Command CIRCULARITY(',i1,',0) missing in input.dat',//)
  125 format('FATAL ERROR! Command FLATRATIO(',i1,',0) missing in input.dat',//)
  126 format('FATAL ERROR! Command SHAPEFACT(',i1,',0) missing in input.dat',//)
  127 format('FATAL ERROR! Command LONGSPHER(',i1,',0) missing in input.dat',//)
  128 format('FATAL ERROR! Command CROSSSPHER(',i1,',0) missing in input.dat',//)
  129 format('FATAL ERROR! Command FRACTDIM(',i1,',0) missing in input.dat',//)
  161 format('FATAL ERROR! Both commands SPHERICITY(',i1,','i2,') and SPHERICITY(',i1,',0) missing in input.dat',//)
  162 format('FATAL ERROR! Both commands COREY(',i1,','i2,') and COREY(',i1,',0) missing in input.dat',//)
  163 format('FATAL ERROR! Both commands VOLEQSPHD(',i1,','i2,') and VOLEQSPHD(',i1,',0) missing in input.dat',//)
  164 format('FATAL ERROR! Both commands CIRCEQARD(',i1,','i2,') and CIRCEQARD(',i1,',0) missing in input.dat',//)
  165 format('FATAL ERROR! Both commands CIRCULARITY(',i1,','i2,') and CIRCULARITY(',i1,',0) missing in input.dat',//)
  166 format('FATAL ERROR! Both commands FLATRATIO(',i1,','i2,') and FLATRATIO(',i1,',0) missing in input.dat',//)
  167 format('WARNING! SHAPEFACTLAW(',i1,') missing in input.dat',/,'Reading input particle shape factor values',//)
  168 format('FATAL ERROR! Both commands SHAPEFACT(',i1,','i2,') and SHAPEFACT(',i1,',0) missing in input.dat',//)
  169 format('FATAL ERROR! Both commands LONGSPHER(',i1,','i2,') and LONGSPHER(',i1,',0) missing in input.dat',//)
  170 format('FATAL ERROR! Both commands CROSSSPHER(',i1,','i2,') and CROSSSPHER(',i1,',0) missing in input.dat',//)
  171 format('FATAL ERROR! Both commands FRACTDIM(',i1,','i2,') and FRACTDIM(',i1,',0) missing in input.dat',//)
  172 format('FATAL ERROR! Cd law selection wrong for component 'i1'.',/,&
     &'Please choose one of the following possibilities',/,&
     &'WARNING! IT IS CASE SENSITIVE!',/,'HAIDLEV (Haider and Levenspiel 1989)',/,&
     &'SWAMOJ (Swamee and Ojha 1991)',/,'GANSER (Ganser 1993)',/,'CHIEN (Chien 1994)',/,&
     &'TRANCONG (Tran-Cong et al. 2004)',/,'DELLINO (Dellino et al. 2005)',/,&
     &'HOLZSOMM (Holzer and Sommerfeld 2008)',/,'DIOGMELE (Dioguardi and Mele 2015)',/,  &
     &'DIOG2017 (Dioguardi et al. 2017)','DIOG2018 (Dioguardi et al. 2018)',//)
  173 format('WARNING! The drag coefficient of the sphere will be calculated for component ',i1,//)
      end subroutine check_drag_input
  
      subroutine exit_pyflow
      implicit none
      write(*,*)''
      write(*,*)'PYFLOW 2.0 is going to be stopped'
      write(52,*)''
      write(52,*)'PYFLOW 2.0 is going to be stopped'
      stop
      end subroutine exit_pyflow
      
