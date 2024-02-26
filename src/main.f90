      USE inoutdata; USE nrtype
      implicit none
      character(len=255) :: cwd, slope_str, outdir_windows, outdir_linux
      integer :: i
      logical :: converged
      
      ! INITIALIZE VARIABLES
      call init_variables
      ! READ DATA FROM input.dat FILE
      call read_data
      ! CHECK DATA FOR THEIR COMPLETENESS AND CORRECTNESS
      call check_data
      ! INITIALIZE slopes
      call init_slopes
      ! CREATE OUTPUT FOLDERS AND FILES
      call outputs
      ! DEFINE GRAINSIZE PARAMETERS
      call grainsize
      if(.not.only_deprates) then
            ! CALL TWOLAYER OR TWOCOMPONENT DEPENDING ON THE USER'S CHOICE #######################
            select case (model)
                  case ('TWOLAYERS')
                  !   Two layer model
                        call twolayer
                  !   Two components model
                  case ('TWOCOMPONENTS')
                        call twocomponent
            end select                     
            do i = 1, size(slopes)
                  slope_ground = slopes(i)
                  output_dir = output_dirs(i)
                  write(flog, *)'**************************************'
                  write(flog, 123)slopes(i)
123               format('Profiles, deposition rates and PDFs calculation at slope ', f6.2,/)
                  ! CALL PROFILES FOR CALCULATING VERTICAL PROFILES OF FLUID-DYNAMIC VARIABLES ######################
                  call profiles(converged)    
                  if(.not.converged) then
                        if(i.eq.size(slopes)) then 
                              call exit_pyflow
                        else
                              cycle
                        endif
                  endif
                  ! CALL DEPOSITION FOR CALCULATING DEPOSITION TIME AND RATES ######################
                  call deposition                            
                  ! WRITE A SUMMARY OF INPUT DATA
                  call write_data_summary
                  ! Write results
                  call write_results
                  ! PROBABILITY FUNCTIONS CALCULATION FOR DYNAMIC PRESSURE AND PARTICLE CONCENTRATION##################################
                  call probfunction      
            enddo
      else
            output_dir = output_dirs(1)
            ! CALL DEPOSITION FOR CALCULATING DEPOSITION TIME AND RATES ######################
            call deposition 
            ! WRITE A SUMMARY OF INPUT DATA
            call write_data_summary
            ! Write results
            call write_results
            ! PROBABILITY FUNCTIONS CALCULATION FOR DYNAMIC PRESSURE AND PARTICLE CONCENTRATION##################################
            call probfunction
      endif      
      end

      subroutine init_slopes
      USE inoutdata; USE nrtype
      implicit none
      integer :: i
      if (slope_ground .eq. undefined) then
            n_slopes = (slope_ground_max - slope_ground_min) / delta_slope + 1
            allocate(slopes(n_slopes))
            slopes(1) = slope_ground_min
            do i = 2, n_slopes
                  slopes(i) = slopes(i-1) + delta_slope
            enddo
      else
            allocate(slopes(1))
            slopes(1) = slope_ground
            n_slopes = 1
      endif
      end subroutine init_slopes
      
      subroutine outputs
            USE inoutdata; USE nrtype
            implicit none
            character(len=255) :: cwd, slope_str, slope_str_temp, outdir, output_file
            integer :: i, i_slope, n_occurrence_linuxsep, n_occurrence_windowssep
            call getcwd(cwd)
            n_occurrence_windowssep = scan(cwd, '\')
            n_occurrence_linuxsep = scan(cwd, '/')
            if(n_occurrence_linuxsep.gt.0.and.n_occurrence_windowssep.eq.0) then
                  path_sep = '/'
            else
                  path_sep = '\'
            endif
            grainsize_dir = trim(cwd)//trim(path_sep)//'grainsize_analysis'
            call system('mkdir '//grainsize_dir)
            if(only_deprates) then
                  allocate(output_dirs(1))
                  output_dirs(1) = trim(cwd)
            else      
                  allocate(output_dirs(n_slopes))         
                  do i_slope = 1, size(slopes)
                        slope_str = 'slope_'   
                        write(slope_str_temp, 40)slopes(i_slope)
40                      format(f6.2)
                        do i=1,len(slope_str_temp)   
                              if(slope_str_temp(i:i).ne.' ')slope_str=trim(slope_str)//trim(slope_str_temp(i:i))   
                        end do
                        outdir = trim(cwd)//trim(path_sep)//'results'
                        call system('mkdir '//outdir)
                        if(slopes(i_slope).ne.undefined) then
                              outdir = trim(cwd)//trim(path_sep)//'results'//trim(path_sep)//slope_str
                        endif
                        call system('mkdir '//outdir)
                        output_dirs(i_slope) = outdir
                        output_file = trim(outdir)//trim(path_sep)//'results.dat'
                        open(fout,file=output_file)
                        write(fout,*) '###PROGRAM PYFLOW 2.4 by Fabio Dioguardi###'
                        close(fout)
                  enddo
            endif
      end subroutine outputs
      
