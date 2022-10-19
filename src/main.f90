      USE inoutdata; USE nrtype
      implicit none

      open (50, file='results.dat')
      open (52, file='log.dat')
      write (*, *) '###PROGRAM PYFLOW 2.3 by Fabio Dioguardi###'
      write (50, *) '###PROGRAM PYFLOW 2.3 by Fabio Dioguardi###'
      write (52, *) '###PROGRAM PYFLOW 2.3 by Fabio Dioguardi###'
      write (52, *) 'LOG FILE'
      write (52, *) ''

! INITIALIZE VARIABLES
      call init_variables
! READ DATA FROM input.dat FILE
      call read_data
! CHECK DATA FOR THEIR COMPLETENESS AND CORRECTNESS
      call check_data
! DEFINE GRAINSIZE PARAMETERS
      call grainsize
      if (only_deprates) goto 123

! CALL TWOLAYER OR TWOCOMPONENT DEPENDING ON THE USER'S CHOICE #######################
      select case (model)
      case ('TWOLAYERS')
!     Two layer model
         call twolayer
!     Two components model
      case ('TWOCOMPONENTS')
         call twocomponent
      end select
! CALL PROFILES FOR CALCULATING VERTICAL PROFILES OF FLUID-DYNAMIC VARIABLES ######################
      call profiles
! CALL DEPOSITION FOR CALCULATING DEPOSITION TIME AND RATES ######################
123   if (deprates) call deposition

! WRITE A SUMMARY OF INPUT DATA
      call write_data_summary

! Write results
      call write_results

! PROBABILITY FUNCTIONS CALCULATION FOR DYNAMIC PRESSURE AND PARTICLE CONCENTRATION##################################
      call probfunction
   end

