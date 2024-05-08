      subroutine read_data
         use inoutdata; use nrtype
         implicit none

! The maximum allowed length of the command
         INTEGER, PARAMETER :: comm_max_length = 80
! Save the entire line of the input.dat file
         CHARACTER(LEN=512) :: comm_string
! Row number in the input file
         INTEGER :: row_num
! Length of command without comments
         INTEGER :: comm_uncomm
         CHARACTER(len=256) :: command
! Integer function which returns COMMENT_INDEX
         INTEGER, EXTERNAL :: search_comment
! Blank line function
         LOGICAL, EXTERNAL :: blank_row
         INTEGER :: IOS

         include 'commands.inc'

         inquire (file='input.dat', exist=lEXISTS)
         IF (.NOT. lEXISTS) THEN
            write (*, *) 'FATAL ERROR! input.dat file is missing'
            write (flog, *) 'FATAL ERROR! input.dat file is missing'
            call exit_pyflow
         ELSE
         END IF

         OPEN (40, FILE='input.dat')
         row_num = 0

         READ_LP: DO
            READ (40, "(A)", IOSTAT=IOS) comm_string
            IF (IOS < 0) EXIT READ_LP
            IF (IOS > 0) THEN
            END IF

            row_num = row_num + 1
            comm_uncomm = search_comment(comm_string, LEN(comm_string)) - 1
            CALL delete_comm(comm_string, comm_uncomm + 1, LEN(comm_string))
            IF (comm_uncomm <= 0) CYCLE READ_LP           ! comment line
            IF (blank_row(comm_string)) CYCLE READ_LP ! blank line

! Standard model input parameters.
            command = ''; command = '&INPUT_DATA '// &
 trim(adjustl(comm_string(1:comm_uncomm)))//'/'
            READ (command, NML=INPUT_DATA, IOSTAT=IOS)
            if (IOS > 0) then
               write (flog, 1254) row_num, command(12:comm_uncomm)
               write (*, 1254) row_num, command(12:comm_uncomm)
1254           format('****INPUT ERROR****', /, 'Unknown command at line ', i3, /, a512, /, &
                      'PYFLOW_2.5 is going to be stopped')
               stop
            else
            end if

         END DO READ_LP

      end

      INTEGER FUNCTION search_comment(row, comm_max_length)
         IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   input data row
         CHARACTER(len=*) row
!
!                   maximum column of input data row to search
         INTEGER comm_max_length
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!                   the number of designated comment characters
         INTEGER, PARAMETER :: comm_types = 2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                   loop indicies
         INTEGER :: L, L2
!
!                   the comment characters
         CHARACTER, DIMENSION(comm_types) :: comm_symbols
!-----------------------------------------------
!
!     The function search_comment returns the index to where a comment
!     character was found in the input data row.  Equals comm_max_length + 1
!     if no-comment characters in the row
!
!
         DATA comm_symbols/'#', '!'/
!
         DO L = 1, comm_max_length
            DO L2 = 1, comm_types
               IF (row(L:L) == comm_symbols(L2)) THEN
                  search_comment = L
                  RETURN
               END IF
            END DO
         END DO
         search_comment = comm_max_length + 1
!
         RETURN
      END FUNCTION search_comment

      SUBROUTINE delete_comm(row, comm_start, comm_max_length)

         IMPLICIT NONE
! Input data row
         CHARACTER(len=*), intent(INOUT) :: row
!Start of comments
         INTEGER, intent(IN) :: comm_start
! Maximum column of input data row to search
         INTEGER, intent(IN) :: comm_max_length

! Local Variables:
!---------------------------------------------------------------------//
! Loop index
         INTEGER :: L

         DO L = comm_start, comm_max_length
            row(L:L) = ' '
         END DO

         RETURN
      END SUBROUTINE delete_comm

      LOGICAL FUNCTION blank_row(row)

         IMPLICIT NONE

         CHARACTER :: row*(*)

         INTEGER :: L

         blank_row = .FALSE.
         DO L = 1, len(row)
            IF (row(L:L) /= ' ' .and. row(L:L) /= '    ') RETURN
         END DO

         blank_row = .TRUE.
         RETURN
      END FUNCTION blank_row
