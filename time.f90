! MODULE: time
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for measuring time and time intervals.
MODULE time
  IMPLICIT NONE

  INTEGER, PARAMETER :: r8=SELECTED_REAL_KIND(12,100)
  REAL (KIND=r8), PRIVATE, SAVE :: time_start

CONTAINS
  FUNCTION get_time() RESULT(sec)
    IMPLICIT NONE
    INTEGER, DIMENSION(8) :: t
    REAL (KIND=r8) :: sec

    CALL DATE_AND_TIME(values=t)
    sec = t(3)*86400.0_r8+t(5)*3600.0_r8+t(6)*60.0_r8+t(7)+t(8)/1000.0_r8
  END FUNCTION get_time

  SUBROUTINE timer_start()
    time_start = get_time()
  END SUBROUTINE timer_start

  FUNCTION timer_end() RESULT(sec)
    IMPLICIT NONE
    REAL (KIND=r8) :: sec
    sec = get_time() - time_start
  END FUNCTION timer_end

  FUNCTION sec_to_str(sec) RESULT(str)
    REAL (KIND=r8), INTENT(IN) :: sec
    CHARACTER (LEN=128) :: str
    INTEGER :: hours, minutes, seconds

    hours = FLOOR(sec/3600.0_r8)
    minutes = FLOOR((sec-hours*3600.0_r8)/60.0_r8)
    seconds = FLOOR((sec-hours*3600.0_r8-minutes*60.0_r8))
    write(str, '(I0,A,I0,A,I0,A)') hours, ' hours ', minutes, ' minutes ',&
         seconds, ' seconds '
  END FUNCTION sec_to_str
END MODULE time
