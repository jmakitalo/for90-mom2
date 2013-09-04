PROGRAM main
  USE interface
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER :: max_threads

  !CALL compare_greenprd()
  !CALL compare_cheby()
  !STOP

  WRITE(*,*) 'Method of Moments solver by Jouni Makitalo'
  WRITE(*,*) 'Tampere University of Technology, Optics Laboratory'
  WRITE(*,*) '---'

  max_threads = OMP_GET_MAX_THREADS()
  CALL OMP_SET_NUM_THREADS(max_threads)
  WRITE(*,'(A,I0,:)') ' Number of threads for OpenMP: ', max_threads

  CALL msgloop()

END PROGRAM main
