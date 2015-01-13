PROGRAM main
  USE interface
  USE OMP_LIB
  !USE strat
  !USE pec
  IMPLICIT NONE

  INTEGER :: max_threads

  !CHARACTER (LEN=256) :: meshname, scalestr, omegastr
  !REAL (KIND=dp) :: scale, omega

  !CALL getarg(1, meshname)
  !CALL getarg(2, scalestr)
  !CALL getarg(3, omegastr)

  !READ(scalestr, '(EN15.3)') scale
  !READ(omegastr, '(EN15.3)') omega

  !WRITE(*,*) 'mesh:  ', TRIM(meshname)
  !WRITE(*,*) 'scale: ', scale
  !WRITE(*,*) 'omega: ', omega

  !CALL pec_modes(TRIM(ADJUSTL(meshname)), scale, omega)
  !STOP

  !CALL plot_integGr()
  !STOP

  WRITE(*,*) 'Method of Moments solver by Jouni Makitalo'
  WRITE(*,*) 'Tampere University of Technology, Optics Laboratory'
  WRITE(*,*) '---'

  max_threads = OMP_GET_MAX_THREADS()
  CALL OMP_SET_NUM_THREADS(max_threads)
  WRITE(*,'(A,I0,:)') ' Number of threads for OpenMP: ', max_threads

  CALL msgloop()

END PROGRAM main
