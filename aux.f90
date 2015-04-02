! MODULE: aux
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Auxiliary routines for various purposes, such as string manipulation, file access and simple
! mathematical tools. This module is a good candidate for refactorization in the future.
MODULE aux
  USE constants

  IMPLICIT NONE

CONTAINS
  ! Finds nval smallest-in-magnitude values in val, whose size is dim.
  ! Returns the linear indices to these values.
  FUNCTION find_smallest(val, dim, nval) RESULT(ind)
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: val
    INTEGER, INTENT(IN) :: dim, nval

    INTEGER, DIMENSION(nval) :: ind
    COMPLEX (KIND=dp), DIMENSION(dim) :: tmpval
    INTEGER :: n
    INTEGER, DIMENSION(1) :: loc
    REAL (KIND=dp) :: mval

    tmpval(:) = val(:)

    mval = MAXVAL(ABS(val))

    DO n=1,nval
       loc = MINLOC(ABS(tmpval))
       tmpval(loc(1)) = mval

       ind(n) = loc(1)
    END DO
  END FUNCTION find_smallest

  ! Returns the three-letter filename extension of the filename.
  FUNCTION getext(filename) RESULT(ext)
    CHARACTER (LEN=*), INTENT(IN) :: filename
    CHARACTER (LEN=3) :: ext
    INTEGER :: n

    n = LEN_TRIM(filename)

    ext = filename((n-2):n)
  END FUNCTION getext

  ! Replaces the three-letter filename extension in filename with ext.
  SUBROUTINE replace_ext(filename, ext)
    CHARACTER (LEN=*), INTENT(INOUT) :: filename
    CHARACTER (LEN=3), INTENT(IN) :: ext
    INTEGER :: n

    n = LEN_TRIM(filename)
    filename((n-2):n) = ext
  END SUBROUTINE replace_ext

  ! Linear interpolation between start and end with t in [0,1].
  FUNCTION linterp(start, end, t) RESULT(res)
    REAL (KIND=dp), INTENT(IN) :: start, end, t
    REAL (KIND=dp) :: res

    res = start + t*(end - start)
  END FUNCTION linterp

  ! Reads a matrix from file.
  SUBROUTINE get_matrix(filename, mat, nrows, ncols)
    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER :: fid = 10, iovar, i
    INTEGER, PARAMETER :: nattempts = 10
    REAL (KIND=dp), DIMENSION(nrows,ncols), INTENT(INOUT) :: mat
    INTEGER, INTENT(IN) :: nrows, ncols

    ! Attempt to open the file for nattempts times before resorting
    ! to failure. This is so that multiple instances of the program
    ! can access the same refractive index files simultaneously.

    DO i=1,nattempts
       OPEN(fid, FILE=TRIM(filename), ACTION='READ', IOSTAT=iovar)
       IF(iovar>0) THEN
          CALL SLEEP(1)
       END IF
    END DO

    IF(i==nattempts .AND. iovar>0) THEN
       WRITE(*,*) 'Could not open matrix data file!'
       STOP
    END IF

    iovar = 0
    DO i=1,nrows
       READ(fid,*) mat(i,1:ncols)
    END DO

    CLOSE(fid)

  END SUBROUTINE get_matrix

  ! Returns a complex refractive index value from the file for the given wavelength wl.
  ! Linear interpolation is used.
  FUNCTION get_refind(filename, wl) RESULT(res)
    CHARACTER (LEN=*), INTENT(IN) :: filename
    REAL (KIND=dp), INTENT(IN) :: wl
    COMPLEX (KIND=dp) :: res
    REAL (KIND=dp) :: wl1, wl2, nr1, ni1, nr2, ni2, t
    INTEGER :: fid = 10, iovar, i
    INTEGER, PARAMETER :: nattempts = 10

    ! Attempt to open the file for nattempts times before resorting
    ! to failure. This is so that multiple instances of the program
    ! can access the same refractive index files simultaneously.

    DO i=1,nattempts
       OPEN(fid, FILE=TRIM(filename), ACTION='READ', IOSTAT=iovar)
       IF(iovar>0) THEN
          CALL SLEEP(1)
       END IF
    END DO

    IF(i==nattempts .AND. iovar>0) THEN
       WRITE(*,*) 'Could not open refractive index lookup file!'
       STOP
    END IF

    iovar = 0
    READ(fid,*) wl1, nr1, ni1

    IF(wl<wl1) THEN
       WRITE(*,*) 'Wavelength below refractive index sample range!'
       WRITE(*,*) 'Using unity refractive index.'
       res = 1.0_dp
       CLOSE(fid)
       RETURN
    END IF

    DO WHILE(iovar==0)
       READ(fid,*) wl2, nr2, ni2

       IF(wl1<=wl .AND. wl<=wl2 .AND. wl1/=wl2) THEN
          t = (wl-wl1)/(wl2 - wl1)
          res = CMPLX(linterp(nr1, nr2, t), linterp(ni1, ni2, t))
          CLOSE(fid)
          RETURN
       ELSE IF(wl==wl1) THEN
          res = CMPLX(nr1, ni1)
          CLOSE(fid)
          RETURN
       ELSE IF(wl==wl2) THEN
          res = CMPLX(nr2, ni2)
          CLOSE(fid)
          RETURN
       ELSE IF(wl1>wl2) THEN
          WRITE(*,*) 'Refractive index data disordered!'
          CLOSE(fid)
          STOP
       END IF

       wl1 = wl2
       nr1 = nr2
       ni1 = ni2
    END DO

    CLOSE(fid)

    WRITE(*,*) 'Wavelength above refractive index sample range!'
    WRITE(*,*) 'Using unity refractive index.'
    res = 1.0_dp
  END FUNCTION get_refind

  ! Finds from list the value that is closest to val and returns index to it.
  FUNCTION find_closest(val, list) RESULT(ind)
    REAL (KIND=dp), INTENT(IN) :: val
    REAL (KIND=dp), DIMENSION(:), INTENT(IN) :: list
    INTEGER :: ind, n
    REAL (KIND=dp) :: dist

    ind = 1
    dist = ABS(val-list(1))

    IF(SIZE(list)==1) THEN
       RETURN
    END IF

    DO n=2,SIZE(list)
       IF(ABS(val-list(n))<dist) THEN
          ind = n
          dist = ABS(val-list(n))
       END IF
    END DO
  END FUNCTION find_closest

  ! Writes matrix data to file.
  SUBROUTINE write_data(filename, data)
    REAL (KIND=dp), DIMENSION(:,:), INTENT(IN) :: data
    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER :: fid = 10, iovar, i, j

    OPEN(fid, FILE=TRIM(filename), ACTION='WRITE', IOSTAT=iovar)
    IF(iovar>0) THEN
       WRITE(*,*) 'Could not open output file' // filename // '!'
       STOP
    END IF

    DO i=1,SIZE(data,1)
       DO j=1,SIZE(data,2)
          WRITE(fid, '(EN15.3)', ADVANCE='NO') data(i,j)
       END DO
       WRITE(fid, '(/)')
    END DO

    CLOSE(fid)
  END SUBROUTINE write_data

  ! Writes complex matrix data to file. Complex numbers are written in the form x+iy.
  SUBROUTINE write_cdata(filename, data)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: data
    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER :: fid = 10, iovar, i, j

    OPEN(fid, FILE=TRIM(filename), ACTION='WRITE', IOSTAT=iovar)
    IF(iovar>0) THEN
       WRITE(*,*) 'Could not open output file' // filename // '!'
       STOP
    END IF

    DO i=1,SIZE(data,1)
       DO j=1,SIZE(data,2)
          WRITE(fid, '(EN15.3,"+",EN15.3,"i")', ADVANCE='NO') REAL(data(i,j)), IMAG(data(i,j))
       END DO
       WRITE(fid, '(/)')
    END DO

    CLOSE(fid)
  END SUBROUTINE write_cdata

  ! Evaluates the Drude model of noble metal dispersion at angular frequency omega.
  FUNCTION drude(omega, epsinf, omegap, gamma) RESULT(eps)
    COMPLEX (KIND=dp), INTENT(IN) :: omega
    REAL (KIND=dp), INTENT(IN) :: epsinf, omegap, gamma
    COMPLEX (KIND=dp) :: eps

    eps = epsinf - omegap**2/(omega*(omega + (0,1)*gamma))
  END FUNCTION drude

  ! Returns the direction vector for spherical polar angles theta and phi.
  FUNCTION get_dir(theta, phi) RESULT(dir)
    REAL (KIND=dp), INTENT(IN) :: theta, phi
    REAL (KIND=dp), DIMENSION(3) :: dir
    
    dir = (/SIN(theta)*COS(phi), SIN(theta)*SIN(phi), COS(theta)/)
  END FUNCTION get_dir

  ! Returns a 'polarization' vector of plane-wave propagating in direction given
  ! by get_dir. The angle psi is the rotation of the polarization along the direction of propagation.
  FUNCTION get_pol(theta, phi, psi) RESULT(pol)
    REAL (KIND=dp), INTENT(IN) :: theta, phi, psi
    REAL (KIND=dp), DIMENSION(3) :: pol

    pol = (/COS(psi)*SIN(phi) - SIN(psi)*COS(theta)*COS(phi),&
         -COS(psi)*COS(phi) - SIN(psi)*COS(theta)*SIN(phi),&
         SIN(psi)*SIN(theta)/)
  END FUNCTION get_pol
END MODULE aux
