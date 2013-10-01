! MODULE: greenprd
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines to import data of pre-computed periodic Green's function
! and its gradient. Implemets also functions for the evaluation of
! these functions based on the imported data.
MODULE greenprd
  USE linalg

  IMPLICIT NONE

  TYPE prdcoef
     COMPLEX (KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: samples
     COMPLEX (KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: samplesz
     REAL (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: rho
     REAL (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: kt
     REAL (KIND=dp) :: E, k0x, k0y, wl, range
     INTEGER :: ic, jc, kc, n, m, np, npz
     COMPLEX (KIND=dp) :: ri, ri0, k
  END TYPE prdcoef

  TYPE prdnfo
     CHARACTER (LEN=256) :: filename
     INTEGER :: type
     INTEGER :: nwl, cwl
     REAL (KIND=dp) :: dx, dy, dz, phi, phasex1, phasex2, phasey
     REAL (KIND=dp) :: cp, sp, pwtheta, pwphi
     LOGICAL :: oblique ! oblique incidence
     TYPE(prdcoef), DIMENSION(:), ALLOCATABLE :: coef
  END TYPE prdnfo

CONTAINS
  ! si: start index for interpolation i.e. use t to interpolate
  !     values corresponding to indices si and si+1.
  ! t: interpolation value in range [0,1]
  ! x1: first value of the range.
  ! x2: last value of the range.
  ! x: the value, for which si and t are requested.
  ! n: number of samples in the range.
  SUBROUTINE interpIndex(x1, x2, x, n, nx, si, t)
    REAL (KIND=dp), INTENT(IN) :: x1, x2
    INTEGER, INTENT(IN) :: n, nx
    REAL (KIND=dp), DIMENSION(nx), INTENT(IN) :: x
    INTEGER, DIMENSION(nx), INTENT(OUT) :: si
    REAL (KIND=dp), DIMENSION(nx), INTENT(OUT) :: t

    REAL (KIND=dp) :: l

    ! Length of one interval.
    l = (x2-x1)/(n - 1)

    si(:) = FLOOR((x(:)-x1)/l) + 1

    ! Check for index over- and underflows.
    !WHERE(si>=n)
    !   si = n-1
    !END WHERE
    !WHERE(si<1)
    !   si = 1
    !END WHERE

    t(:) = (x(:) - (x1 + l*(si(:)-1)))/l
  END SUBROUTINE interpIndex

  SUBROUTINE vGp(ro, rp, prd, nrp, near, g)
    INTEGER, INTENT(IN) :: nrp
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: ro
    REAL (KIND=dp), DIMENSION(3,nrp), INTENT(IN) :: rp
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near

    COMPLEX (KIND=dp), DIMENSION(nrp), INTENT(INOUT) :: g

    CALL vGprd_interp_1d(ro, rp, prd, nrp, near, g)
  END SUBROUTINE vGp

  SUBROUTINE vGprd_interp_1d(ro, rp, prd, nrp, near, g)
    INTEGER, INTENT(IN) :: nrp
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: ro
    REAL (KIND=dp), DIMENSION(3,nrp), INTENT(IN) :: rp
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near

    COMPLEX (KIND=dp), DIMENSION(nrp), INTENT(INOUT) :: g
    COMPLEX (KIND=dp), DIMENSION(nrp) :: phase, phase_inc,&
         phase1, phase2, phase1n, phase2m
    COMPLEX (KIND=dp) :: phase_inc2
    REAL (KIND=dp), DIMENSION(3,nrp) :: r
    REAL (KIND=dp), DIMENSION(2) :: rho, kt
    REAL (KIND=dp), DIMENSION(nrp) :: t, tz, Rnm
    REAL (KIND=dp) :: range
    INTEGER :: n, m, npoints, npointsz, nn, mm
    INTEGER, DIMENSION(nrp) :: ti, tzi

    g(:) = 0.0_dp

    r(1,:) = ro(1) - rp(1,:)
    r(2,:) = ro(2) - rp(2,:)
    r(3,:) = ro(3) - rp(3,:)

    npoints = prd%coef(prd%cwl)%np
    npointsz = prd%coef(prd%cwl)%npz
    range = prd%coef(prd%cwl)%range

    CALL interpIndex(-prd%dz, prd%dz, r(3,:), npointsz, nrp, tzi, tz)

    IF(prd%oblique) THEN
       phase_inc = EXP((0,1)*(prd%coef(prd%cwl)%k0x*r(1,:) + prd%coef(prd%cwl)%k0y*r(2,:)))
    ELSE
       phase_inc(:) = 1.0_dp
    END IF

    phase1 = EXP((0,1)*prd%phasex1*r(1,:))
    phase2 = EXP((0,1)*(prd%phasex2*r(1,:) + prd%phasey*r(2,:)))

    phase1n = phase1**(-prd%coef(prd%cwl)%n)

    DO nn=-prd%coef(prd%cwl)%n,prd%coef(prd%cwl)%n

       phase2m = phase2**(-prd%coef(prd%cwl)%m)

       DO mm=-prd%coef(prd%cwl)%m,prd%coef(prd%cwl)%m

          n = prd%coef(prd%cwl)%n + nn + 1
          m = prd%coef(prd%cwl)%m + mm + 1
          
          rho = prd%coef(prd%cwl)%rho(n,m,:)
          kt = prd%coef(prd%cwl)%kt(n,m,:)
        
          Rnm = SQRT((r(1,:)-rho(1))**2 + (r(2,:)-rho(2))**2 + r(3,:)**2)

          phase = phase_inc*phase1n*phase2m
        
          CALL interpIndex(0.0_dp, range, Rnm, npoints, nrp, ti, t)
        
          ! Interpolate spatial series.
          g = g + prd%coef(prd%cwl)%samples(n,m,ti,1)*(1-t) &
               + prd%coef(prd%cwl)%samples(n,m,ti+1,1)*t

          ! Interpolate spectral series.
          g = g + (prd%coef(prd%cwl)%samplesz(n,m,tzi,1)*(1-tz) &
               + prd%coef(prd%cwl)%samplesz(n,m,tzi+1,1)*tz)*phase

          ! For far elements, add singular parts.
          IF(near==.FALSE. .AND. ABS(nn)<2 .AND. ABS(mm)<2) THEN
             IF(prd%oblique) THEN
                phase_inc2 = EXP((0,1)*(prd%coef(prd%cwl)%k0x*rho(1) &
                     + prd%coef(prd%cwl)%k0y*rho(2)))
             ELSE
                phase_inc2 = 1.0_dp
             END IF
             g = g - (-1/(4*pi*Rnm) + (prd%coef(prd%cwl)%k**2)*Rnm/(8*pi))*phase_inc2
          END IF

          phase2m = phase2m*phase2
       END DO

       phase1n = phase1n*phase1
    END DO
  END SUBROUTINE vGprd_interp_1d

  SUBROUTINE vgradGp(ro, rp, prd, nrp, near, gg)
    INTEGER, INTENT(IN) :: nrp
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: ro
    REAL (KIND=dp), DIMENSION(3,nrp), INTENT(IN) :: rp
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near

    COMPLEX (KIND=dp), DIMENSION(3,nrp), INTENT(INOUT) :: gg

    CALL vgradGprd_interp_1d(ro, rp, prd, nrp, near, gg)
    !CALL vgradGprd_interp_1d_fd(ro, rp, prd, nrp, near, gg)
  END SUBROUTINE vgradGp

  SUBROUTINE vgradGprd_interp_1d_fd(ro, rp, prd, nrp, near, gg)
    INTEGER, INTENT(IN) :: nrp
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: ro
    REAL (KIND=dp), DIMENSION(3,nrp), INTENT(IN) :: rp
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near

    COMPLEX (KIND=dp), DIMENSION(3,nrp), INTENT(INOUT) :: gg
    COMPLEX (KIND=dp), DIMENSION(nrp) :: g1, g2
    REAL (KIND=dp), PARAMETER :: d = 1D-9

    CALL vGprd_interp_1d(ro-(/d,0.0_dp,0.0_dp/), rp, prd, nrp, near, g1)
    CALL vGprd_interp_1d(ro+(/d,0.0_dp,0.0_dp/), rp, prd, nrp, near, g2)
    gg(1,:) = g2 - g1

    CALL vGprd_interp_1d(ro-(/0.0_dp,d,0.0_dp/), rp, prd, nrp, near, g1)
    CALL vGprd_interp_1d(ro+(/0.0_dp,d,0.0_dp/), rp, prd, nrp, near, g2)
    gg(2,:) = g2 - g1

    CALL vGprd_interp_1d(ro-(/0.0_dp,0.0_dp,d/), rp, prd, nrp, near, g1)
    CALL vGprd_interp_1d(ro+(/0.0_dp,0.0_dp,d/), rp, prd, nrp, near, g2)
    gg(3,:) = g2 - g1

    gg = -gg/(2*d)

  END SUBROUTINE vgradGprd_interp_1d_fd

  SUBROUTINE vgradGprd_interp_1d(ro, rp, prd, nrp, near, gg)
    INTEGER, INTENT(IN) :: nrp
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: ro
    REAL (KIND=dp), DIMENSION(3,nrp), INTENT(IN) :: rp
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near

    COMPLEX (KIND=dp), DIMENSION(3,nrp), INTENT(INOUT) :: gg
    COMPLEX (KIND=dp), DIMENSION(nrp) :: phase, cspat_int, cspect_int,&
         cspecz_int, phase_inc, phase1, phase2, phase1n, phase2m
    COMPLEX (KIND=dp) :: phase_inc2
    REAL (KIND=dp), DIMENSION(3,nrp) :: r
    REAL (KIND=dp), DIMENSION(2) :: rho, kt
    REAL (KIND=dp) :: range
    REAL (KIND=dp), DIMENSION(nrp) :: t, tz, Rnm, xnm, ynm, znm
    INTEGER :: n, m, npoints, npointsz, nn, mm
    INTEGER, DIMENSION(nrp) :: ti, tzi

    r(1,:) = ro(1) - rp(1,:)
    r(2,:) = ro(2) - rp(2,:)
    r(3,:) = ro(3) - rp(3,:)

    gg(:,:) = 0.0_dp

    npoints = prd%coef(prd%cwl)%np
    npointsz = prd%coef(prd%cwl)%npz
    range = prd%coef(prd%cwl)%range
    
    CALL interpIndex(-prd%dz, prd%dz, r(3,:), npointsz, nrp, tzi, tz)

    IF(prd%oblique) THEN
       phase_inc = EXP((0,1)*(prd%coef(prd%cwl)%k0x*r(1,:) + prd%coef(prd%cwl)%k0y*r(2,:)))
    ELSE
       phase_inc(:) = 1.0_dp
    END IF

    phase1 = EXP((0,1)*prd%phasex1*r(1,:))
    phase2 = EXP((0,1)*(prd%phasex2*r(1,:) + prd%phasey*r(2,:)))

    phase1n = phase1**(-prd%coef(prd%cwl)%n)

    DO nn=-prd%coef(prd%cwl)%n,prd%coef(prd%cwl)%n

       phase2m = phase2**(-prd%coef(prd%cwl)%m)

       DO mm=-prd%coef(prd%cwl)%m,prd%coef(prd%cwl)%m

          n = prd%coef(prd%cwl)%n + nn + 1
          m = prd%coef(prd%cwl)%m + mm + 1
          
          rho = prd%coef(prd%cwl)%rho(n,m,:)
          kt = prd%coef(prd%cwl)%kt(n,m,:)
          
          Rnm = SQRT((r(1,:)-rho(1))**2 + (r(2,:)-rho(2))**2 + r(3,:)**2)

          phase = phase_inc*phase1n*phase2m
          
          CALL interpIndex(0.0_dp, range, Rnm, npoints, nrp, ti, t)

          WHERE(Rnm/=0)          
             xnm = (r(1,:)-rho(1))/Rnm
             ynm = (r(2,:)-rho(2))/Rnm
             znm = r(3,:)/Rnm
          END WHERE

          ! Interpolate spatial series.          
          cspat_int = prd%coef(prd%cwl)%samples(n,m,ti,2)*(1-t) +&
               prd%coef(prd%cwl)%samples(n,m,ti+1,2)*t

          WHERE(Rnm==0)
             cspat_int = 0.0_dp
          END WHERE

          ! For far elements, add singular parts.
          IF(near==.FALSE. .AND. ABS(nn)<2 .AND. ABS(mm)<2) THEN
             IF(prd%oblique) THEN
                phase_inc2 = EXP((0,1)*(prd%coef(prd%cwl)%k0x*rho(1) &
                     + prd%coef(prd%cwl)%k0y*rho(2)))
             ELSE
                phase_inc2 = 1.0_dp
             END IF

             cspat_int = cspat_int + (1/(4*pi*(Rnm**3))&
                  + (prd%coef(prd%cwl)%k**2)/(8*pi*Rnm))*phase_inc2*Rnm
          END IF

          ! Interpolate spectral series.
          cspect_int = prd%coef(prd%cwl)%samplesz(n,m,tzi,2)*(1-tz) + prd%coef(prd%cwl)%samplesz(n,m,tzi+1,2)*tz
          cspecz_int = prd%coef(prd%cwl)%samplesz(n,m,tzi,3)*(1-tz) + prd%coef(prd%cwl)%samplesz(n,m,tzi+1,3)*tz
          
          gg(1,:) = gg(1,:) + cspat_int*xnm + cspect_int*phase*kt(1)
          gg(2,:) = gg(2,:) + cspat_int*ynm + cspect_int*phase*kt(2)
          gg(3,:) = gg(3,:) + cspat_int*znm + cspecz_int*phase

          phase2m = phase2m*phase2
       END DO

       phase1n = phase1n*phase1
    END DO

    !WHERE(r(1,:)==0 .AND. r(2,:)==0 .AND. r(3,:)==0)
    !   gg(1,:) = 0.0_dp
    !   gg(2,:) = 0.0_dp
    !   gg(3,:) = 0.0_dp
    !END WHERE
  END SUBROUTINE vgradGprd_interp_1d

  SUBROUTINE import_pgfw_interp_1d(filename, prd)
    CHARACTER (LEN=*), INTENT(IN) :: filename
    TYPE(prdnfo), INTENT(INOUT) :: prd
    INTEGER :: fid=10, n, wlind
    CHARACTER (LEN=32) :: tag
    COMPLEX (KIND=dp) :: k
    REAL (KIND=dp) :: k0
    INTEGER :: a,b,c,d,aa,bb

    WRITE(*,*) 'Importing periodic Green function data'

    OPEN(fid, FILE=filename, ACTION='READ')

    ! Read header.

    READ(fid,*) tag, tag
    IF(tag/='2dinterp1d') THEN
       WRITE(*,*) 'Invalid periodicity type ', tag, '!'
       CLOSE(fid)
       STOP
    ELSE
       WRITE(*,*) '2D periodicity'
    END IF

    prd%type = prd_2d

    READ(fid,*) tag, prd%dx, prd%dy, prd%dz, prd%phi
    READ(fid,*) tag, prd%pwtheta, prd%pwphi
    READ(fid,*) tag, prd%nwl

    WRITE(*,'(A,F6.1,A)') 'dx: ', prd%dx*1e9, ' nm'
    WRITE(*,'(A,F6.1,A)') 'dy: ', prd%dy*1e9, ' nm'
    WRITE(*,'(A,F6.1,A)') 'dz: ', prd%dz*1e9, ' nm'
    WRITE(*,'(A,F6.1,A)') 'phi: ', prd%phi*180/pi, ' deg'
    WRITE(*,'(A,F6.1,A)') 'pwtheta: ', prd%pwtheta*180/pi, ' deg'
    WRITE(*,'(A,F6.1,A)') 'pwphi: ', prd%pwphi*180/pi, ' deg'
    WRITE(*,'(I0,A)') prd%nwl, ' wavelengths'

    prd%cp = COS(prd%phi)
    prd%sp = SIN(prd%phi)

    IF(ABS(prd%pwtheta)>0.001_dp) THEN
       prd%oblique = .TRUE.
       WRITE(*,*) 'Oblique plane-wave incidence'
    ELSE
       prd%oblique = .FALSE.
       WRITE(*,*) 'Normal plane-wave incidence'
    END IF

    ! Precompute phase arguments.
    prd%phasex1 = 2*pi/(prd%dx*prd%cp)
    prd%phasex2 = -2*pi*prd%sp/(prd%dy*prd%cp)
    prd%phasey = 2*pi/prd%dy

    ALLOCATE(prd%coef(prd%nwl))

    ! Read PGF for each wavelength.
    DO n=1,prd%nwl
       READ(fid,*) tag, wlind
       READ(fid,*) tag, prd%coef(n)%wl
       READ(fid,*) tag, prd%coef(n)%ri
       READ(fid,*) tag, prd%coef(n)%ri0
       READ(fid,*) tag, prd%coef(n)%E
       READ(fid,*) tag, prd%coef(n)%np
       READ(fid,*) tag, prd%coef(n)%npz
       READ(fid,*) tag, prd%coef(n)%n, prd%coef(n)%m
       READ(fid,*) tag, prd%coef(n)%range

       ! Precompute wavenumbers and vectors.
       k = prd%coef(n)%ri*2.0_dp*pi/prd%coef(n)%wl
       prd%coef(n)%k = k

       k0 = prd%coef(n)%ri0*2.0_dp*pi/prd%coef(n)%wl
       prd%coef(n)%k0x = SIN(prd%pwtheta)*COS(prd%pwphi)*k0
       prd%coef(n)%k0y = SIN(prd%pwtheta)*SIN(prd%pwphi)*k0

       ! Allocate data for interpolants.
       ALLOCATE(prd%coef(n)%samples(prd%coef(n)%n*2+1, prd%coef(n)%m*2+1, prd%coef(n)%np, 2))
       ALLOCATE(prd%coef(n)%samplesz(prd%coef(n)%n*2+1, prd%coef(n)%m*2+1, prd%coef(n)%npz, 3))

       ! Allocate data for precomputed variables.
       ALLOCATE(prd%coef(n)%rho(1+2*prd%coef(n)%n, 1+2*prd%coef(n)%m, 2))
       ALLOCATE(prd%coef(n)%kt(1+2*prd%coef(n)%n, 1+2*prd%coef(n)%m, 2))

       ! Precompute rho and kt.
       DO a=-prd%coef(n)%n,prd%coef(n)%n
          DO b=-prd%coef(n)%m,prd%coef(n)%m
             aa = prd%coef(n)%n + a + 1
             bb = prd%coef(n)%m + b + 1

             prd%coef(n)%rho(aa,bb,:) = (/a*prd%dx*prd%cp, b*prd%dy + a*prd%dx*prd%sp/)
             prd%coef(n)%kt(aa,bb,:) = (/prd%coef(n)%k0x + 2*pi*(a/(prd%dx*prd%cp)&
                  - b*prd%sp/(prd%dy*prd%cp)),&
                  prd%coef(n)%k0y + 2*pi*b/prd%dy/)
          END DO
       END DO

       ! Read spatial part samples of PGF.
       DO a=1,2
          DO b=1,prd%coef(n)%np
             DO c=1,(2*prd%coef(n)%m+1)
                DO d=1,(2*prd%coef(n)%n+1)
                   READ(fid,*) prd%coef(n)%samples(d,c,b,a)
                END DO
             END DO
          END DO
       END DO

       ! Read spectral part samples of PGF.
       DO a=1,3
          DO b=1,prd%coef(n)%npz
             DO c=1,(2*prd%coef(n)%m+1)
                DO d=1,(2*prd%coef(n)%n+1)
                   READ(fid,*) prd%coef(n)%samplesz(d,c,b,a)
                END DO
             END DO
          END DO
       END DO

    END DO

    CLOSE(fid)
  END SUBROUTINE import_pgfw_interp_1d

  SUBROUTINE clear_prd(prd)
    TYPE(prdnfo), INTENT(INOUT) :: prd
    INTEGER :: n
    
    IF(ALLOCATED(prd%coef)) THEN
       DO n=1,prd%nwl
          IF(ALLOCATED(prd%coef(n)%samples)) THEN
             DEALLOCATE(prd%coef(n)%samples)
          END IF

          IF(ALLOCATED(prd%coef(n)%samplesz)) THEN
             DEALLOCATE(prd%coef(n)%samplesz)
          END IF

          IF(ALLOCATED(prd%coef(n)%rho)) THEN
             DEALLOCATE(prd%coef(n)%rho)
          END IF

          IF(ALLOCATED(prd%coef(n)%kt)) THEN
             DEALLOCATE(prd%coef(n)%kt)
          END IF
       END DO
       DEALLOCATE(prd%coef)
    END IF
  END SUBROUTINE clear_prd
END MODULE greenprd
