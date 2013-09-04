! MODULE: diffr
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for computing diffracted power in periodic problems.
MODULE diffr
  USE source
  USE nfields
  USE common

  IMPLICIT NONE

  INTEGER, PARAMETER :: max_prdsrc = 4
  INTEGER, PARAMETER :: fresnel_nquad = 25

CONTAINS
  FUNCTION diffracted_power(b, wlindex, dindex, r0, xorder, yorder) RESULT(power)
    TYPE(batch), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: wlindex, dindex, xorder, yorder
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0

    REAL (KIND=dp) :: power, wl, omega
    REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: qwx, ptx, qwy, pty
    COMPLEX (KIND=dp), DIMENSION(3) :: e, h, einc, hinc
    COMPLEX (KIND=dp), DIMENSION(2) :: eff
    COMPLEX (KIND=dp) :: ri, ri1
    REAL (KIND=dp), DIMENSION(3) :: pt
    REAL (KIND=dp), DIMENSION(2) :: kt
    REAL (KIND=dp) :: irrinc, hdx, hdy, A
    INTEGER :: n, m, nx, ny
    REAL (KIND=dp), DIMENSION(3) :: xaxis, yaxis
    TYPE(prdnfo), POINTER :: prd

    prd => b%prd(b%domains(dindex)%gf_index)

    wl = b%sols(wlindex)%wl
    omega = 2.0_dp*pi*c0/wl
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri
    ri1 = b%media(b%domains(1)%medium_index)%prop(wlindex)%ri

    ! Select the number of integration points based on wavelength and period.
    !nx = NINT(prd%dx/b%sols(wlindex)%wl*20)
    !ny = NINT(prd%dy/b%sols(wlindex)%wl*20)
    nx = 51
    ny = 51

    ! Make sure that the numbers are odd.
    IF(MOD(nx,2)==0) THEN
       nx = nx + 1
    END IF

    IF(MOD(ny,2)==0) THEN
       ny = ny + 1
    END IF

    ALLOCATE(qwx(1:nx), ptx(1:nx), qwy(1:ny), pty(1:ny))

    hdx = prd%dx*0.5_dp
    hdy = prd%dy*0.5_dp

    xaxis = (/prd%cp, prd%sp, 0.0_dp/)
    yaxis = (/0.0_dp, 1.0_dp, 0.0_dp/)

    ! Compute weights and nodes from Simpson's rule.
    CALL get_simpsons_weights(-hdx, hdx, nx-1, qwx)
    CALL get_simpsons_points(-hdx, hdx, nx-1, ptx)
    CALL get_simpsons_weights(-hdy, hdy, ny-1, qwy)
    CALL get_simpsons_points(-hdy, hdy, ny-1, pty)

    eff(:) = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(nx,ny,prd,xorder,yorder,r0,xaxis,ptx,yaxis,pty,dindex,b,wlindex,qwx,qwy,eff,omega,ri)&
    !$OMP PRIVATE(m,n,kt,pt,e,h,einc,hinc)
    !$OMP DO REDUCTION(+:eff) SCHEDULE(STATIC)
    DO m=1,ny
       DO n=1,nx

          ! Lattice vector.
          kt = (/prd%coef(prd%cwl)%k0x + 2.0_dp*PI*(xorder/(prd%dx*prd%cp) - yorder*prd%sp/(prd%dy*prd%cp)) ,&
               prd%coef(prd%cwl)%k0y + 2.0_dp*PI*yorder/prd%dy/)

          pt = r0 + xaxis*ptx(n) + yaxis*pty(m)

          CALL scat_fields(b%domains(dindex)%mesh, b%ga, b%sols(wlindex)%x, b%mesh%nedges,&
               omega, ri, prd, pt, e, h)

          IF(dindex==1) THEN
             CALL src_fields(b%src, omega, ri, pt, einc, hinc)

             e = e + einc
             h = h + hinc
          END IF

          eff = eff + qwx(n)*qwy(m)*(/e(1), e(2)/)*EXP(-(0,1)*dotr(kt,pt(1:2)))
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    ! Area of the unit cell.
    A = prd%dx*prd%dy*prd%cp

    ! cp is the Jacobian of the area integration.
    eff = eff*prd%cp/A

    ! The relative power diffracted to 0th order in the given domain.
    power = (ABS(eff(1))**2 + ABS(eff(2))**2)*ri/ri1

    DEALLOCATE(qwx, ptx, qwy, pty)

  END FUNCTION diffracted_power
END MODULE diffr
