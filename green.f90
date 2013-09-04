! MODULE: green
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implements the fundamental Green's function and its gradient
! for the Helmholtz operator. Also the smoothed (singularity subtracted)
! versions are implemented here.
MODULE green
  USE linalg

  IMPLICIT NONE

    REAL (KIND=dp), DIMENSION(7), PARAMETER :: expcoef = (/1.0000000000000000_dp,&
         5.0000000000000000D-001,&
         1.6666666666666666D-001,&
         4.1666666666666664D-002,&
         8.3333333333333332D-003,&
         1.3888888888888889D-003,&
         1.9841269841269841D-004/)

    COMPLEX (KIND=dp), PARAMETER :: oofpi = 1.0_dp/(4.0_dp*pi),&
         iofpi = (0,1)/(4.0_dp*pi)

CONTAINS
  PURE FUNCTION Gs(r0, rp, k) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0, rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp) :: R
    COMPLEX (KIND=dp) :: res

    R = normr(r0-rp)

    IF(R==0) THEN
       res = k*iofpi
    ELSE
       res = oofpi*((EXP((0,1)*k*R)-1.0_dp)/R + 0.5_dp*(k**2)*R)
    END IF
  END FUNCTION Gs

  SUBROUTINE vGs(r0, rp, k, N, res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0
    REAL (KIND=dp), DIMENSION(3,N), INTENT(IN) :: rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    INTEGER, INTENT(IN) :: N
    COMPLEX (KIND=dp), DIMENSION(N), INTENT(INOUT) :: res

    ! Store distances temporarily to result.
    res = SQRT((r0(1)-rp(1,:))**2 + (r0(2)-rp(2,:))**2 + (r0(3)-rp(3,:))**2)

    WHERE(res==0)
       res = k*iofpi
    ELSEWHERE
       res = oofpi*((EXP((0,1)*k*res)-1.0_dp)/res + 0.5_dp*(k**2)*res)
    END WHERE
  END SUBROUTINE vGs

  PURE FUNCTION gradGs(r0, rp, k) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0, rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp) :: R, RR
    COMPLEX (KIND=dp), DIMENSION(3) :: res

    R = normr(r0-rp)

    IF(R==0) THEN
       res(:) = 0.0_dp
    ELSE
       RR = R*R
       res = oofpi*(EXP((0,1)*k*R)*(1.0_dp/(R*RR) -&
            (0,1)*k/(RR)) - 1.0_dp/(R*RR) - (k**2)/(2.0_dp*R))*(r0-rp)
    END IF
  END FUNCTION gradGs

  SUBROUTINE vgradGs(r0, rp, k, N, res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0
    REAL (KIND=dp), DIMENSION(3,N), INTENT(IN) :: rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    INTEGER, INTENT(IN) :: N
    REAL (KIND=dp), DIMENSION(N) :: R, RR
    COMPLEX (KIND=dp), DIMENSION(N) :: phasor
    COMPLEX (KIND=dp), DIMENSION(3,N), INTENT(INOUT) :: res

    R = SQRT((r0(1)-rp(1,:))**2 + (r0(2)-rp(2,:))**2 + (r0(3)-rp(3,:))**2)

    WHERE(R==0)
       res(1,:) = 0.0_dp
       res(2,:) = 0.0_dp
       res(3,:) = 0.0_dp
    ELSEWHERE
       RR = R*R
       phasor = oofpi*(EXP((0,1)*k*R)*(1.0_dp/(R*RR) -&
            (0,1)*k/(RR)) - 1.0_dp/(R*RR) - (k**2)/(2.0_dp*R))
       res(1,:) = phasor*(r0(1)-rp(1,:))
       res(2,:) = phasor*(r0(2)-rp(2,:))
       res(3,:) = phasor*(r0(3)-rp(3,:))
    END WHERE
  END SUBROUTINE vgradGs

  PURE FUNCTION Gf(r0, rp, k) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0, rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp) :: R
    COMPLEX (KIND=dp) :: res

    R = normr(r0-rp)

    res = oofpi*EXP((0,1)*k*R)/R
  END FUNCTION Gf

  SUBROUTINE vGf(r0, rp, k, N, res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0
    REAL (KIND=dp), DIMENSION(3,N), INTENT(IN) :: rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    INTEGER, INTENT(IN) :: N
    COMPLEX (KIND=dp), DIMENSION(N), INTENT(INOUT) :: res

    ! Store distances temporarily to result.
    res = SQRT((r0(1)-rp(1,:))**2 + (r0(2)-rp(2,:))**2 + (r0(3)-rp(3,:))**2)

    res = oofpi*EXP((0,1)*k*res)/res
  END SUBROUTINE vGf

  PURE FUNCTION gradGf(r0, rp, k) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0, rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp) :: R, RR
    COMPLEX (KIND=dp), DIMENSION(3) :: res

    R = normr(r0-rp)
    RR = R*R

    res = (r0-rp)*oofpi*EXP((0,1)*k*R)*(1.0_dp/(R*RR) - (0,1)*k/(RR))
  END FUNCTION gradGf

  SUBROUTINE vgradGf(r0, rp, k, N, res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0
    REAL (KIND=dp), DIMENSION(3,N), INTENT(IN) :: rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    INTEGER, INTENT(IN) :: N
    REAL (KIND=dp), DIMENSION(N) :: R, RR
    COMPLEX (KIND=dp), DIMENSION(N) :: phasor
    COMPLEX (KIND=dp), DIMENSION(3,N), INTENT(INOUT) :: res

    R = SQRT((r0(1)-rp(1,:))**2 + (r0(2)-rp(2,:))**2 + (r0(3)-rp(3,:))**2)
    RR = R*R
    
    phasor = oofpi*EXP((0,1)*k*R)*(1.0_dp/(R*RR) - (0,1)*k/(RR))
    res(1,:) = phasor*(r0(1)-rp(1,:))
    res(2,:) = phasor*(r0(2)-rp(2,:))
    res(3,:) = phasor*(r0(3)-rp(3,:))
  END SUBROUTINE vgradGf

  ! Green's function, where the first term of Taylor series
  ! is subtracted. To be used with Mueller method, where
  ! this term is cancelled.
  PURE FUNCTION Gs1(r0, rp, k) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0, rp
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp) :: R
    COMPLEX (KIND=dp) :: res, x
    INTEGER :: n

    R = normr(r0-rp)

    IF(R==0) THEN
       res = k*iofpi
    !ELSE IF((ABS(k)*R)<0.1_dp) THEN
    !   res = 0.0_dp
    !   x = (0,1)*k*R
    !   DO n=7,1,-1
    !      res = res + expcoef(n)*(x**n)
    !   END DO
    !   res = res/(4.0_dp*pi*R)
    ELSE
       res = oofpi*(EXP((0,1)*k*R)-1.0_dp)/R
    END IF
  END FUNCTION Gs1

  SUBROUTINE vGs1(r0, rp, k, N, res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r0
    REAL (KIND=dp), DIMENSION(3,N), INTENT(IN) :: rp
    INTEGER, INTENT(IN) :: N
    COMPLEX (KIND=dp), INTENT(IN) :: k
    COMPLEX (KIND=dp), DIMENSION(N), INTENT(INOUT) :: res

    res = SQRT((r0(1)-rp(1,:))**2 + (r0(2)-rp(2,:))**2 + (r0(3)-rp(3,:))**2)

    WHERE(res==0)
       res = k*iofpi
    ELSEWHERE
       res = oofpi*(EXP((0,1)*k*res)-1.0_dp)/res
    END WHERE
  END SUBROUTINE vGs1
END MODULE green
