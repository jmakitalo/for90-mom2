! MODULE: bessel
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implements the Bessel function.
MODULE bessel
  USE constants

  IMPLICIT NONE

CONTAINS
  Subroutine besselj(m,n,e,x,y)
    !Labels: 100, 200, 210
    real*8 a,b0,b1,b2,e,x,y
    integer i,m,n
    a = 1.d0
    if (n <= 1) goto 100
    ! Calculate N! 
    do i = 1, n
       a = a * dfloat(i)
    end do
100 a = 1 / a
    if (n.eq.0) goto 200
    ! Calculate multiplying term
    do i = 1, n
       a = a * x / 2.d0
    end do
200 b0 = 1.d0
    b2 = 1.d0
    m = 0
    !Assemble series sum
210 m = m + 1
    b1 = -(x * x * b0) / (m * (m + n) * 4.d0)
    b2 = b2 + b1
    b0 = b1
    ! Test for convergence
    if (dabs(b1) > e)  goto 210
    ! form final answer
    y = a * b2
    return
  end Subroutine besselj

  FUNCTION factorial(n) RESULT(res)
    INTEGER, INTENT(IN) :: n
    REAL (KIND=dp) :: res
    INTEGER :: i

    res = 1.0_dp

    DO i=2,n
       res = res*REAL(i,KIND=dp)
    END DO
  END FUNCTION factorial

  ! Return Bessel function of first kind and integer order nu>=0 for
  ! complex argument z. The result is only accurate for small |z|
  ! and nu.
  FUNCTION zbessj(z, nu) RESULT(res)
    COMPLEX (KIND=dp), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: nu
    COMPLEX (KIND=dp) :: res, z2

    INTEGER, PARAMETER :: maxn = 20
    INTEGER :: n

    z2 = z/2.0_dp
    
    res = 0.0_dp

    DO n=0,maxn
       res = res + ((-1)**n)*(z2**(2*n+nu))/( factorial(n)*factorial(n+nu) )
    END DO

  END FUNCTION zbessj

END MODULE bessel
