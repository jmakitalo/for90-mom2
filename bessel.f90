! MODULE: bessel
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implements the Bessel function.
MODULE bessel
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
END MODULE bessel
