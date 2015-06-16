! MODULE: multipole
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines to project scattered near-fields to multipole moments.
MODULE multipole
  USE nfields

  IMPLICIT NONE

CONTAINS
  ! Calculates the total number of multipole moments A(n,m) if the principal
  ! index n goes up to nmp. The zeroth moment is omitted.
  FUNCTION mp_number(nmp) RESULT(res)
    INTEGER, INTENT(IN) :: nmp
    INTEGER :: n, res

    res = 0

    DO n=1,nmp
      res = res + 2*n + 1
    END DO

  END FUNCTION mp_number

  ! Calculates transverse-electric (TE) and transverse-magnetic (TM)
  ! multipole moments from solution x in the specified domain.
  ! mesh: boundary mesh of the domain.
  ! ga: symmetry groups actions.
  ! x: solution vector (columns: expansion coef, group rep, source).
  ! nedgestot: total number of edges in the problem geometry.
  ! omega: angular frequency.
  ! ri: refractive index in the domain.
  ! qd: quadrature routine data.
  ! nmp: number of multipoles (mnp>=1).
  ! d: distance from the origin at which the fields are evaluated.
  !    The scatterer must be enclosed by a sphere of this radius.
  ! te,tm: multipole coefficients divided by h^(1)_n(k*d) (spherical Hankel).
  !        These are linear arrays, where A(n,m) multipoles are placed as
  !        (A(1,-1),A(1,0),A(1,1),A(2,-2),A(2,-1),A(2,0),A(2,1),A(2,2),...).
  !        The total size is given by mp_number(nmp).
  SUBROUTINE mp_moments(mesh, ga, x, nedgestot, omega, ri, qd, nmp, d, te, tm)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega, d
    INTEGER, INTENT(IN) :: nedgestot, nmp
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: x
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: te, tm

    TYPE(prdnfo), POINTER :: prd
    INTEGER :: ind, n, m
    COMPLEX (KIND=dp) :: coef, eta

    INTEGER, PARAMETER :: mptype_tm = 1, mptype_te = 2

    ! Multipoles cannot be calculated for periodic scatterers.
    prd => NULL()

    eta = eta0/ri

    ind = 0

    DO n=1,nmp
      coef = (ri*omega/c0)/SQRT(n*(n+1.0_dp))

      DO m=-n,n
        CALL asqz2(int_te, 0.0_dp, pi, 0.0_dp, 2*pi, 1D-4, 10, te(ind))
        CALL asqz2(int_tm, 0.0_dp, pi, 0.0_dp, 2*pi, 1D-4, 10, tm(ind))

        te(ind) = te(ind)*coef
        tm(ind) = -tm(ind)*coef/eta

        ind = ind + 1
      END DO
    END DO

    CONTAINS
      FUNCTION int_te(theta, phi) RESULT(res)
        REAL (KIND=dp), INTENT(IN) :: theta, phi

        COMPLEX (KIND=dp) :: res

        res = int_moment(theta, phi, mptype_te)
      END FUNCTION int_te

      FUNCTION int_tm(theta, phi) RESULT(res)
        REAL (KIND=dp), INTENT(IN) :: theta, phi

        COMPLEX (KIND=dp) :: res

        res = int_moment(theta, phi, mptype_tm)
      END FUNCTION int_tm

      FUNCTION int_moment(theta, phi, mptype) RESULT(res)
        REAL (KIND=dp), INTENT(IN) :: theta, phi
        INTEGER, INTENT(IN) :: mptype

        COMPLEX (KIND=dp) :: res
        COMPLEX (KIND=dp), DIMENSION(3,SIZE(x,3),2) :: eh
        REAL (KIND=dp), DIMENSION(3) :: r
        REAL (KIND=dp) :: sp, cp, st, ct

        sp = SIN(phi)
        cp = COS(phi)
        st = SIN(theta)
        ct = COS(theta)

        r = (/st*cp, st*sp, ct/)*d

        CALL scat_fields(mesh, ga, x, nedgestot, omega, ri, prd, r, qd, eh(:,:,1), eh(:,:,2))

        res = CONJG(sph(ct, phi, n, m))*dotc(CMPLX(r,KIND=dp), eh(:,1,mptype))*st
      END FUNCTION int_moment

  END SUBROUTINE mp_moments

  ! Spherical harmonic function, Jackson, Eq. (3.53).
  FUNCTION sph(ct, phi, n, m) RESULT(res)
    REAL (KIND=dp), INTENT(IN) :: ct, phi
    INTEGER, INTENT(IN) :: n, m
    REAL (KIND=dp), DIMENSION(1,0:n) :: leg
    COMPLEX (KIND=dp) :: res

    CALL pmn_polynomial_value(1, n, m, (/ct/), leg)

    res = leg(1,n)*EXP((0,1)*m*phi)/SQRT(2*pi)
  END FUNCTION sph

  function r8_factorial ( n )

  !*****************************************************************************80
  !
  !! R8_FACTORIAL computes the factorial of N.
  !
  !  Discussion:
  !
  !    factorial ( N ) = product ( 1 <= I <= N ) I
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 January 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the argument of the factorial function.
  !    If N is less than 1, the function value is returned as 1.
  !
  !    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
  !
    implicit none

    real ( kind = 8 ) r8_factorial
    integer ( kind = 4 ) i
    integer ( kind = 4 ) n

    r8_factorial = 1.0D+00

    do i = 1, n
      r8_factorial = r8_factorial * real ( i, kind = 8 )
    end do

    return
  end function r8_factorial

  subroutine pmn_polynomial_value ( mm, n, m, x, cx )
  !*****************************************************************************80
  !
  !! PMN_POLYNOMIAL_VALUE: normalized Legendre polynomial Pmn(n,m,x).
  !
  !  Discussion:
  !
  !    The unnormalized associated Legendre functions P_N^M(X) have
  !    the property that
  !
  !      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX
  !      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
  !
  !    By dividing the function by the square root of this term,
  !    the normalized associated Legendre functions have norm 1.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 March 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Milton Abramowitz, Irene Stegun,
  !    Handbook of Mathematical Functions,
  !    National Bureau of Standards, 1964,
  !    ISBN: 0-486-61272-4,
  !    LC: QA47.A34.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) MM, the number of evaluation points.
  !
  !    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
  !    function, which must be at least 0.
  !
  !    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
  !    which must be at least 0, and no greater than N.
  !
  !    Input, real ( kind = 8 ) X(MM), the evaluation points.
  !
  !    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
  !
    implicit none

    integer ( kind = 4 ) mm
    integer ( kind = 4 ) n

    real ( kind = 8 ) cx(mm,0:n)
    real ( kind = 8 ) factor
    integer ( kind = 4 ) j
    integer ( kind = 4 ) m
    !real ( kind = 8 ) r8_factorial
    real ( kind = 8 ) x(mm)

    if ( m < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PMN_POLYNOMIAL_VALUE - Fatal error!'
      write ( *, '(a,i8)' ) '  Input value of M is ', m
      write ( *, '(a)' ) '  but M must be nonnegative.'
      stop 1
    end if

    if ( n < m ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PMN_POLYNOMIAL_VALUE - Fatal error!'
      write ( *, '(a,i8)' ) '  Input value of M = ', m
      write ( *, '(a,i8)' ) '  Input value of N = ', n
      write ( *, '(a)' ) '  but M must be less than or equal to N.'
      stop 1
    end if

    cx(1:mm,0:n) = 0.0D+00

    if ( m <= n ) then
      cx(1:mm,m) = 1.0D+00
      factor = 1.0D+00
      do j = 1, m
        cx(1:mm,m) = - cx(1:mm,m) * factor * sqrt ( 1.0D+00 - x(1:mm)**2 )
        factor = factor + 2.0D+00
      end do
    end if

    if ( m + 1 <= n ) then
      cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
    end if

    do j = m + 2, n
      cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                   + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &
                   / real (     j - m,     kind = 8 )
    end do
  !
  !  Normalization.
  !
    do j = m, n
      factor = sqrt ( ( real ( 2 * j + 1, kind = 8 ) * r8_factorial ( j - m ) ) &
        / ( 2.0D+00 * r8_factorial ( j + m ) ) )
      cx(1:mm,j) = cx(1:mm,j) * factor
    end do

    return
  end subroutine pmn_polynomial_value
END MODULE multipole
