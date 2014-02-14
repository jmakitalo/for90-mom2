! MODULE: quad
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Contains weights and points for mainly Gauss-Legendre quadrature over
! triangles and tetrahedra.
MODULE quad
  USE mesh

  IMPLICIT NONE

  REAL (KIND=dp), DIMENSION(4), PARAMETER :: qw4 = (/-0.56250000_dp,&
       0.52083333_dp,0.52083333_dp,0.52083333_dp/)

  REAL (KIND=dp), DIMENSION(4,3), PARAMETER :: quadFactors4 = RESHAPE((/&
       1.0_dp/3.0_dp,&
       0.6_dp,&
       0.2_dp,&
       0.2_dp,&

       1.0_dp/3.0_dp,&
       0.2_dp,&
       0.6_dp,&
       0.2_dp,&

       1.0_dp/3.0_dp,&
       0.2_dp,&
       0.2_dp,&
       0.6_dp/), (/4,3/))

  REAL (KIND=dp), DIMENSION(7), PARAMETER :: qw7 = (/0.22500000_dp,&
       0.13239415_dp, 0.13239415_dp, 0.13239415_dp, 0.12593918_dp,&
       0.12593918_dp, 0.12593918_dp/)

  REAL (KIND=dp), DIMENSION(7,3), PARAMETER :: quadFactors7 = RESHAPE((/&
       0.33333333_dp,&
       0.05971587_dp,&
       0.47014206_dp,&
       0.47014206_dp,&
       0.79742698_dp,&
       0.10128650_dp,&
       0.10128650_dp,&

       0.33333333_dp,&
       0.47014206_dp,&
       0.05971587_dp,&
       0.47014206_dp,&
       0.10128650_dp,&
       0.79742698_dp,&
       0.10128650_dp,&

       0.33333333_dp,&
       0.47014206_dp,&
       0.47014206_dp,&
       0.05971587_dp,&
       0.10128650_dp,&
       0.10128650_dp,&
       0.79742698_dp&
       /), (/7,3/))

  REAL (KIND=dp), DIMENSION(6), PARAMETER :: qw6 = (/&
       0.109951743655322_dp,&
       0.109951743655322_dp,&
       0.109951743655322_dp,&
       0.223381589678011_dp,&
       0.223381589678011_dp,&
       0.223381589678011_dp/)

  REAL (KIND=dp), DIMENSION(6,3), PARAMETER :: quadFactors6 = RESHAPE((/&
       0.816847572980459_dp,&
       0.091576213509771_dp,&
       0.091576213509771_dp,&
       0.108103018168070_dp,&
       0.445948490915965_dp,&
       0.445948490915965_dp,&

       0.091576213509771_dp,&
       0.816847572980459_dp,&
       0.091576213509771_dp,&
       0.445948490915965_dp,&
       0.108103018168070_dp,&
       0.445948490915965_dp,&

       0.091576213509771_dp,&
       0.091576213509771_dp,&
       0.816847572980459_dp,&
       0.445948490915965_dp,&
       0.445948490915965_dp,&
       0.108103018168070_dp/),(/6,3/))

  REAL (KIND=dp), DIMENSION(3), PARAMETER :: qw3 = (/&
       0.333333333333333,&
       0.333333333333333,&
       0.333333333333333/)

  REAL (KIND=dp), DIMENSION(3,3), PARAMETER :: quadFactors3 = RESHAPE((/&
       0.666666666666667,&
       0.166666666666667,&
       0.166666666666667,&

       0.166666666666667,&
       0.666666666666667,&
       0.166666666666667,&

       0.166666666666667,&
       0.166666666666667,&
       0.666666666666667/),(/3,3/))

  REAL (KIND=dp), DIMENSION(13), PARAMETER :: qw13 = (/&
       -0.149570044467670,&
       0.175615257433204,&
       0.175615257433204,&
       0.175615257433204,&
       0.053347235608839,&
       0.053347235608839,&
       0.053347235608839,&
       0.077113760890257,&
       0.077113760890257,&
       0.077113760890257,&
       0.077113760890257,&
       0.077113760890257,&
       0.077113760890257/)

  REAL (KIND=dp), DIMENSION(13,3), PARAMETER :: quadFactors13 = RESHAPE((/&
       0.333333333333333,&
       0.479308067841923,&
       0.260345966079038,&
       0.260345966079038,&
       0.869739794195568,&
       0.065130102902216,&
       0.065130102902216,&
       0.638444188569809,&
       0.638444188569809,&
       0.312865496004875,&
       0.312865496004875,&
       0.048690315425316,&
       0.048690315425316,&

       0.333333333333333,&
       0.260345966079038,&
       0.260345966079038,&
       0.479308067841923,&
       0.065130102902216,&
       0.065130102902216,&
       0.869739794195568,&
       0.048690315425316,&
       0.312865496004875,&
       0.048690315425316,&
       0.638444188569809,&
       0.312865496004875,&
       0.638444188569809,&

       0.333333333333333,&
       0.260345966079038,&
       0.479308067841923,&
       0.260345966079038,&
       0.065130102902216,&
       0.869739794195568,&
       0.065130102902216,&
       0.312865496004875,&
       0.048690315425316,&
       0.638444188569809,&
       0.048690315425316,&
       0.638444188569809,&
       0.312865496004875/),(/13,3/))

  REAL (KIND=dp), DIMENSION(13), PARAMETER :: qw = qw13
  REAL (KIND=dp), DIMENSION(13,3), PARAMETER :: quadFactors = quadFactors13
!  REAL (KIND=dp), DIMENSION(4), PARAMETER :: qw = qw4
!  REAL (KIND=dp), DIMENSION(4,3), PARAMETER :: quadFactors = quadFactors4

  REAL (KIND=dp), DIMENSION(1), PARAMETER :: volQw1 = (/1.0_dp/)

  REAL (KIND=dp), DIMENSION(1,4), PARAMETER :: volQuadFactors1 = RESHAPE((/&
       0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp/), (/1,4/))

  REAL (KIND=dp), DIMENSION(4), PARAMETER :: volQw4 = (/&
       0.25_dp,&
       0.25_dp,&
       0.25_dp,&
       0.25_dp/)

  REAL (KIND=dp), DIMENSION(4,4), PARAMETER :: volQuadFactors4 = RESHAPE((/&
       0.585410196624969,&
       0.138196601125011,&
       0.138196601125011,&
       0.138196601125011,&

       0.138196601125011,&
       0.585410196624969,&
       0.138196601125011,&
       0.138196601125011,&

       0.138196601125011,&
       0.138196601125011,&
       0.585410196624969,&
       0.138196601125011,&

       0.138196601125011,&
       0.138196601125011,&
       0.138196601125011,&
       0.585410196624969/),(/4,4/))

  REAL (KIND=dp), DIMENSION(11), PARAMETER :: volQw11 = (/&
       -0.013155555555556,&

       0.007622222222222,&
       0.007622222222222,&
       0.007622222222222,&
       0.007622222222222,&

       0.024888888888889,&
       0.024888888888889,&
       0.024888888888889,&
       0.024888888888889,&
       0.024888888888889,&
       0.024888888888889/)

  REAL (KIND=dp), DIMENSION(11,4), PARAMETER :: volQuadFactors11 = RESHAPE((/&
       0.250000000000000,&
       0.785714285714286,&
       0.071428571428571,&
       0.071428571428571,&
       0.071428571428571,&
       0.399403576166799,&
       0.100596423833201,&
       0.399403576166799,&
       0.100596423833201,&
       0.399403576166799,&
       0.100596423833201,&

       0.250000000000000,&
       0.071428571428571,&
       0.785714285714286,&
       0.071428571428571,&
       0.071428571428571,&
       0.399403576166799,&
       0.100596423833201,&
       0.100596423833201,&
       0.399403576166799,&
       0.100596423833201,&
       0.399403576166799,&

       0.250000000000000,&
       0.071428571428571,&
       0.071428571428571,&
       0.785714285714286,&
       0.071428571428571,&
       0.100596423833201,&
       0.399403576166799,&
       0.399403576166799,&
       0.100596423833201,&
       0.100596423833201,&
       0.399403576166799,&

       0.250000000000000,&
       0.071428571428571,&
       0.071428571428571,&
       0.071428571428571,&
       0.785714285714286,&
       0.100596423833201,&
       0.399403576166799,&
       0.100596423833201,&
       0.399403576166799,&
       0.399403576166799,&
       0.100596423833201/),(/11,4/))

  REAL (KIND=dp), DIMENSION(4), PARAMETER :: volQw = volQw4
  REAL (KIND=dp), DIMENSION(4,4), PARAMETER :: volQuadFactors = volQuadFactors4
!  REAL (KIND=dp), DIMENSION(1), PARAMETER :: volQw = volQw1
!  REAL (KIND=dp), DIMENSION(1,4), PARAMETER :: volQuadFactors = volQuadFactors1

CONTAINS
  FUNCTION GLquad_points(faceind, mesh) RESULT(res)
    INTEGER, INTENT(IN) :: faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: res
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, p3
    INTEGER :: n

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
    p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
    p3 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p

    DO n=1,SIZE(qw)
       res(:,n) = quadFactors(n,3)*p1 + quadFactors(n,1)*p2 + quadFactors(n,2)*p3
    END DO
  END FUNCTION GLquad_points

  FUNCTION GLquad_points_ext(faceind, mesh, fact, dim) RESULT(res)
    INTEGER, INTENT(IN) :: faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: dim
    REAL (KIND=dp), DIMENSION(dim,3), INTENT(IN) :: fact
    REAL (KIND=dp), DIMENSION(3,dim) :: res
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, p3
    INTEGER :: n

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
    p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
    p3 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p

    DO n=1,dim
       res(:,n) = fact(n,3)*p1 + fact(n,1)*p2 + fact(n,2)*p3
    END DO
  END FUNCTION GLquad_points_ext

  FUNCTION GLsolid_quad_points(solidind, mesh) RESULT(res)
    INTEGER, INTENT(IN) :: solidind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3,SIZE(volQw)) :: res
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, p3, p4
    INTEGER :: n

    p1 = mesh%nodes(mesh%solids(solidind)%node_indices(1))%p
    p2 = mesh%nodes(mesh%solids(solidind)%node_indices(2))%p
    p3 = mesh%nodes(mesh%solids(solidind)%node_indices(3))%p
    p4 = mesh%nodes(mesh%solids(solidind)%node_indices(4))%p

    DO n=1,SIZE(volQw)
       res(:,n) = volQuadFactors(n,1)*p1 + volQuadFactors(n,2)*p2 + volQuadFactors(n,3)*p3 + volQuadFactors(n,4)*p4
    END DO
  END FUNCTION GLsolid_quad_points

  ! n is the number of subintervals -> number of nodes is n+1.
  ! n must be even.
  SUBROUTINE get_simpsons_weights(a, b, n, weights)
    REAL (KIND=dp), INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: n
    REAL (KIND=dp), DIMENSION(n+1), INTENT(OUT) :: weights
    REAL (KIND=dp) :: h

    IF(MOD(n,2)/=0) THEN
       WRITE(*,*) "Number of subintervals for Simpson's rule must be even!"
       STOP
    END IF

    h = (b-a)/n

    weights((/1,n+1/)) = h/3.0_dp
    weights(3:(n-1):2) = 2.0_dp*h/3.0_dp
    weights(2:n:2) = 4.0_dp*h/3.0_dp
  END SUBROUTINE get_simpsons_weights

  ! n is the number of subintervals -> number of nodes is n+1.
  ! n must be even.
  SUBROUTINE get_simpsons_points(a, b, n, points)
    REAL (KIND=dp), INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: n
    REAL (KIND=dp), DIMENSION(n+1), INTENT(OUT) :: points
    REAL (KIND=dp) :: h
    INTEGER :: i

    IF(MOD(n,2)/=0) THEN
       WRITE(*,*) "Number of subintervals for Simpson's rule must be even!"
       STOP
    END IF

    h = (b-a)/n

    points(1:(n+1)) = a + h*(/(i,i=0,n)/)
  END SUBROUTINE get_simpsons_points

  RECURSIVE FUNCTION asqz_aux(f, a, b, eps, s, fa, fb, fc, level, maxDepth) RESULT(res)
    COMPLEX (KIND=dp), EXTERNAL :: f
    REAL (KIND=dp), INTENT(IN) :: a, b, eps
    COMPLEX (KIND=dp), INTENT(IN) :: s, fa, fb, fc
    INTEGER, INTENT(IN) :: level, maxDepth

    COMPLEX (KIND=dp) :: res, fd, fe, sleft, sright, s2
    REAL (KIND=dp) :: c, h, d, e

    c = (a + b)/2
    h = b - a

    d = (a + c)/2
    e = (c + b)/2

    fd = f(d)
    fe = f(e)

    sleft = (h/12)*(fa + 4*fd + fc)
    sright = (h/12)*(fc + 4*fe + fb)

    s2 = sleft + sright

    !IF(bottom<=0 .OR. (ABS(s2 - s)<=15*eps .AND. bottom<4) ) THEN
    IF(level>=maxDepth .OR. (ABS(s2 - s)<=15*eps*ABS(s2) .AND. level>1) ) THEN
       res = s2 + (s2 - s)/15
    ELSE
       res = asqz_aux(f, a, c, eps/2, sleft, fa, fc, fd, level+1, maxDepth) +&
            asqz_aux(f, c, b, eps/2, sright, fc, fb, fe, level+1, maxDepth)
    END IF
  END FUNCTION asqz_aux

  ! Adaptive Simpson's method based on listing at
  ! http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method#C
  FUNCTION asqz(f, a, b, eps, maxDepth) RESULT(res)
    COMPLEX (KIND=dp), EXTERNAL :: f
    REAL (KIND=dp), INTENT(IN) :: a, b, eps
    INTEGER, INTENT(IN) :: maxDepth

    COMPLEX (KIND=dp) :: fa, fb, fc, s, res
    REAL (KIND=dp) :: c, h

    c = (a + b)/2
    h = b - a

    fa = f(a)
    fb = f(b)
    fc = f(c)

    s = (h/6)*(fa + 4*fc + fb)

    res = asqz_aux(f, a, b, eps, s, fa, fb, fc, 0, maxDepth)
  END FUNCTION asqz

  ! Integrates f(x,y) over [x1,x2]x[y1,y2].
  SUBROUTINE asqz2(f, x1, x2, y1, y2, eps, maxDepth, res)
    COMPLEX (KIND=dp), EXTERNAL :: f
    REAL (KIND=dp), INTENT(IN) :: x1, x2, y1, y2, eps
    INTEGER, INTENT(IN) :: maxDepth
    COMPLEX (KIND=dp), INTENT(INOUT) :: res

    ! Auxiliary variable to make f(x,y) to appear a
    ! function of single variable for two nested
    ! integration routines.
    REAL (KIND=dp) :: gy

    res = asqz(fnested, y1, y2, eps, maxDepth)

  CONTAINS
    ! Evaluates f with y fixed to global value gy.
    FUNCTION fproxy(x) RESULT(z)
      REAL (KIND=dp), INTENT(IN) :: x
      COMPLEX (KIND=dp) :: z

      z = f(x,gy)
    END FUNCTION fproxy

    FUNCTION fnested(y) RESULT(z)
      REAL (KIND=dp), INTENT(IN) :: y
      COMPLEX (KIND=dp) :: z
      
      gy = y

      z = asqz(fproxy, x1, x2, eps, maxDepth)
    END FUNCTION fnested

  END SUBROUTINE asqz2
END MODULE quad
