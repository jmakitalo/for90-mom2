! MODULE: dipole
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines to calculate the excitation vector due to a near-field dipole.
MODULE dipole
  USE srcint
  USE aux
  USE nfields

  IMPLICIT NONE

CONTAINS
  SUBROUTINE srcvec_dipole(mesh, nedgestot, omega, ri, ga, dpos, dmom, q)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: dpos
    COMPLEX (KIND=dp), DIMENSION(3), INTENT(IN) :: dmom
    REAL (KIND=dp), INTENT(IN) :: omega

    COMPLEX (KIND=dp), DIMENSION(:,:) :: q

    INTEGER :: m, p, index, nf, ns
    COMPLEX (KIND=dp), DIMENSION(3,3) :: int1, int2
    COMPLEX (KIND=dp) :: k, c1, c2

    TYPE(prdnfo), POINTER :: prd

    prd => NULL()

    k = ri*omega/c0
    c1 = (0,1)*omega*mu0
    c2 = 1.0_dp/((0,1)*omega*eps0*(ri**2))

    q(:,:) = 0.0_dp

    ! Assumes only identity representation.
    nf = 1
    ns = 1

    DO m=1,mesh%nfaces
       int1 = intK2(dpos, m, mesh, k, ga(ns), prd, .TRUE.)
       int2 = intK3(dpos, m, mesh, k, ga(ns), prd, .TRUE.)

       DO p=1,3
          index = mesh%faces(m)%edge_indices(p)
          index = mesh%edges(index)%parent_index
        
          q(index,nf) = q(index,nf) + dotc(dmom, int1(:,p)*c1 + int2(:,p)*c2)
       END DO
    END DO

    DO m=1,mesh%nfaces
       int1 = intK4(dpos, m, mesh, k, ga(ns), -1, prd, .TRUE.)

       DO p=1,3
          index = mesh%faces(m)%edge_indices(p)
          index = mesh%edges(index)%parent_index
        
          q(nedgestot + index,nf) = q(nedgestot + index,nf) - dotc(dmom, int1(:,p))
       END DO
    END DO
  END SUBROUTINE srcvec_dipole

  SUBROUTINE dipoleff(r, dpos, dmom, omega, ri, E, H)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r, dpos
    COMPLEX (KIND=dp), DIMENSION(3), INTENT(IN) :: dmom
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(OUT) :: E, H
    COMPLEX (KIND=dp) :: k, eta
    REAL (KIND=dp) :: rm
    REAL (KIND=dp), DIMENSION(3) :: dir

    k = ri*omega/c0
    eta = eta0/ri

    rm = normr(r)
    dir = r/rm

    E = (0,1)*omega*mu0*EXP((0,1)*k*rm)/(4*pi*rm)*EXP(-(0,1)*k*dotr(dpos,dir))*dmom

    H = crossc(CMPLX(r-dpos,KIND=dp),E)/(rm*eta)
  END SUBROUTINE dipoleff

  ! Computes the polarizability alpha of a particle. The dipole moment p
  ! is given by p = epsilon0*alpha*Einc, where Einc is the incident field at the location
  ! of the dipole. It is assumed that Einc = 1.
  FUNCTION polarizability(mesh, ga, x, nedgestot, omega, ri, prd, a) RESULT(alpha)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega, a
    INTEGER, INTENT(IN) :: nedgestot
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: x

    COMPLEX (KIND=dp), DIMENSION(3) :: alpha
    INTEGER, PARAMETER :: ntheta = 80
    INTEGER, PARAMETER :: nphi = 160
    INTEGER :: n, nt, np
    REAL (KIND=dp), DIMENSION(ntheta) :: thetaw, thetap
    REAL (KIND=dp), DIMENSION(nphi) :: phiw, phip
    REAL (KIND=dp), DIMENSION(3,3) :: map
    REAL (KIND=dp), DIMENSION(3) :: r
    COMPLEX (KIND=dp), DIMENSION(3,1) :: e, h
    COMPLEX (KIND=dp), DIMENSION(3) :: integ
    COMPLEX (KIND=dp) :: f, k

    k = ri*omega/c0

    CALL get_simpsons_weights(0.0_dp, PI, ntheta, thetaw)
    CALL get_simpsons_points(0.0_dp, PI, ntheta, thetap)

    CALL get_simpsons_weights(0.0_dp, 2*PI, nphi, phiw)
    CALL get_simpsons_points(0.0_dp, 2*PI, nphi, phip)

    integ(:) = 0.0_dp

    ! x,y,z components of polarizability.
    DO n=1,3

       ! Integration over spherical surface.
       DO nt=1,ntheta
          DO np=1,nphi
             r = get_dir(thetap(nt), phip(np))*a

             ! Expressions compute the z-component naturally.
             ! For x and y components a rotation mapping is used.
             IF(n==1) THEN
                map(1,:) = (/0,0,1/)
                map(2,:) = (/0,1,0/)
                map(3,:) = (/-1,0,0/)
             ELSE IF(n==2) THEN
                map(1,:) = (/1,0,0/)
                map(2,:) = (/0,0,1/)
                map(3,:) = (/0,-1,0/)
             ELSE
                map = id33r
             END IF

             CALL scat_fields(mesh, ga, x, nedgestot, omega, ri, prd,&
                  MATMUL(map, r), e, h)

             ! Map the evaluated fields.
             e = MATMUL(CMPLX(TRANSPOSE(map),KIND=dp), e)
             h = MATMUL(CMPLX(TRANSPOSE(map),KIND=dp), h)

             integ(n) = integ(n) +&
                  thetaw(nt)*phiw(np)*COS(thetap(nt))*SIN(thetap(nt))*dotc(CMPLX(r,KIND=dp), e(:,1))
          END DO
       END DO
    END DO

    ! Coefficient of spherical harmonic Y_10 (z-dipole).
    integ = integ*SQRT(3/(4*pi))

    f = SQRT(4*pi/3)*EXP((0,1)*k*a)*(1/(a**2) - (0,1)*k/a)/(2*pi*eps0)

    alpha = integ/(f*eps0)

  END FUNCTION polarizability
END MODULE dipole
