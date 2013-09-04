! MODULE: dipole
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines to calculate the excitation vector due to a near-field dipole.
MODULE dipole
  USE srcint

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
END MODULE dipole
