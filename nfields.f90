! MODULE: nfields
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines to compute the scattered fields in the vicinity
! of scatterers.
MODULE nfields
  USE srcint

  IMPLICIT NONE

  TYPE nfield_plane
     REAL (KIND=dp), DIMENSION(3) :: origin, v1, v2
     REAL (KIND=dp) :: d1, d2
     INTEGER :: n1, n2
  END TYPE nfield_plane

CONTAINS
  SUBROUTINE scat_fields(mesh, ga, x, nedgestot, omega, ri, prd, r, e, h)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    INTEGER, INTENT(IN) :: nedgestot
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: x
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3,SIZE(x,3)), INTENT(INOUT) :: e, h
    COMPLEX (KIND=dp), DIMENSION(3,SIZE(ga),SIZE(x,3)) :: e2, h2
    INTEGER :: nf, m

    e(:,:) = 0.0_dp
    h(:,:) = 0.0_dp

    CALL scat_fields_frags(mesh, ga, x, nedgestot, omega, ri, prd, r, e2, h2)

    DO nf=1,SIZE(ga)
       DO m=1,SIZE(x,3)
          e(:,m) = e(:,m) + e2(:,nf,m)
          h(:,m) = h(:,m) + h2(:,nf,m)
       END DO
    END DO
  END SUBROUTINE scat_fields

  SUBROUTINE scat_fields_ga(mesh, ga, x, nedgestot, omega, ri, prd, r, e, h)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    INTEGER, INTENT(IN) :: nedgestot
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: x
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3,SIZE(ga),SIZE(x,3)), INTENT(INOUT) :: e, h
    COMPLEX (KIND=dp), DIMENSION(3,SIZE(ga),SIZE(x,3)) :: e2, h2
    COMPLEX (KIND=dp) :: gae
    INTEGER :: nf, na, ns

    e(:,:,:) = 0.0_dp
    h(:,:,:) = 0.0_dp

    CALL scat_fields_frags(mesh, ga, x, nedgestot, omega, ri, prd, r, e2, h2)

    DO na=1,SIZE(ga)
       DO nf=1,SIZE(ga)
          gae = ga(na)%ef(nf)
          
          DO ns=1,SIZE(x,3)
             e(:,na,ns) = e(:,na,ns) + gae*MATMUL(ga(na)%j, e2(:,nf,ns))
             h(:,na,ns) = h(:,na,ns) + ga(na)%detj*gae*MATMUL(ga(na)%j, h2(:,nf,ns))
          END DO
       END DO
    END DO
  END SUBROUTINE scat_fields_ga

!!$  SUBROUTINE scat_fields_invmap(mesh, ga, x, nedgestot, omega, ri, prd, gai, r, e, h)
!!$    TYPE(mesh_container), INTENT(IN) :: mesh
!!$    COMPLEX (KIND=dp), INTENT(IN) :: ri
!!$    REAL (KIND=dp), INTENT(IN) :: omega
!!$    INTEGER, INTENT(IN) :: nedgestot, gai
!!$    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
!!$    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
!!$    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
!!$    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
!!$
!!$    COMPLEX (KIND=dp), DIMENSION(3), INTENT(INOUT) :: e, h
!!$    COMPLEX (KIND=dp), DIMENSION(3,SIZE(ga)) :: e2, h2
!!$    COMPLEX (KIND=dp) :: cgae
!!$    INTEGER :: nf
!!$
!!$    e(:) = 0.0_dp
!!$    h(:) = 0.0_dp
!!$
!!$    CALL scat_fields_frags(mesh, ga, x, nedgestot, omega, ri, prd, r, e2, h2)
!!$
!!$    DO nf=1,SIZE(ga)
!!$       cgae = CONJG(ga(gai)%ef(nf))
!!$
!!$       e = e + cgae*MATMUL(TRANSPOSE(ga(gai)%j), e2(:,nf))
!!$       h = h + ga(gai)%detj*cgae*MATMUL(TRANSPOSE(ga(gai)%j), h2(:,nf))
!!$    END DO
!!$  END SUBROUTINE scat_fields_invmap

  SUBROUTINE scat_fields_frags(mesh, ga, x, nedgestot, omega, ri, prd, r, e, h)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: x
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3,SIZE(ga),SIZE(x,3)), INTENT(INOUT) :: e, h

    INTEGER :: n, q, index, ns, nf, m
    COMPLEX (KIND=dp) :: c1, c2, k, gae
    COMPLEX (KIND=dp), DIMENSION(3,3) :: int1, int2, int3
    REAL (KIND=dp), DIMENSION(3) :: diff
    REAL (KIND=dp) :: thresholdsq
    LOGICAL :: near

    e(:,:,:) = 0.0_dp
    h(:,:,:) = 0.0_dp

    k = ri*omega/c0

    thresholdsq = (mesh%avelen*3)**2

    ! Coefficients of partial integrals.
    c1 = (0,1)*omega*mu0
    c2 = 1.0_dp/((0,1)*omega*eps0*(ri**2))
    
    DO n=1,mesh%nfaces

       DO ns=1,SIZE(ga)

          diff = r - MATMUL(ga(ns)%j, mesh%faces(n)%cp)
          IF(SUM(diff*diff)<thresholdsq) THEN
             near = .TRUE.
          ELSE
             near = .FALSE.
          END IF

          int1 = intK2(r, n, mesh, k, ga(ns), prd, near)
          int2 = intK3(r, n, mesh, k, ga(ns), prd, near)
          int3 = intK4(r, n, mesh, k, ga(ns), 0, prd, near)
       
          DO q=1,3
             index = mesh%faces(n)%edge_indices(q)
             index = mesh%edges(index)%parent_index

             DO nf=1,SIZE(ga)
                
                gae = ga(ns)%ef(nf)
                
                DO m=1,SIZE(x,3)
                   e(:,nf,m) = e(:,nf,m) + x(index,nf,m)*gae*(c1*int1(:,q) + c2*int2(:,q)) +&
                        gae*ga(ns)%detj*x(index + nedgestot,nf,m)*int3(:,q)
                   
                   h(:,nf,m) = h(:,nf,m) + x(index + nedgestot,nf,m)*gae*ga(ns)%detj*&
                        (int1(:,q)/c2 + int2(:,q)/c1) -&
                        gae*x(index,nf,m)*int3(:,q)
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE scat_fields_frags
END MODULE nfields
