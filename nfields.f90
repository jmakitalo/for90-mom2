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
  ! If multidomain description is used, domain is the domain index referencing
  ! b%domains. Otherwise domain==1 corresponds to exterior and domain==2 to interior
  ! domain.
  SUBROUTINE scat_fields(mesh, ga, x, nedgestot, omega, ri, prd, r, e, h)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    INTEGER, INTENT(IN) :: nedgestot
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(INOUT) :: e, h
    COMPLEX (KIND=dp), DIMENSION(3) :: e2, h2
    INTEGER :: nf

    e(:) = 0.0_dp
    h(:) = 0.0_dp

    DO nf=1,SIZE(ga)
       CALL scat_fields_frag(mesh, ga, nf, x(:,nf), nedgestot, omega, ri, prd, r, e, h)

       e = e + e2
       h = h + h2
    END DO
  END SUBROUTINE scat_fields

  SUBROUTINE scat_fields_frag(mesh, ga, nf, x, nedgestot, omega, ri, prd, r, e, h)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nf, nedgestot
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: x
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(INOUT) :: e, h

    INTEGER :: n, q, index, ns
    COMPLEX (KIND=dp) :: c1, c2, k, gae
    COMPLEX (KIND=dp), DIMENSION(3,3) :: int1, int2, int3
    COMPLEX (KIND=dp), DIMENSION(nedgestot) :: alpha, beta

    e(:) = 0.0_dp
    h(:) = 0.0_dp

    k = ri*omega/c0

    alpha = x(1:nedgestot)
    beta = x((nedgestot+1):(2*nedgestot))

    ! Coefficients of partial integrals.
    c1 = (0,1)*omega*mu0
    c2 = 1.0_dp/((0,1)*omega*eps0*(ri**2))
    
    DO n=1,mesh%nfaces

       DO ns=1,SIZE(ga)
          int1 = intK2(r, n, mesh, k, ga(ns), prd, .TRUE.)
          int2 = intK3(r, n, mesh, k, ga(ns), prd, .TRUE.)
          int3 = intK4(r, n, mesh, k, ga(ns), 0, prd, .TRUE.)
       
          DO q=1,3
             index = mesh%faces(n)%edge_indices(q)
             index = mesh%edges(index)%parent_index
          
             gae = ga(ns)%ef(nf)
             
             e = e + alpha(index)*gae*(c1*int1(:,q) + c2*int2(:,q)) +&
                  gae*ga(ns)%detj*beta(index)*int3(:,q)

             h = h + beta(index)*gae*ga(ns)%detj*(int1(:,q)/c2 + int2(:,q)/c1) -&
                  gae*alpha(index)*int3(:,q)
          END DO
       END DO
    END DO
  END SUBROUTINE scat_fields_frag

  SUBROUTINE field_mesh(name, mesh, scale, nedgestot, x, ga, omega, ri)
    CHARACTER (LEN=*), INTENT(IN) :: name
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: scale, omega
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    COMPLEX (KIND=dp), INTENT(IN) :: ri

    INTEGER :: n, n2, q, index, nf, nga, na
    COMPLEX (KIND=dp), DIMENSION(3,mesh%nfaces*SIZE(ga)) :: ef, hf
    COMPLEX (KIND=dp), DIMENSION(3) :: et, ht
    COMPLEX (KIND=dp) :: en, hn
    COMPLEX (KIND=dp) :: eps, gae
    REAL (KIND=dp), DIMENSION(3) :: fn, nor
    REAL (KIND=dp) :: detj
    CHARACTER (LEN=256) :: oname, numstr
    TYPE(mesh_container) :: mesh2

    WRITE(*,*) 'Computing near fields on particle mesh.'

    eps = (ri**2)*eps0

    nga = SIZE(ga)

    mesh2%nnodes = mesh%nnodes*nga
    mesh2%nfaces = mesh%nfaces*nga
    ALLOCATE(mesh2%nodes(mesh2%nnodes))
    ALLOCATE(mesh2%faces(mesh2%nfaces))

    DO na=1,nga
       DO n=1,mesh%nnodes
          mesh2%nodes(n + mesh%nnodes*(na-1))%p = MATMUL(ga(na)%j, mesh%nodes(n)%p)
       END DO

       DO n=1,mesh%nfaces
          mesh2%faces(n + mesh%nfaces*(na-1))%node_indices(:) = &
               mesh%faces(n)%node_indices(:) + mesh%nnodes*(na-1)

          mesh2%faces(n + mesh%nfaces*(na-1))%n = MATMUL(ga(na)%j, mesh%faces(n)%n)
       END DO
    END DO

    DO na=1,nga
       detj = ga(na)%detj

       DO n=1,mesh%nfaces

          n2 = n + mesh%nfaces*(na-1)

          et(:) = 0.0_dp
          en = 0.0_dp

          ht(:) = 0.0_dp
          hn = 0.0_dp

          nor = MATMUL(ga(na)%j, mesh%faces(n)%n)
          
          DO nf=1,nga
             gae = ga(na)%ef(nf)

             DO q=1,3
                index = mesh%faces(n)%edge_indices(q)
                index = mesh%edges(index)%parent_index

                fn = MATMUL(ga(na)%j, crossr(mesh%faces(n)%n, rwg(mesh%faces(n)%cp, n, q, mesh)))

                ht = ht - fn*x(index, nf)*gae*detj
                hn = hn + rwgDiv(n, q, mesh)*x(nedgestot + index, nf)*gae*detj
                
                et = et + fn*x(nedgestot + index, nf)*gae
                en = en + rwgDiv(n, q, mesh)*x(index, nf)*gae
             END DO
          END DO

          en = en/((0,1)*omega*eps)
          hn = hn/((0,1)*omega*mu0)

          ef(:,n2) = et + nor*en
          hf(:,n2) = ht + nor*hn
          
       END DO

    END DO

    CALL save_vector_fields_msh(name, mesh2, ef, hf, scale)

    WRITE(*,*) 'Maximum of |E| is ',&
         MAXVAL(SQRT(ABS(ef(1,:))**2 + ABS(ef(2,:))**2 + ABS(ef(3,:))**2))


    !CALL save_field_msh(name, mesh2, en, eta, scale)

    DEALLOCATE(mesh2%nodes, mesh2%faces)

  END SUBROUTINE field_mesh
END MODULE nfields
