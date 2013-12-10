! MODULE: cs
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for computing scattering and absoprtion cross-sections.
MODULE cs
  USE source

  IMPLICIT NONE

CONTAINS
  SUBROUTINE cs_prtsrf(mesh, nedgestot, omega, ri, ga, x, src, csca, cabs)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(srcdata), INTENT(IN) :: src

    REAL (KIND=dp), INTENT(OUT) :: csca, cabs

    REAL (KIND=dp) :: Sn
    INTEGER :: n, nweights, p, t, ns, nf, ns2, index
    COMPLEX (KIND=dp), DIMENSION(3) :: Et, Ht, einc, hinc,&
         Eti, Hti, Ets, Hts, Et2, Ht2, einc2, hinc2
    REAL (KIND=dp) :: An, int1, int2
    REAL (KIND=dp), DIMENSION(3) :: fn, nor
    REAL (KIND=dp) :: c
    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: qpn
    COMPLEX (KIND=dp) :: gae

    nweights = SIZE(qw)

    c = 2.0_dp/(eps0*c0*REAL(ri,KIND=dp))

    csca = 0.0_dp
    cabs = 0.0_dp

    DO n=1,mesh%nfaces

       qpn = GLquad_points(n, mesh)
       An = mesh%faces(n)%area

       int1 = 0.0_dp
       int2 = 0.0_dp
       DO t=1,nweights

          DO ns=1,SIZE(ga)

             Et(:) = 0.0_dp
             Ht(:) = 0.0_dp
             Ets(:) = 0.0_dp
             Hts(:) = 0.0_dp

             nor = MATMUL(ga(ns)%j, mesh%faces(n)%n)

             ! Compute the action of O_g to nxH and nxE.
             DO nf=1,SIZE(ga)
                gae = ga(ns)%ef(nf)

                Et2(:) = 0.0_dp
                Ht2(:) = 0.0_dp
                DO p=1,3
                   fn = MATMUL(ga(ns)%j, crossr(mesh%faces(n)%n, rwg(qpn(:,t),n,p,mesh)))
                   index = mesh%faces(n)%edge_indices(p)
                   index = mesh%edges(index)%parent_index

                   Ht2 = Ht2 - x(index, nf)*fn*gae*ga(ns)%detj
                   Et2 = Et2 + x(index + nedgestot, nf)*fn*gae
                END DO

                Et = Et + Et2
                Ht = Ht + Ht2
             END DO

             ! Compute the incident fields.
             CALL src_fields(src, omega, ri, MATMUL(ga(ns)%j,qpn(:,t)), einc, hinc)

             ! Compute the scattered fields.
             Hti = crossc(crossc(CMPLX(nor,KIND=dp), hinc), CMPLX(nor,KIND=dp))
             Hts = Ht - Hti
                
             Eti = crossc(crossc(CMPLX(nor,KIND=dp), einc), CMPLX(nor,KIND=dp))
             Ets = Et - Eti

             Sn = 0.5_dp*REAL(dotc(crossc(Et, CONJG(Ht)), CMPLX(nor,KIND=dp)), KIND=dp)
             int1 = int1 + qw(t)*Sn

             Sn = 0.5_dp*REAL(dotc(crossc(Ets, CONJG(Hts)), CMPLX(nor,KIND=dp)), KIND=dp)
             int2 = int2 + qw(t)*Sn

          END DO
       END DO
       cabs = cabs + int1*An
       csca = csca + int2*An
    END DO

    cabs = -c*cabs
    csca = c*csca
  END SUBROUTINE cs_prtsrf

END MODULE cs
