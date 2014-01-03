MODULE nlbulk
  USE nfields
  
  IMPLICIT NONE

  TYPE medium_nlb
     COMPLEX (KIND=dp) :: delta_prime
     COMPLEX (KIND=dp) :: gamma
     COMPLEX (KIND=dp) :: chi2zzz
  END type medium_nlb

CONTAINS
  ! Computes the action on nonlinear bulk dipole polarization Pnlb.
  FUNCTION Pnlb_dipole_ga(mesh, nedgestot, omegaff, ri, x, ga, gai, nlb, r) RESULT(Pnlb)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: omegaff
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot, gai
    TYPE(medium_nlb), INTENT(IN) :: nlb
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3) :: Pnlb, e, h
    COMPLEX (KIND=dp) :: Pnlbz
    INTEGER :: nf

    CALL scat_fields_ga(mesh, ga, x, nedgestot, omegaff, ri, NULL(), gai, r, e, h)

    Pnlbz = eps0*(nlb%chi2zzz)*(e(3)**2)

    Pnlb = (/0.0_dp, 0.0_dp, 1.0_dp/)*Pnlbz
  END FUNCTION Pnlb_dipole_ga

  ! Computes the excitation vector for a second-order dipolar bulk source for
  ! all group representations.
  SUBROUTINE srcvec_nlbulk_dipole(mesh, nedgestot, omegaff, riff, rish, xff, ga, nlb, src)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp), INTENT(IN) :: omegaff
    COMPLEX (KIND=dp), INTENT(IN) :: riff, rish
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: xff
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(medium_nlb), INTENT(IN) :: nlb

    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: src

    REAL (KIND=dp) :: omegash, vol
    COMPLEX (KIND=dp), DIMENSION(3,SIZE(volQw),SIZE(ga),mesh%nsolids) :: Jb
    COMPLEX (KIND=dp), DIMENSION(3) :: jacJb
    REAL (KIND=dp), DIMENSION(3) ::rp
    COMPLEX (KIND=dp) :: int1, int2, int3, c1, c2, k, cgae
    INTEGER :: nedges, nweights, m, n, p, q, t, na, na2, nr, index
    REAL (KIND=dp), DIMENSION(3,SIZE(volQw)) :: qp
    COMPLEX (KIND=dp), DIMENSION(3,3,SIZE(volQw),mesh%nfaces) :: intaux1, intaux2, intaux3

    nedges = mesh%nedges
    nweights = SIZE(volQw)

    ! Second-harmonic frequency.
    omegash = 2.0_dp*omegaff
    k = rish*omegash/c0

    c1 = (0,1)*omegash*mu0
    c2 = 1.0_dp/((0,1)*omegash*eps0*(rish**2))

    src(:,:) = 0.0_dp

    ! Pre-compute the current sources for all representations.

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(mesh,nweights,ga,omegash,nedgestot,omegaff,riff,xff,nlb,Jb)&
    !$OMP PRIVATE(m,na,t,qp)
    !$OMP DO SCHEDULE(STATIC)
    DO m=1,mesh%nsolids

       qp = GLsolid_quad_points(m, mesh)

       DO na=1,SIZE(ga)
          DO t=1,nweights
             Jb(:,t,na,m) = -(0,1)*omegash*Pnlb_dipole_ga(mesh, nedgestot, omegaff,&
                  riff, xff, ga, na, nlb, qp(:,t))
          END DO
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    DO na=1,SIZE(ga)
       DO na2=1,SIZE(ga)
          DO m=1,mesh%nsolids

             qp = GLsolid_quad_points(m, mesh)
             vol = mesh%solids(m)%volume

             ! Pre-compute field integrals in parallel.
             
             !$OMP PARALLEL DEFAULT(NONE)&
             !$OMP SHARED(mesh,nweights,ga,intaux1,intaux2,intaux3,k,na,na2,qp)&
             !$OMP PRIVATE(n,t,rp)
             !$OMP DO SCHEDULE(STATIC)
             DO n=1,mesh%nfaces
                DO t=1,nweights
                   rp = MATMUL(ga(na2)%j, qp(:,t))
                   
                   intaux1(:,:,t,n) = intK2(rp, n, mesh, k, ga(na), NULL(), .TRUE.)
                   intaux2(:,:,t,n) = intK3(rp, n, mesh, k, ga(na), NULL(), .TRUE.)
                   intaux3(:,:,t,n) = intK4(rp, n, mesh, k, ga(na), -1, NULL(), .TRUE.)
                END DO
             END DO
             !$OMP END DO
             !$OMP END PARALLEL

             DO nr=1,SIZE(ga)
                
                cgae = CONJG(ga(na)%ef(nr))
                
                DO n=1,mesh%nfaces
                   
                   DO q=1,3
                      
                      int1 = 0.0_dp
                      int2 = 0.0_dp
                      int3 = 0.0_dp
                      
                      DO t=1,nweights
                         int1 = int1 + volQw(t)*dotc(Jb(:,t,na2,m), intaux1(:,q,t,n))
                         int2 = int2 + volQw(t)*dotc(Jb(:,t,na2,m), intaux2(:,q,t,n))
                         int3 = int3 + volQw(t)*dotc(Jb(:,t,na2,m), intaux3(:,q,t,n))
                      END DO
                      
                      int1 = c1*int1*vol
                      int2 = c2*int2*vol
                      int3 = int3*vol
                      
                      index = mesh%faces(n)%edge_indices(q)
                      
                      src(index,nr) = src(index,nr) + cgae*(int1 - int2)
                      src(index + nedges,nr) = src(index + nedges,nr) + cgae*ga(na)%detj*int3
                   END DO
                END DO

             END DO
          END DO

       END DO
    END DO

    src(:,:) = src(:,:)/REAL(SIZE(ga))
  END SUBROUTINE srcvec_nlbulk_dipole

END MODULE nlbulk
