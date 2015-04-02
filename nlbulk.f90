! MODULE: nlbulk
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Bulk nonlinear sources.
MODULE nlbulk
  USE nfields
  
  IMPLICIT NONE

  ! Properties of material bulk nonlinearity.
  TYPE medium_nlb
     COMPLEX (KIND=dp) :: delta_prime
     COMPLEX (KIND=dp) :: gamma
     !COMPLEX (KIND=dp) :: chi2zzz ! deprecated

     ! Maps (Ex**2,Ey**2,Ez**2,EyEz,ExEz,ExEy) to (Px,Py,Pz) in crystal frame.
     COMPLEX (KIND=dp), DIMENSION(3,6) :: chi2

     ! T maps from crystal frame to lab frame. invT = inverse(T).
     REAL (KIND=dp), DIMENSION(3,3) :: T, invT
  END type medium_nlb

CONTAINS
  ! Computes the action on nonlinear bulk dipole polarization Pnlb.
  FUNCTION Pnlb_dipole_ga(mesh, nedgestot, omegaff, ri, x, ga, nlb, r, qd) RESULT(Pnlb)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: omegaff
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    TYPE(medium_nlb), INTENT(IN) :: nlb
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp), DIMENSION(3,SIZE(ga),SIZE(x,3)) :: Pnlb, e, h
    COMPLEX (KIND=dp), DIMENSION(3) :: e2, Pnlb2
    COMPLEX (KIND=dp) :: Pnlbz
    INTEGER :: na, ns

    CALL scat_fields_ga(mesh, ga, x, nedgestot, omegaff, ri, NULL(), r, qd, e, h)

    DO na=1,SIZE(ga)
       DO ns=1,SIZE(x,3)
          ! Map e from lab frame to crystal frame.
          e2 = MATMUL(nlb%invT, e(:,na,ns))

          ! Compute polarization in crystal frame.
          Pnlb2 = eps0*MATMUL(nlb%chi2, (/e2(1)**2, e2(2)**2, e2(3)**2,&
               e2(2)*e2(3), e2(1)*e2(3), e2(1)*e2(2)/))

          ! Map Pnlb to lab frame.
          Pnlb(:,na,ns) = MATMUL(nlb%T, Pnlb2)

          !Pnlbz = eps0*(nlb%chi2zzz)*(e(3,na,ns)**2)
          !Pnlb(:,na,ns) = (/0.0_dp, 0.0_dp, 1.0_dp/)*Pnlbz
       END DO
    END DO
  END FUNCTION Pnlb_dipole_ga

  ! Computes the excitation vector for a second-order dipolar bulk source for
  ! all group representations.
  SUBROUTINE srcvec_nlbulk_dipole(mesh, nedgestot, omegaff, riff, rish, xff, ga, nlb,&
       qd_tri, qd_tetra, src)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp), INTENT(IN) :: omegaff
    COMPLEX (KIND=dp), INTENT(IN) :: riff, rish
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: xff
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(medium_nlb), INTENT(IN) :: nlb
    TYPE(quad_data), INTENT(IN) :: qd_tri, qd_tetra

    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: src

    REAL (KIND=dp) :: omegash, vol, thresholdsq
    COMPLEX (KIND=dp), DIMENSION(3,qd_tetra%num_nodes,SIZE(ga),SIZE(xff,3)) :: Jb
    COMPLEX (KIND=dp), DIMENSION(3) :: jacJb
    REAL (KIND=dp), DIMENSION(3) ::rp, diff
    COMPLEX (KIND=dp) :: int1, int2, int3, c1, c2, k, cgae
    INTEGER :: nedges, nweights, m, n, p, q, t, na, na2, nr, index, ns, progress
    REAL (KIND=dp), DIMENSION(3,qd_tetra%num_nodes) :: qp
    COMPLEX (KIND=dp), DIMENSION(3,3,qd_tetra%num_nodes,mesh%nfaces,SIZE(ga),SIZE(ga)) :: intaux1,&
         intaux2, intaux3
    LOGICAL :: near

    nedges = mesh%nedges
    nweights = qd_tetra%num_nodes

    thresholdsq = (mesh%avelen*3)**2

    ! Second-harmonic frequency.
    omegash = 2.0_dp*omegaff
    k = rish*omegash/c0

    c1 = (0,1)*omegash*mu0
    c2 = 1.0_dp/((0,1)*omegash*eps0*(rish**2))

    src(:,:,:) = 0.0_dp

    progress = 0

    DO m=1,mesh%nsolids
       
       qp = quad_tetra_points(qd_tetra, m, mesh)
       vol = mesh%solids(m)%volume

       ! Pre-compute field integrals in parallel.
       DO na=1,SIZE(ga)
          DO na2=1,SIZE(ga)
             
             !$OMP PARALLEL DEFAULT(NONE)&
             !$OMP SHARED(mesh,nweights,ga,intaux1,intaux2,intaux3,k,na,na2,qp,thresholdsq,qd_tri)&
             !$OMP PRIVATE(n,t,rp,diff,near)
             !$OMP DO SCHEDULE(STATIC)
             DO n=1,mesh%nfaces
                DO t=1,nweights
                   rp = MATMUL(ga(na2)%j, qp(:,t))

                   diff = rp - MATMUL(ga(na)%j, mesh%faces(n)%cp)
                   IF(SUM(diff*diff)<thresholdsq) THEN
                      near = .TRUE.
                   ELSE
                      near = .FALSE.
                   END IF
                   
                   intaux1(:,:,t,n,na2,na) = intK2(rp, n, mesh, k, ga(na), NULL(), near, qd_tri)
                   intaux2(:,:,t,n,na2,na) = intK3(rp, n, mesh, k, ga(na), NULL(), near, qd_tri)
                   intaux3(:,:,t,n,na2,na) = intK4(rp, n, mesh, k, ga(na), -1, NULL(), near, qd_tri)
                END DO
             END DO
             !$OMP END DO
             !$OMP END PARALLEL
          END DO
       END DO

       ! Pre-compute the current sources for all representations.

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,Jb,omegash,mesh,nedgestot,omegaff,riff,xff,ga,nlb,qp,qd_tri)&
       !$OMP PRIVATE(t)
       !$OMP DO SCHEDULE(STATIC)
       DO t=1,nweights
          Jb(:,t,:,:) = -(0,1)*omegash*Pnlb_dipole_ga(mesh, nedgestot, omegaff,&
               riff, xff, ga, nlb, qp(:,t), qd_tri)
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(xff,nweights,omegash,mesh,nedgestot,omegaff,riff,ga,nlb,qp,qd_tetra,intaux1,intaux2,intaux3,c1,c2,vol,src,m,nedges,Jb)&
       !$OMP PRIVATE(ns,t,na,na2,nr,cgae,n,q,int1,int2,int3,index)
       !$OMP DO SCHEDULE(STATIC)
       DO ns=1,SIZE(xff,3)
          DO na=1,SIZE(ga)
             DO na2=1,SIZE(ga)
                   
                DO n=1,mesh%nfaces
                   
                   DO q=1,3
                      
                      int1 = 0.0_dp
                      int2 = 0.0_dp
                      int3 = 0.0_dp
                      
                      DO t=1,nweights
                         int1 = int1 + qd_tetra%weights(t)*dotc(Jb(:,t,na2,ns),&
                              intaux1(:,q,t,n,na2,na))
                         int2 = int2 + qd_tetra%weights(t)*dotc(Jb(:,t,na2,ns),&
                              intaux2(:,q,t,n,na2,na))
                         int3 = int3 + qd_tetra%weights(t)*dotc(Jb(:,t,na2,ns),&
                              intaux3(:,q,t,n,na2,na))
                      END DO
                      
                      int1 = c1*int1*vol
                      int2 = c2*int2*vol
                      int3 = int3*vol
                      
                      index = mesh%faces(n)%edge_indices(q)

                      DO nr=1,SIZE(ga)
                         cgae = CONJG(ga(na)%ef(nr))
                         
                         src(index,nr,ns) = src(index,nr,ns) + cgae*(int1 - int2)
                         src(index + nedges,nr,ns) = src(index + nedges,nr,ns) &
                              + cgae*ga(na)%detj*int3
                      END DO
                   END DO
                   
                END DO
             END DO
             
          END DO
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       IF(progress>mesh%nsolids/10) THEN
          WRITE(*,'(A,I0,A)') 'Computed ', NINT(100*REAL(m)/mesh%nsolids), ' percent of sources'
          progress = 0
       END IF

       progress = progress + 1

    END DO

    src(:,:,:) = src(:,:,:)/REAL(SIZE(ga))
  END SUBROUTINE srcvec_nlbulk_dipole

END MODULE nlbulk
