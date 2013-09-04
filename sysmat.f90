! MODULE: sysmat
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for calculating approximate matrix representations of
! three boundary integral operators in space spanned by RWG basis.
! Also contains routines for calculating PMCHWT system matrix from
! these operator representations.
MODULE sysmat
  USE srcint
  USE common

  IMPLICIT NONE

CONTAINS
  ! Computes the PMCHWT system matrix for a problem described by
  ! struct b and for wavelength of given index wlind. Matrices are
  ! computed for all group representations.
  SUBROUTINE sysmat_pmchwt(b, wlind, A)
    TYPE(batch), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: wlind

    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: A
    INTEGER :: N, nf, ns, nd, nind
    INTEGER, DIMENSION(1:b%mesh%nedges) :: ind
    REAL (KIND=dp) :: omega, detj
    COMPLEX (KIND=dp) :: Zsi, ri, gae
    TYPE(prdnfo), POINTER :: prd

    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: M

    N = b%mesh%nedges

    omega = 2.0_dp*pi*c0/b%sols(wlind)%wl
    
    A(:,:,:) = 0.0_dp

    ALLOCATE(M(N,N))

    ! Loop through domains.
    DO nd=1,SIZE(b%domains)
       ri = b%media(b%domains(nd)%medium_index)%prop(wlind)%ri
       Zsi = (ri**2)/(eta0**2)

       ! Determine global edges indices to index basis functions.
       nind = b%domains(nd)%mesh%nedges
       ind(1:nind) = b%domains(nd)%mesh%edges(:)%parent_index

       ! Choose Green's function.
       IF(b%domains(nd)%gf_index==-1) THEN
          prd => NULL()
       ELSE
          prd => b%prd(b%domains(nd)%gf_index)
       END IF

       ! Loop through group actions.
       DO ns=1,SIZE(b%ga)
          detj = b%ga(ns)%detj
       
          ! Compute block matrix using action of index ns.
          ! Matrices are independent of field actions.

          CALL computeD(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, M(1:nind,1:nind))

          ! Add the block matrix to different sub-problems multiplied by proper
          ! field actions.
          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A(ind(1:nind),ind(1:nind),nf) = A(ind(1:nind),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A((ind(1:nind)+N),(ind(1:nind)+N),nf) = A((ind(1:nind)+N),(ind(1:nind)+N),nf)&
                  + gae*detj*Zsi*M(1:nind,1:nind)
          END DO

          CALL computeH(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, M(1:nind,1:nind))

          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A(ind(1:nind),ind(1:nind),nf) = A(ind(1:nind),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A((ind(1:nind)+N),(ind(1:nind)+N),nf) = A((ind(1:nind)+N),(ind(1:nind)+N),nf)&
                  + gae*detj*Zsi*M(1:nind,1:nind)
          END DO

          CALL computeK(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, M(1:nind,1:nind))

          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A((ind(1:nind)+N),ind(1:nind),nf) = A((ind(1:nind)+N),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A(ind(1:nind),(ind(1:nind)+N),nf) = A(ind(1:nind),(ind(1:nind)+N),nf)&
                  - gae*detj*M(1:nind,1:nind)
          END DO

       END DO
    END DO

    DEALLOCATE(M)
  END SUBROUTINE sysmat_pmchwt

  ! Computes the PMCHWT system matrix for second-harmonic problem.
  ! Adds products of src_coef and the submatrices to src_vec.
  ! Note that omega and ri must correspond to the SH case.
  ! Also the electric field group actions (gae) are the squares of
  ! actions in the linear problem.
  SUBROUTINE sysmat_pmchwt_nls(b, wlind, A, src_coef, src_vec)
    TYPE(batch), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: wlind
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: src_coef
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: src_vec

    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: A
    INTEGER :: N, nf, ns, nd, nind
    INTEGER, DIMENSION(1:b%mesh%nedges) :: ind
    REAL (KIND=dp) :: omega, detj
    COMPLEX (KIND=dp) :: Zsi, ri, gae
    TYPE(prdnfo), POINTER :: prd

    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: M

    N = b%mesh%nedges

    ! Second-harmonic frequency.
    omega = 4.0_dp*pi*c0/b%sols(wlind)%wl
    
    A(:,:,:) = 0.0_dp

    ALLOCATE(M(N,N))

    ! Loop through domains.
    DO nd=1,SIZE(b%domains)
       ri = b%media(b%domains(nd)%medium_index)%prop(wlind)%shri
       Zsi = (ri**2)/(eta0**2)

       ! Determine global edges indices to index basis functions.
       nind = b%domains(nd)%mesh%nedges
       ind(1:nind) = b%domains(nd)%mesh%edges(:)%parent_index

       ! Choose Green's function.
       IF(b%domains(nd)%gf_index==-1) THEN
          prd => NULL()
       ELSE
          prd => b%prd(b%domains(nd)%gf_index)
       END IF

       ! Loop through actions.
       DO ns=1,SIZE(b%ga)
          detj = b%ga(ns)%detj
       
          ! Compute block matrix using action of index ns.
          ! Matrices are independent of field actions.

          CALL computeD(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, M(1:nind,1:nind))

          ! Add the block matrix to different sub-problems multiplied by proper
          ! field actions.
          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)**2

             A(ind(1:nind),ind(1:nind),nf) = A(ind(1:nind),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A((ind(1:nind)+N),(ind(1:nind)+N),nf) = A((ind(1:nind)+N),(ind(1:nind)+N),nf)&
                  + gae*detj*Zsi*M(1:nind,1:nind)

             src_vec(ind(1:nind),nf) = src_vec(ind(1:nind),nf) - &
                  gae*MATMUL(M(1:nind,1:nind), src_coef((nind+1):(2*nind),nd,nf))

             src_vec((ind(1:nind)+N),nf) = src_vec((ind(1:nind)+N),nf) - &
                  gae*detj*MATMUL(M(1:nind,1:nind), src_coef(1:nind,nd,nf))*Zsi
          END DO

          CALL computeH(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, M(1:nind,1:nind))

          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)**2

             A(ind(1:nind),ind(1:nind),nf) = A(ind(1:nind),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A((ind(1:nind)+N),(ind(1:nind)+N),nf) = A((ind(1:nind)+N),(ind(1:nind)+N),nf)&
                  + gae*detj*Zsi*M(1:nind,1:nind)

             src_vec(ind(1:nind),nf) = src_vec(ind(1:nind),nf) - &
                  gae*MATMUL(M(1:nind,1:nind), src_coef((nind+1):(2*nind),nd,nf))

             src_vec((ind(1:nind)+N),nf) = src_vec((ind(1:nind)+N),nf) - &
                  gae*detj*MATMUL(M(1:nind,1:nind), src_coef(1:nind,nd,nf))*Zsi
          END DO

          CALL computeK(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, M(1:nind,1:nind))

          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)**2

             A((ind(1:nind)+N),ind(1:nind),nf) = A((ind(1:nind)+N),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A(ind(1:nind),(ind(1:nind)+N),nf) = A(ind(1:nind),(ind(1:nind)+N),nf)&
                  - gae*detj*M(1:nind,1:nind)

             src_vec(ind(1:nind),nf) = src_vec(ind(1:nind),nf) + &
                  gae*MATMUL(M(1:nind,1:nind), src_coef(1:nind,nd,nf))

             src_vec((ind(1:nind)+N),nf) = src_vec((ind(1:nind)+N),nf) - &
                  gae*MATMUL(M(1:nind,1:nind), src_coef((nind+1):(2*nind),nd,nf))
          END DO

       END DO
    END DO

    DEALLOCATE(M)
  END SUBROUTINE sysmat_pmchwt_nls

  ! Determine if faces of indices n and m are considered near,
  ! so as to use singularity subtraction, which is time-consuming.
  FUNCTION near_faces(mesh, prd, n, m) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    INTEGER, INTENT(IN) :: n, m

    LOGICAL :: res
    REAL (KIND=dp), DIMENSION(3) :: diff
    REAL (KIND=dp) :: distsq, threshold

    threshold = mesh%avelen*3

    diff = mesh%faces(m)%cp - mesh%faces(n)%cp
    distsq = SUM(diff*diff)

    IF(distsq<threshold**2) THEN
       res = .TRUE.
    ELSE IF(ASSOCIATED(prd)) THEN
       IF(prd%type==prd_2d .AND. distsq>(MINVAL((/prd%dx,prd%dy/))-threshold)**2) THEN
          res = .TRUE.
       ELSE
          res = .FALSE.
       END IF
    ELSE
       res = .FALSE.
    END IF

  END FUNCTION near_faces

  ! Computes the D-matrix elements
  ! D_mng = -i*omega_mu*int_Sm dS fm(r). int_Sn' dS' J_g*fn(r')*O'_gG(r,r').
  ! IN:
  ! omega: The angular frequency.
  ! ri: Refractive index of the domain.
  ! mesh: Surface mesh of the scatterer.
  ! ga: Group action identifiers.
  ! prd: The periodic Green's function.
  ! OUT:
  ! D: The D-matrix.
  SUBROUTINE computeD(omega, ri, mesh, ga, prd, D)
    ! Input variables
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd

    ! Internal variables
    COMPLEX (KIND=dp), INTENT(INOUT), DIMENSION(mesh%nedges,mesh%nedges) :: D
    COMPLEX (KIND=dp) :: c1, k
    INTEGER :: nweights, n, m, p, q, r, index1, index2
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: qpm
    COMPLEX (KIND=dp) :: int1
    COMPLEX (KIND=dp), DIMENSION(3,3,SIZE(qw),mesh%nfaces) :: intaux
    REAL (KIND=dp), DIMENSION(3,SIZE(qw),3) :: fmv
    LOGICAL :: near

    WRITE(*,*) 'Building a D-matrix'

    nweights = SIZE(qw)

    D(:,:) = 0.0_dp

    k = ri*omega/c0

    ! Coefficients of partial integrals.
    c1 = -(0,1)*omega*mu0

    DO m=1,mesh%nfaces

       qpm = GLquad_points(m, mesh)
       Am = mesh%faces(m)%area

       DO p=1,3
          CALL vrwg(qpm,m,p,mesh,fmv(:,:,p))
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intaux,qpm,mesh,k,ga,prd,m)&
       !$OMP PRIVATE(n,near,r)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m)

          DO r=1,nweights
             intaux(:,:,r,n) = intK2(qpm(:,r), n, mesh, k, ga, prd, near)
          END DO                   
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       DO n=1,mesh%nfaces

          DO p=1,3
             DO q=1,3

                int1 = 0.0_dp
                DO r=1,nweights
                   int1 = int1 + qw(r)*dotc(CMPLX(fmv(:,r,p),KIND=dp), intaux(:,q,r,n))
                END DO
                int1 = int1*Am

                index1 = mesh%faces(m)%edge_indices(p)
                index2 = mesh%faces(n)%edge_indices(q)

                D(index1,index2) = D(index1,index2) + c1*int1
             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE computeD

  ! Computes the H-matrix elements
  ! H_mng = -1/(i*omega*epsilon)*int_Sm dS div_S(fm(r)) int_Sn' dS' div'_S(fn(r'))*O'_gG(r,r').
  ! IN:
  ! omega: The angular frequency.
  ! ri: Refractive index of the domain.
  ! mesh: Surface mesh of the scatterer.
  ! ga: Group action identifiers.
  ! prd: The periodic Green's function.
  ! OUT:
  ! H: The H-matrix.
  SUBROUTINE computeH(omega, ri, mesh, ga, prd, H)
    ! Input variables
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd

    ! Internal variables
    COMPLEX (KIND=dp), INTENT(INOUT), DIMENSION(mesh%nedges,mesh%nedges) :: H
    COMPLEX (KIND=dp) :: c2, k
    INTEGER :: nweights, n, m, p, q, r, index1, index2
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: qpm
    COMPLEX (KIND=dp) :: int1
    REAL (KIND=dp), DIMENSION(3) :: fmDiv
    COMPLEX (KIND=dp), DIMENSION(3,SIZE(qw),mesh%nfaces) :: intaux
    LOGICAL :: near

    WRITE(*,*) 'Building an H-matrix'

    nweights = SIZE(qw)

    H(:,:) = 0.0_dp

    k = ri*omega/c0

    ! Coefficients of partial integrals.
    c2 = -1.0_dp/((0,1)*omega*eps0*(ri**2))

    DO m=1,mesh%nfaces

       qpm = GLquad_points(m, mesh)
       Am = mesh%faces(m)%area

       DO p=1,3
          fmDiv(p) = rwgDiv(m, p, mesh)
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intaux,qpm,mesh,k,ga,prd,m)&
       !$OMP PRIVATE(r,n,near)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m)

          DO r=1,nweights
             intaux(:,r,n) = intK1(qpm(:,r), n, mesh, k, ga, prd, near)
          END DO
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
    
       DO n=1,mesh%nfaces

          DO p=1,3
             DO q=1,3

                int1 = 0.0_dp
                DO r=1,nweights
                   int1 = int1 + qw(r)*intaux(q,r,n)
                END DO
                int1 = int1*Am*fmDiv(p)

                index1 = mesh%faces(m)%edge_indices(p)
                index2 = mesh%faces(n)%edge_indices(q)

                H(index1,index2) = H(index1,index2) + c2*int1
             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE computeH

  ! Computes the K-matrix elements
  ! K_mng = int_Sm dS fm(r). int_Sn' dS' [O'_g grad'G(r,r')]x(J_g*fn(r')).
  ! IN:
  ! omega: The angular frequency.
  ! ri: Refractive index of the domain.
  ! mesh: Surface mesh of the scatterer.
  ! ga: Group action identifiers.
  ! prd: The periodic Green's function.
  ! OUT:
  ! A: The K-matrix.
  SUBROUTINE computeK(omega, ri, mesh, ga, prd, A)
    ! Input variables
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd

    ! Internal variables
    COMPLEX (KIND=dp), INTENT(INOUT), DIMENSION(mesh%nedges,mesh%nedges) :: A
    COMPLEX (KIND=dp) :: k, int1
    INTEGER :: nweights, n, m, p, q, r, index1, index2
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: qpm
    COMPLEX (KIND=dp), DIMENSION(3,3,SIZE(qw),mesh%nfaces) :: intaux
    REAL (KIND=dp), DIMENSION(3,SIZE(qw),3) :: fmv
    LOGICAL :: near

    WRITE(*,*) 'Building a K-matrix'

    nweights = SIZE(qw)

    A(:,:) = 0.0_dp

    k = ri*omega/c0

    DO m=1,mesh%nfaces

       qpm = GLquad_points(m, mesh)
       Am = mesh%faces(m)%area

       DO p=1,3
          CALL vrwg(qpm,m,p,mesh,fmv(:,:,p))
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intaux,qpm,mesh,k,ga,m,prd)&
       !$OMP PRIVATE(r,n,near)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m)

          DO r=1,nweights
             intaux(:,:,r,n) = intK4(qpm(:,r), n, mesh, k, ga, m, prd, near)
          END DO
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
    
       DO n=1,mesh%nfaces


          DO q=1,3
             DO p=1,3

                int1 = 0.0_dp
                DO r=1,nweights
                   int1 = int1 + qw(r)*dotc(CMPLX(fmv(:,r,p),KIND=dp), intaux(:,q,r,n))
                END DO
                int1 = int1*Am

                index1 = mesh%faces(m)%edge_indices(p)
                index2 = mesh%faces(n)%edge_indices(q)

                A(index1,index2) = A(index1,index2) + int1

             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE computeK
END MODULE sysmat
