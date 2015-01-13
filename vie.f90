MODULE vie
  USE mesh
  USE rwgf
  USE quad
  USE green
  USE greenprd
  USE symmetry
  USE source

  IMPLICIT NONE

  ! Gauss-Legendre weights and nodes for 1-D quadrature over [-1,1].
  REAL (KIND=dp), DIMENSION(6), PARAMETER :: GL1Dw =&
       (/0.467913934572691_dp, 0.467913934572691_dp,&
       0.360761573048139_dp, 0.360761573048139_dp,&
       0.171324492379170_dp, 0.171324492379170_dp/)
  REAL (KIND=dp), DIMENSION(6), PARAMETER :: GL1Dn =&
       (/0.238619186083197_dp, -0.238619186083197_dp,&
       0.661209386466265_dp, -0.661209386466265_dp,&
       0.932469514203152_dp, -0.932469514203152_dp/)

  ! Index permutation for last three indices of 4 sub-tetrahedra.
  INTEGER, DIMENSION(3,4), PARAMETER :: subsolidind = (/1,2,3,  2,3,4, 3,4,1,  1,2,4/)

CONTAINS
  SUBROUTINE rcs_vie(mesh, k, ga, prd, qd_tetra, xi, x, ntheta, nphi, rcsdata)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd_tetra
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ntheta, nphi

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    REAL (KIND=dp), DIMENSION(1:ntheta,1:nphi), INTENT(OUT) :: rcsdata

    INTEGER :: n, m
    REAL (KIND=dp) :: theta, phi
    COMPLEX (KIND=dp), DIMENSION(3) :: scatamp, dir

    dir = get_dir(theta, phi)

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(ntheta,nphi,mesh,k,ga,prd,qd_tetra,x,rcsdata,dir)&
    !$OMP PRIVATE(n,m,theta,phi,scatamp)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,ntheta
       DO m=1,nphi
          theta = pi*REAL((n-1),KIND=dp)/REAL((ntheta-1),KIND=dp)
          phi = 2.0_dp*pi*REAL((m-1),KIND=dp)/REAL((nphi-1),KIND=dp)

          scatamp = scatamp_vie(mesh, k, ga, prd, qd_tetra, xi, x, theta, phi)

          rcsdata(n,m) = normc(crossc(dir, scatamp))/(4*pi)
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  END SUBROUTINE rcs_vie

  FUNCTION cext_vie(mesh, k, ga, prd, qd_tetra, xi, x, src) RESULT(cext)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd_tetra
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: x
    TYPE(srcdata), INTENT(IN) :: src

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    REAL (KIND=dp) :: cext
    COMPLEX (KIND=dp), DIMENSION(3) :: pol
    COMPLEX (KIND=dp), DIMENSION(3) :: scatamp

    scatamp = scatamp_vie(mesh, k, ga, prd, qd_tetra, xi, x, src%theta, src%phi)

    pol = get_pol(src%theta, src%phi, src%psi)

    cext = AIMAG(dotc(scatamp, pol))/REAL(k,KIND=dp)
  END FUNCTION cext_vie

  FUNCTION scatamp_vie(mesh, k, ga, prd, qd_tetra, xi, x, theta, phi) RESULT(scatamp)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd_tetra
    REAL (KIND=dp), INTENT(IN) :: theta, phi
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: x

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    COMPLEX (KIND=dp), DIMENSION(3) :: scatamp
    REAL (KIND=dp), DIMENSION(3,qd_tetra%num_nodes,4) :: fmv
    INTEGER :: m, p, r, nweights, index
    REAL (KIND=dp), DIMENSION(3,qd_tetra%num_nodes) :: qpm
    REAL (KIND=dp), DIMENSION(3) :: dir
    REAL (KIND=dp) :: Vm
    COMPLEX (KIND=dp) :: phasor
    COMPLEX (KIND=dp), DIMENSION(3) :: aux, int1

    dir = get_dir(theta, phi)

    nweights = qd_tetra%num_nodes

    scatamp(:) = 0.0_dp

    DO m=1,mesh%nsolids
       qpm = quad_tetra_points(qd_tetra, m, mesh)
       Vm = mesh%solids(m)%volume

       DO p=1,4
          CALL vsolid_rwg(qpm(:,:),m,p,mesh,fmv(:,:,p))
       END DO
       
       DO p=1,4
          index = mesh%solids(m)%solid_face_indices(p)

          int1 = 0.0_dp
          DO r=1,nweights
             phasor = EXP(-(0,1)*k*dotr(dir, qpm(:,r)))
             aux = MATMUL(xi(qpm(:,r),m), fmv(:,r,p))*x(index)

             int1 = int1 + qd_tetra%weights(r)*phasor*(aux - dir*dotc(CMPLX(dir,KIND=dp), aux))
          END DO
          
          scatamp = scatamp + Vm*int1
       END DO
    END DO

    scatamp = scatamp*(k**2)/eps0
  END FUNCTION scatamp_vie

  SUBROUTINE vie_matrix(mesh, k, ga, prd, qd_tri, qd_tetra, xi, A)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd_tri, qd_tetra

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: A
    TYPE(quad_data) :: qd_tetra_test

    qd_tetra_test = tetra_quad_data('tetra_gl1')

    ! Declare matrix to zero.
    A(:,:) = 0.0_dp

    ! Add Green operator elements.
    CALL vieGreenMatrix(mesh, k, ga, prd, qd_tri, qd_tetra, qd_tetra_test, .TRUE., xi, A)

    CALL delete_quad_data(qd_tetra_test)

    ! Move result to "LHS".
    A(:,:) = -A(:,:)

    ! Add identity operator elements.
    CALL vieIdMatrix(mesh, k, ga, prd, qd_tetra, xi, A)
  END SUBROUTINE vie_matrix

  SUBROUTINE vie_eigen(mesh, k, qd_tri, qd_tetra, epsr)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(quad_data), INTENT(IN) :: qd_tri, qd_tetra
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: epsr

    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: A, idmat, idmatInv, eigvec
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval
    INTEGER, DIMENSION(SIZE(epsr)) :: minind
    TYPE(quad_data) :: qd_tetra_test
    INTEGER :: dim
    TYPE(prdnfo), POINTER :: prd
    TYPE(group_action), DIMENSION(:), ALLOCATABLE :: ga

    CALL group_id(ga)

    prd => NULL()

    dim = mesh%nsolid_faces

    qd_tetra_test = tetra_quad_data('tetra_gl4')

    ALLOCATE(A(1:dim,1:dim), idmat(1:dim, 1:dim), idmatInv(1:dim, 1:dim))

    ! Declare matrix to zero.
    A(:,:) = 0.0_dp

    WRITE(*,*) 'Computing matrices'

    ! Add Green operator elements.
    CALL vieGreenMatrix(mesh, k, ga(1), prd, qd_tri, qd_tetra, qd_tetra_test, .TRUE., xi_id, A)

    CALL delete_quad_data(qd_tetra_test)

    idmat(:,:) = 0.0_dp

    ! Add identity operator elements.
    CALL vieIdMatrix(mesh, k, ga(1), prd, qd_tetra, xi_zero, idmat)

    ALLOCATE(eigvec(1:dim,1:dim), eigval(1:dim))


    WRITE(*,*) 'Computing eigenvalues'

    !CALL matrix_eigenvalues_gen(A, idmat, eigval, eigvec)
    CALL matrix_inverse(idmat, idmatInv)
    A = MATMUL(idmatInv, A)
    CALL matrix_eigenvalues(A, eigval, eigvec)

    ! eigval to epsr
    !eigval = 1.0_dp/eigval(:) + 1.0_dp
    !minind = find_smallest(eigval, dim, SIZE(epsr))
    !epsr = eigval(minind(:))
    epsr = 1.0_dp/eigval(:) + 1.0_dp

    DEALLOCATE(A, idmat, idmatInv, eigvec, eigval, ga)

  CONTAINS
    FUNCTION xi_id(pos, s) RESULT(xires)
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
      INTEGER, INTENT(IN) :: s
      COMPLEX (KIND=dp) :: eps, diag
      COMPLEX, DIMENSION(3,3) :: xires

      diag = 1.0_dp

      xires(:,:) = 0.0_dp

      xires(1,1) = diag
      xires(2,2) = diag
      xires(3,3) = diag
    END FUNCTION xi_id

    FUNCTION xi_zero(pos, s) RESULT(xires)
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
      INTEGER, INTENT(IN) :: s
      COMPLEX, DIMENSION(3,3) :: xires

      xires(:,:) = 0.0_dp
    END FUNCTION xi_zero

  END SUBROUTINE vie_eigen

  SUBROUTINE vie_eigen_spec()
    REAL (KIND=dp) :: scale, omega1, omega2, omega, t
    TYPE(quad_data) :: qd_tri, qd_tetra
    TYPE(mesh_container) :: mesh
    INTEGER :: nepsr, n
    COMPLEX (KIND=dp) :: k
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: epsr
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: data_re, data_im

    scale = 10D-9

    mesh = load_mesh('sphere.msh')
    
    CALL build_mesh(mesh, scale)

    IF(mesh%nsolid_faces==0) THEN
       CALL build_solid_faces(mesh)
       CALL compute_basis_data(mesh)
    END IF

    ! Quadrature rules.
    qd_tri = tri_quad_data('tri_gl4')
    qd_tetra = tetra_quad_data('tetra_gl1')

    nepsr = mesh%nsolid_faces
    omega = 2*pi*c0/6D-7

    ALLOCATE(epsr(1:nepsr), data_re(1:nepsr,1), data_im(1:nepsr,1))

    k = omega/c0
    
    CALL vie_eigen(mesh, k, qd_tri, qd_tetra, epsr)

    data_re(:,1) = REAL(epsr,KIND=dp)
    data_im(:,1) = AIMAG(epsr)

    CALL delete_quad_data(qd_tri)
    CALL delete_quad_data(qd_tetra)
    CALL delete_mesh(mesh)

    CALL write_data('epsr_re.dat', data_re)
    CALL write_data('epsr_im.dat', data_im)

    DEALLOCATE(epsr, data_re, data_im)
  END SUBROUTINE vie_eigen_spec

  SUBROUTINE vie_srcvec(mesh, omega, ri, ga, qd_tetra, src, b)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(quad_data), INTENT(IN) :: qd_tetra
    TYPE(srcdata), INTENT(IN) :: src

    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: b
    COMPLEX (KIND=dp), DIMENSION(3,qd_tetra%num_nodes) :: dinc
    REAL (KIND=dp), DIMENSION(3,qd_tetra%num_nodes,4) :: fmv
    INTEGER :: m, p, r, nweights, index
    REAL (KIND=dp), DIMENSION(3,qd_tetra%num_nodes) :: qpm
    REAL (KIND=dp) :: Vm
    COMPLEX (KIND=dp) :: int1
    COMPLEX (KIND=dp), DIMENSION(3) :: einc, hinc

    b(:) = 0.0_dp

    nweights = qd_tetra%num_nodes

    DO m=1,mesh%nsolids
       qpm = quad_tetra_points(qd_tetra, m, mesh)
       Vm = mesh%solids(m)%volume

       DO p=1,4
          CALL vsolid_rwg(qpm(:,:),m,p,mesh,fmv(:,:,p))
       END DO

       DO r=1,nweights
          CALL src_fields(src, omega, ri, qpm(:,r), einc, hinc)

          dinc(:,r) = eps0*einc
       END DO
       
       DO p=1,4
          int1 = 0.0_dp
          DO r=1,nweights
             int1 = int1 + qd_tetra%weights(r)*dotc(CMPLX(fmv(:,r,p),KIND=dp), dinc(:,r))
          END DO

          index = mesh%solids(m)%solid_face_indices(p)
          
          b(index) = b(index) + Vm*int1
       END DO
    END DO

  END SUBROUTINE vie_srcvec

  ! Does not zero matrix A at start.
  SUBROUTINE vieIdMatrix(mesh, k, ga, prd, qd_tetra, xi, A)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd_tetra

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: A
    REAL (KIND=dp), DIMENSION(3,qd_tetra%num_nodes,4) :: fmv
    INTEGER :: m, p, q, r, nweights, index1, index2
    REAL (KIND=dp), DIMENSION(3,qd_tetra%num_nodes) :: qpm
    REAL (KIND=dp) :: Vm
    COMPLEX (KIND=dp) :: int1

    nweights = qd_tetra%num_nodes

    DO m=1,mesh%nsolids
       qpm = quad_tetra_points(qd_tetra, m, mesh)
       Vm = mesh%solids(m)%volume

       DO p=1,4
          CALL vsolid_rwg(qpm(:,:),m,p,mesh,fmv(:,:,p))
       END DO
       
       DO p=1,4
          DO q=1,4
             
             int1 = 0.0_dp
             DO r=1,nweights
                int1 = int1 + qd_tetra%weights(r)*dotc(CMPLX(fmv(:,r,p),KIND=dp),&
                     MATMUL((id33r-xi(qpm(:,r),m)), CMPLX(fmv(:,r,q),KIND=dp)))
             END DO

             index1 = mesh%solids(m)%solid_face_indices(p)
             index2 = mesh%solids(m)%solid_face_indices(q)
             
             A(index1,index2) = A(index1,index2) + Vm*int1
          END DO
       END DO
    END DO
  END SUBROUTINE vieIdMatrix

  ! Does not zero matrix A at start.
  SUBROUTINE vieGreenMatrix(mesh, k, ga, prd, qd_tri, qd_tetra, qd_tetra_test, regularize, xi, A)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd_tri, qd_tetra, qd_tetra_test
    LOGICAL, INTENT(IN) :: regularize

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: A
    REAL (KIND=dp), DIMENSION(3) :: nor
    REAL (KIND=dp), DIMENSION(3,qd_tetra_test%num_nodes,4) :: fmv
    REAL (KIND=dp), DIMENSION(3,qd_tri%num_nodes,4) :: fmv_bnd
    INTEGER :: m, n, p, q, r, nweights, index1, index2, nweights_tri, bndfaceind, bnd
    COMPLEX (KIND=dp), DIMENSION(3,4,qd_tetra_test%num_nodes,mesh%nsolids) :: intaux1
    COMPLEX (KIND=dp), DIMENSION(4,qd_tetra_test%num_nodes,mesh%nsolids) :: intaux2
    COMPLEX (KIND=dp), DIMENSION(4,qd_tri%num_nodes,mesh%nsolids) :: intaux3
    COMPLEX (KIND=dp) :: ksq, int1, int2
    REAL (KIND=dp), DIMENSION(3,qd_tetra_test%num_nodes) :: qpm
    REAL (KIND=dp), DIMENSION(3,qd_tri%num_nodes) :: qpm_bnd
    REAL (KIND=dp) :: Vm, fmDiv, Am

    nweights = qd_tetra_test%num_nodes
    nweights_tri = qd_tri%num_nodes

    ksq = k**2

    DO m=1,mesh%nsolids
       qpm = quad_tetra_points(qd_tetra_test, m, mesh)
       Vm = mesh%solids(m)%volume

       DO p=1,4
          CALL vsolid_rwg(qpm(:,:),m,p,mesh,fmv(:,:,p))
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intaux1,intaux2,qpm,mesh,k,ga,prd,m,qd_tetra,regularize)&
       !$OMP PRIVATE(n,r)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nsolids

          IF(n==m .AND. regularize) THEN
             DO r=1,nweights
                intaux1(:,:,r,n) = intVolGSelf(qpm(:,r), n, mesh, k, ga, prd, qd_tetra, xi)

                intaux2(:,r,n) = intVolGradGSelf(qpm(:,r), n, mesh, k, ga, prd, qd_tetra, xi)
             END DO
          ELSE
             DO r=1,nweights
                intaux1(:,:,r,n) = intVolG(qpm(:,r), n, mesh, k, ga, prd, qd_tetra, xi)

                intaux2(:,r,n) = intVolGradG(qpm(:,r), n, mesh, k, ga, prd, qd_tetra, xi)
             END DO
          END IF
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       ! Volume-volume integrals.
       DO n=1,mesh%nsolids

          DO p=1,4
             fmDiv = solid_rwgDiv(m, p, mesh)

             DO q=1,4

                int1 = 0.0_dp
                int2 = 0.0_dp
                DO r=1,nweights
                   int1 = int1 + qd_tetra_test%weights(r)*dotc(CMPLX(fmv(:,r,p),KIND=dp), intaux1(:,q,r,n))

                   int2 = int2 + qd_tetra_test%weights(r)*fmDiv*intaux2(q,r,n)
                END DO

                index1 = mesh%solids(m)%solid_face_indices(p)
                index2 = mesh%solids(n)%solid_face_indices(q)

                A(index1,index2) = A(index1,index2) + Vm*(ksq*int1 + int2)
             END DO
          END DO
       END DO

       ! Boundary-volume integrals.
       DO bnd=1,4
          bndfaceind = mesh%solid_faces(mesh%solids(m)%solid_face_indices(bnd))%face_index
          IF(bndfaceind==-1) THEN
             CYCLE
          END IF

          qpm_bnd = quad_tri_points(qd_tri, bndfaceind, mesh)
          Am = mesh%faces(bndfaceind)%area
          nor = mesh%faces(bndfaceind)%n

          DO n=1,mesh%nsolids
             
             IF(n==m .AND. regularize) THEN
                DO r=1,nweights_tri
                   intaux3(:,r,n) = intVolGradGSelf(qpm_bnd(:,r), n, mesh, k, ga, prd, qd_tetra, xi)
                END DO
             ELSE
                DO r=1,nweights_tri
                   intaux3(:,r,n) = intVolGradG(qpm_bnd(:,r), n, mesh, k, ga, prd, qd_tetra, xi)
                END DO
             END IF
          END DO

          DO p=1,4
             CALL vsolid_rwg(qpm_bnd(:,:),m,p,mesh,fmv_bnd(:,:,p))
          END DO

          DO n=1,mesh%nsolids
             
             DO p=1,4
                DO q=1,4
                   
                   int1 = 0.0_dp
                   DO r=1,nweights_tri
                      int1 = int1 + qd_tri%weights(r)*dotr(fmv_bnd(:,r,p),nor)*intaux3(q,r,n)
                   END DO
                   
                   index1 = mesh%solids(m)%solid_face_indices(p)
                   index2 = mesh%solids(n)%solid_face_indices(q)
                   
                   A(index1,index2) = A(index1,index2) - Am*int1
                END DO
             END DO
          END DO

       END DO

    END DO

  END SUBROUTINE vieGreenMatrix

  ! Integrate 3x4 complex matrix valued function of three real variables
  ! over [-1,1]^3 by Gauss-Legendre rule.
  FUNCTION quadGL3D_34(f) RESULT(res)
    INTERFACE
       FUNCTION f(x,y,z) RESULT(fres)
         DOUBLE PRECISION, INTENT(IN) :: x,y,z
         DOUBLE COMPLEX, DIMENSION(3,4) :: fres
       END FUNCTION f
    END INTERFACE
    INTEGER :: i, j, k, nw
    COMPLEX, DIMENSION(3,4) :: res

    nw = SIZE(GL1Dw)

    res(:,:) = 0.0_dp

    DO i=1,nw
       DO j=1,nw
          DO k=1,nw
             res = res + GL1Dw(i)*GL1Dw(j)*GL1Dw(k)*f(GL1Dn(i),GL1Dn(j),GL1Dn(k))
          END DO
       END DO
    END DO
  END FUNCTION quadGL3D_34

  ! Performs the same as intVolG, but assumes that r is within the
  ! given solid. Duffy type regularization is applied.
  FUNCTION intVolGSelf(r, solidind, mesh, k, ga, prd, qd, xi) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: solidind
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    COMPLEX (KIND=dp), DIMENSION(3,4) :: res
    REAL (KIND=dp), DIMENSION(3,4) :: nodes
    INTEGER :: n, m
    REAL (KIND=dp) :: tetraJac

    ! This node is where singularity occurs.
    nodes(:,4) = r

    res(:,:) = 0.0_dp

    ! Split tetrahedron at r into four sub-tetrahedra.
    DO n=1,4
       ! Form the first three solid nodes.
       DO m=1,3
          nodes(:,m) = mesh%nodes(mesh%solids(solidind)%node_indices(subsolidind(m,n)))%p
       END DO
       
       ! Jacobian from canonical tetrahedron to current world tetrahedron.
       tetraJac = ABS(dotr(crossr(nodes(:,1)-nodes(:,4), nodes(:,2)-nodes(:,4)), nodes(:,3)-nodes(:,4)))

       IF(tetraJac/=0.0_dp) THEN
          res = res + quadGL3D_34(integrand)
       END IF
    END DO

    CONTAINS
      FUNCTION integrand(xi1, xi2, xi3) RESULT(integ)
        DOUBLE PRECISION, INTENT(IN) :: xi1, xi2, xi3
        DOUBLE COMPLEX, DIMENSION(3,4) :: integ

        REAL (KIND=dp) :: rho, phi, theta, rhoLim, sp, cp, st, ct, u, v, w, jac
        REAL (KIND=dp), DIMENSION(3) :: fn, pos
        COMPLEX (KIND=dp) :: g
        INTEGER :: t

        ! Compute spherical polar coords.
        phi = 0.25_dp*pi*(xi2+1.0_dp)
        theta = 0.25_dp*pi*(xi3+1.0_dp)

        sp = SIN(phi)
        cp = COS(phi)
        st = SIN(theta)
        ct = COS(theta)

        rhoLim = 1.0_dp/(cp*st + sp*st + ct)
        rho = 0.5_dp*(xi1+1.0_dp)*rhoLim

        ! Compute canonical tetrahedron coords.
        u = rho*cp*st
        v = rho*sp*st
        w = rho*ct

        ! Position in currently constructed tetrahedron in the lab frame (x,y,z).
        pos = MATMUL(nodes, (/u, v, w, 1.0_dp-u-v-w/))

        ! Full jacobian from (xi1,xi2,xi3) -> (x,y,z).
        jac = 0.5_dp*rhoLim*(rho**2)*st*((0.25_dp*pi)**2)*tetraJac

        g = gf(r, pos, k)

        DO t=1,4
           fn = solid_rwg(pos,solidind,t,mesh)

           integ(:,t) = jac*g*MATMUL(xi(pos,solidind), fn)
        END DO

      END FUNCTION integrand

  END FUNCTION intVolGSelf

  ! Computes int_Vn G(r,r')*xi(r')*fn(r')dV'
  ! where xi is the contrast tensor function (I-epsr^-1)
  ! and fn is SWG basis function. Result is returned for all four
  ! fn supported by the solid.
  FUNCTION intVolG(r, solidind, mesh, k, ga, prd, qd, xi) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: solidind
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    COMPLEX (KIND=dp), DIMENSION(3,4) :: res

    INTEGER :: t, nweights, faceind
    REAL (KIND=dp) :: Vn
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn
    COMPLEX (KIND=dp), DIMENSION(qd%num_nodes) :: gv
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: fv
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes) :: aux

    nweights = qd%num_nodes

    qpn = quad_tetra_points(qd, solidind, mesh)
    Vn = mesh%solids(solidind)%volume

    CALL vGf(r, qpn, k, nweights, gv)

    gv = gv*qd%weights*Vn

    DO faceind=1,4
       CALL vsolid_rwg(qpn(:,:),solidind,faceind,mesh,fv)

       ! Product xi*fn at quad points.
       DO t=1,nweights
          aux(:,t) = MATMUL(xi(qpn(:,t),solidind), fv(:,t))
       END DO

       res(:,faceind) = MATMUL(aux, gv)
    END DO

  END FUNCTION intVolG

  ! Computes int_Vn grad'G(r,r').(xi(r')*fn(r'))dV'
  ! where xi is the contrast tensor function (I-epsr^-1)
  ! and fn is SWG basis function. Result is returned for all four
  ! fn supported by the solid.
  FUNCTION intVolGradG(r, solidind, mesh, k, ga, prd, qd, xi) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: solidind
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    COMPLEX (KIND=dp), DIMENSION(4) :: res

    INTEGER :: t, nweights, faceind
    REAL (KIND=dp) :: Vn
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes) :: ggv
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: fv
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes) :: aux

    nweights = qd%num_nodes

    qpn = quad_tetra_points(qd, solidind, mesh)
    Vn = mesh%solids(solidind)%volume

    CALL vgradGf(r, qpn, k, nweights, ggv)

    res(:) = 0.0_dp

    DO faceind=1,4
       CALL vsolid_rwg(qpn(:,:),solidind,faceind,mesh,fv)

       ! Product xi*fn at quad points.
       DO t=1,nweights
          aux(:,t) = MATMUL(xi(qpn(:,t),solidind), fv(:,t))
          res(faceind) = res(faceind) + dotc(ggv(:,t), aux(:,t))*qd%weights(t)*Vn
       END DO
    END DO

  END FUNCTION intVolGradG

  ! Integrate 1x4 complex matrix valued function of three real variables
  ! over [-1,1]^3 by Gauss-Legendre rule.
  FUNCTION quadGL3D_14(f) RESULT(res)
    INTERFACE
       FUNCTION f(x,y,z) RESULT(fres)
         DOUBLE PRECISION, INTENT(IN) :: x,y,z
         DOUBLE COMPLEX, DIMENSION(4) :: fres
       END FUNCTION f
    END INTERFACE
    INTEGER :: i, j, k, nw
    COMPLEX, DIMENSION(4) :: res

    nw = SIZE(GL1Dw)

    res(:) = 0.0_dp

    DO i=1,nw
       DO j=1,nw
          DO k=1,nw
             res = res + GL1Dw(i)*GL1Dw(j)*GL1Dw(k)*f(GL1Dn(i),GL1Dn(j),GL1Dn(k))
          END DO
       END DO
    END DO
  END FUNCTION quadGL3D_14

  ! Performs the same as intVolGradG, but assumes that r is within the
  ! given solid. Duffy type regularization is applied.
  FUNCTION intVolGradGSelf(r, solidind, mesh, k, ga, prd, qd, xi) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: solidind
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd

    INTERFACE
       FUNCTION xi(pos, s) RESULT(xires)
         DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
         INTEGER, INTENT(IN) :: s
         COMPLEX, DIMENSION(3,3) :: xires
       END FUNCTION xi
    END INTERFACE

    COMPLEX (KIND=dp), DIMENSION(4) :: res
    REAL (KIND=dp), DIMENSION(3,4) :: nodes
    INTEGER :: n, m
    REAL (KIND=dp) :: tetraJac

    ! This node is where singularity occurs.
    nodes(:,4) = r

    res(:) = 0.0_dp

    ! Split tetrahedron at r into four sub-tetrahedra.
    DO n=1,4
       ! Form the first three solid nodes.
       DO m=1,3
          nodes(:,m) = mesh%nodes(mesh%solids(solidind)%node_indices(subsolidind(m,n)))%p
       END DO
       
       ! Jacobian from canonical tetrahedron to current world tetrahedron.
       tetraJac = ABS(dotr(crossr(nodes(:,1)-nodes(:,4), nodes(:,2)-nodes(:,4)), nodes(:,3)-nodes(:,4)))

       IF(tetraJac/=0.0_dp) THEN       
          res = res + quadGL3D_14(integrand)
       END IF
    END DO

    CONTAINS
      FUNCTION integrand(xi1, xi2, xi3) RESULT(integ)
        DOUBLE PRECISION, INTENT(IN) :: xi1, xi2, xi3
        DOUBLE COMPLEX, DIMENSION(4) :: integ

        REAL (KIND=dp) :: rho, phi, theta, rhoLim, sp, cp, st, ct, u, v, w, jac
        REAL (KIND=dp), DIMENSION(3) :: fn, pos
        COMPLEX (KIND=dp), DIMENSION(3) :: gg
        INTEGER :: t

        ! Compute spherical polar coords.
        phi = 0.25_dp*pi*(xi2+1.0_dp)
        theta = 0.25_dp*pi*(xi3+1.0_dp)

        sp = SIN(phi)
        cp = COS(phi)
        st = SIN(theta)
        ct = COS(theta)

        rhoLim = 1.0_dp/(cp*st + sp*st + ct)
        rho = 0.5_dp*(xi1+1.0_dp)*rhoLim

        ! Compute canonical tetrahedron coords.
        u = rho*cp*st
        v = rho*sp*st
        w = rho*ct

        ! Position in currently constructed tetrahedron in the lab frame (x,y,z).
        pos = MATMUL(nodes, (/u, v, w, 1.0_dp-u-v-w/))

        ! Full jacobian from (xi1,xi2,xi3) -> (x,y,z).
        jac = 0.5_dp*rhoLim*(rho**2)*st*((0.25_dp*pi)**2)*tetraJac

        gg = gradGf(r, pos, k)

        DO t=1,4
           fn = solid_rwg(pos,solidind,t,mesh)

           integ(t) = jac*dotc(gg, MATMUL(xi(pos,solidind), fn))
        END DO

      END FUNCTION integrand

    END FUNCTION intVolGradGSelf
END MODULE vie
