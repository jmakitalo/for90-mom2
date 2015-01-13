MODULE pec
  USE sysmat
  USE source
  USE bc

  IMPLICIT NONE

CONTAINS
  SUBROUTINE pec_modes(meshfile, scale, omega)
    CHARACTER (LEN=*), INTENT(IN) :: meshfile
    REAL (KIND=dp), INTENT(IN) :: scale, omega
    TYPE(quad_data) :: qd_tri, qd_tetra
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: D, H, Z, R, X
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval, eigvecexp
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec, curdens
    TYPE(mesh_container) :: mesh
    COMPLEX (KIND=dp) :: ri
    TYPE(srcdata), ALLOCATABLE :: src
    TYPE(prdnfo), POINTER :: prd
    TYPE(group_action), DIMENSION(:), ALLOCATABLE :: ga
    TYPE(domain), DIMENSION(2) :: domains
    TYPE(medium), DIMENSION(2) :: media
    INTEGER :: n, m, dim
    CHARACTER (LEN=256) :: numstr, oname
    INTEGER, DIMENSION(6) :: minind
    INTEGER, DIMENSION(:), ALLOCATABLE :: id
    
    ri = 1
    
    prd => NULL()
    
    CALL group_id(ga)

    WRITE(oname, '(A,A)') meshfile, '.msh'

    mesh = load_mesh(oname)
    
    CALL build_mesh(mesh, scale)
    
    ALLOCATE(D(1:mesh%nedges, 1:mesh%nedges))
    ALLOCATE(H(1:mesh%nedges, 1:mesh%nedges))
    ALLOCATE(id(1:mesh%nedges))
     
    ! Quadrature rules.
    qd_tri = tri_quad_data('tri_gl4')
    qd_tetra = tetra_quad_data('tetra_gl1')
     
    CALL computeD(omega, ri, mesh, ga(1), prd, qd_tri, D)
    CALL computeH(omega, ri, mesh, ga(1), prd, qd_tri, H)

    ! Find boundary edges.
    dim = 0
    DO n=1,mesh%nedges
       IF(mesh%edges(n)%face_indices(1)/=-1 .AND. mesh%edges(n)%face_indices(2)/=-1) THEN
          dim = dim + 1
          id(dim) = n
       END IF
    END DO

    WRITE(*,*) 'dim ', dim

    ALLOCATE(Z(1:dim, 1:dim))
    ALLOCATE(R(1:dim, 1:dim))
    ALLOCATE(X(1:dim, 1:dim))
    ALLOCATE(eigval(1:dim), eigvec(1:dim, 1:dim), eigvecexp(1:mesh%nedges))
    ALLOCATE(curdens(3,1:mesh%nfaces))

    Z = D(id(1:dim),id(1:dim)) + H(id(1:dim),id(1:dim))

    DEALLOCATE(D, H)

    R = REAL(Z, KIND=dp)
    X = AIMAG(Z)

    WRITE(*,*) 'Computing eigenvalues'

    CALL matrix_eigenvalues_gen(X, R, eigval, eigvec)

    !CALL rwg_moments(mesh, qd_tri, D)
    !CALL matrix_inverse(D(id(1:dim),id(1:dim)), R)
    !Z = MATMUL(R, Z)
    !CALL matrix_eigenvalues(Z, eigval, eigvec)

    minind = find_smallest(eigval, dim, 6)

    WRITE(*,*) eigval(minind(:))

    WRITE(*,*) 'Plotting modes'

    DO n=1,6
       eigvecexp(:) = 0.0_dp
       eigvecexp(id(1:dim)) = eigvec(1:dim, minind(n))

       DO m=1,mesh%nfaces
          curdens(:,m) = rwg_exp(mesh%faces(m)%cp, mesh, m, eigvecexp)
       END DO

       WRITE(numstr, '(I0)') n
       oname = TRIM(meshfile) // '-mode' // TRIM(ADJUSTL(numstr)) // '.msh'

       CALL save_vector_fields_msh(oname, mesh, curdens, curdens, scale)
    END DO
    
    DEALLOCATE(ga, Z, R, X, eigval, eigvec, curdens, id, eigvecexp)
    
    CALL delete_mesh(mesh)
    CALL delete_quad_data(qd_tri)
    CALL delete_quad_data(qd_tetra)
  END SUBROUTINE pec_modes
END MODULE pec
