! MODULE: bc
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for detecting and enforcing boundary conditions for
! the surface current densities. Boundary conditions may appear
! from symmetry and periodicity (translational symmetry).
MODULE bc
  USE symmetry
  USE mesh

  IMPLICIT NONE

CONTAINS
  ! Determines possible boundary conditions for edge elements.
  ! mesh: Surface mesh of the problem.
  ! ga: An array of group actions that describe the symmetries.
  ! phdx: Phase shift over periodic unit cell in x-direction.
  ! phdy: Phase shift over periodic unit cell in y-direction.
  ! nr: Group representation index.
  ! id: Identifier of boundary conditions for J and M.
  !  id==-1 means that corresponding row and column of the matrix
  !  must be removed.
  !  id==0 means that no boundary condition is obtained.
  !  id>0 means that column of index id is multiplied by given phase
  !  and added to corresponding column. Then row and column of index
  !  id must be removed.
  ! phase: A phase factor of a boundary condition.
  SUBROUTINE edge_bc(mesh, ga, phdx, phdy, nr, id, phase)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nr
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), INTENT(IN) :: phdx, phdy

    INTEGER, DIMENSION(:), INTENT(INOUT) :: id
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: phase

    INTEGER :: k, n, nedges, bnd, couple_coef
    COMPLEX (KIND=dp) :: gae, gah
    REAL (KIND=dp) :: detj
    REAL (KIND=dp), DIMENSION(3,3) :: j

    nedges = mesh%nedges

    id(:) = 0
    phase(:) = 1.0_dp

    DO k=1,nedges
       
       IF(mesh%edges(k)%couple_index>0) THEN
          couple_coef = mesh%edges(k)%couple_index
       END IF
       
       bnd = mesh%edges(k)%bnd
       
       ! Check if pseudo-periodic symmetry leads to boundary conditions.
       
       ! Remove coefficients related to "secondary" boundary edges.
       IF(bnd==mesh_bnd_prdx2 .OR. bnd==mesh_bnd_prdy2) THEN
          id(k) = -1
          id(k + nedges) = -1
       END IF
       
       ! Enforce Bloch phase-shift to boundary edges.
       IF(bnd==mesh_bnd_prdx1) THEN
          id(k) = couple_coef
          phase(k) = phdx
          
          id(k + nedges) = couple_coef + nedges
          phase(k + nedges) = phdx
       END IF
          
       IF(bnd==mesh_bnd_prdy1) THEN
          id(k) = couple_coef
          phase(k) = phdy
          
          id(k + nedges) = couple_coef + nedges
          phase(k + nedges) = phdy
       END IF
          
       ! Check if any group actions lead to boundary conditions.
       DO n=1,SIZE(ga)
          j = ga(n)%j
          detj = ga(n)%detj
          
          gae = ga(n)%ef(nr)
          gah = ga(n)%ef(nr)*detj
          
          ! Check for BC for plane mirror symmetries.
          IF((ga(n)%id==gid_mxp .AND. bnd==mesh_bnd_xplane) .OR.&
               (ga(n)%id==gid_myp .AND. bnd==mesh_bnd_yplane) .OR.&
               (ga(n)%id==gid_mzp .AND. bnd==mesh_bnd_zplane)) THEN
             
             IF(REAL(-gae)<0.99_dp) THEN
                id(k) = -1
             END IF
             
             IF(REAL(-gah)<0.99_dp) THEN
                id(k + nedges) = -1
             END IF
          END IF
          
          ! Check for BC for rotation by the z-axis.
          IF(ga(n)%id==gid_rz) THEN
             IF(bnd==mesh_bnd_rz1) THEN
                id(k) = couple_coef
                phase(k) = gah
                
                id(k + nedges) = couple_coef + nedges
                phase(k + nedges) = gae
             ELSE IF(bnd==mesh_bnd_rz2) THEN
                id(k) = -1
                id(k + nedges) = -1
             END IF
          END IF
          
       END DO
    END DO
  END SUBROUTINE edge_bc

  ! Returns indices to true elements in given logical array.
  ! Is is assumed that n == SIZE(loc) == COUNT(array).
  FUNCTION array_locations(array, n) RESULT(loc)
    LOGICAL, DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: n
    INTEGER, DIMENSION(n) :: loc
    INTEGER :: m, k

    k = 0

    DO m=1,SIZE(array)
       IF(array(m)) THEN
          k = k + 1
          loc(k) = m
       END IF
    END DO
  END FUNCTION array_locations

  SUBROUTINE resolve_system_dependencies(A, b, id, phase)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: A
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: b
    INTEGER, DIMENSION(:), INTENT(IN) :: id
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: phase

    INTEGER :: n

    ! If the symmetry implies linear dependence between elements, resolve these
    ! dependencies by adding a column of the matrix, multiplied by appropriate phase,
    ! to another column. Note that id(n)>0 => id(id(n)) == -1.
    DO n=1,SIZE(A,2)
       IF(id(n)>0) THEN
          A(:,n) = A(:,n) + A(:,id(n))*phase(n)
       END IF
    END DO

    DO n=1,SIZE(A,1)
       IF(id(n)>0) THEN
          A(n,:) = A(n,:) + A(id(n),:)*CONJG(phase(n))

          b(n,:) = b(n,:) + b(id(n),:)*CONJG(phase(n))
       END IF
    END DO
  END SUBROUTINE resolve_system_dependencies

  ! System matrix A and source vector b are shrinked by dimension according
  ! to given BC identifiers.
  ! A: The system matrix.
  ! b: The source vector.
  ! dim: The reduced dimension. Valid ranges of A and b are A(1:dim,1:dim) and b(1:dim).
  ! id: Boundary condition identifiers. SIZE(id) == nbasis*2.
  SUBROUTINE reduce_system(A, b, dim, id)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: A
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: b
    INTEGER, INTENT(OUT) :: dim
    INTEGER, DIMENSION(:), INTENT(IN) :: id

    INTEGER :: n
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind

    ! The reduced dimension is the number of non-negative coefficients.
    dim = COUNT(id/=-1)

    ! Allocate array to store indices to non-zero coefficients.
    ALLOCATE(ind(dim))

    ! Obtain indices to non-zero coefficients.
    ind = array_locations(id/=-1, dim)

    ! Reduce the system matrix and the source vector.
    A(1:dim,1:dim) = A(ind, ind)
    b(1:dim,:) = b(ind,:)

    DEALLOCATE(ind)
  END SUBROUTINE reduce_system

  ! Given a reduced form of the solution vector b and the id and phase data given by
  ! subroutine reduce_system, the solution is expanded into a form with dimension nbasis*2.
  ! This solution may then be used directly by post processing routines, which assume this
  ! full-dimensional form.
  SUBROUTINE expand_solution(dim, id, phase, b)
    INTEGER, INTENT(IN) :: dim
    INTEGER, DIMENSION(:), INTENT(IN) :: id
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: phase
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: b

    COMPLEX (KIND=dp), DIMENSION(SIZE(b,1),SIZE(b,2)) :: btmp
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind
    INTEGER :: n

    ! Make a copy of the solution.
    btmp(:,:) = b(:,:)

    ! Make zero the default value in the full-dimensional vector.
    b(:,:) = 0.0_dp

    ! Allocate array to store indices to non-zero coefficients.
    ALLOCATE(ind(dim))

    ! Obtain indices to non-zero coefficients.
    ind = array_locations(id/=-1, dim)

    ! Add non-zero solution coefficients to proper places in the expanded form.
    b(ind,:) = btmp(1:dim,:)

    ! Assign solution coefficients to linearly dependent couples.
    DO n=1,SIZE(b,1)
       IF(id(n)>0) THEN
          b(id(n),:) = b(n,:)*phase(n)
       END IF
    END DO

    DEALLOCATE(ind)
  END SUBROUTINE expand_solution

  ! Returns true if the group action is admissible for mapping a boundary
  ! integral operator defined over the surface mesh, which is a symmetry unit cell.
  ! Otherwise returns false.
  FUNCTION admissible_ga(mesh, ga, exterior) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    TYPE(group_action), INTENT(IN) :: ga
    LOGICAL, INTENT(IN) :: exterior
    LOGICAL :: res

    res = .TRUE.

    IF(exterior) THEN
       RETURN
    ELSE IF(BTEST(ga%genbits, gid_mxp) .AND. (has_mesh_bnd(mesh, mesh_bnd_xplane) .EQV. .FALSE.)) THEN
       res = .FALSE.
    ELSE IF(BTEST(ga%genbits, gid_myp) .AND. (has_mesh_bnd(mesh, mesh_bnd_yplane) .EQV. .FALSE.)) THEN
       res = .FALSE.
    ELSE IF(BTEST(ga%genbits, gid_mzp) .AND. (has_mesh_bnd(mesh, mesh_bnd_zplane) .EQV. .FALSE.)) THEN
       res = .FALSE.
    ELSE IF(BTEST(ga%genbits, gid_rz) .AND. (has_mesh_bnd(mesh, mesh_bnd_rz1) .EQV. .FALSE.)) THEN
       res = .FALSE.
    END IF

  END FUNCTION admissible_ga
END MODULE bc
