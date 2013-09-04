! MODULE: linalg
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Simple low-level linear algebraic functions mainly for
! operations in R^3. Contains also some wrappers for
! LAPACK routines.
MODULE linalg
  USE constants

  IMPLICIT NONE

  INTRINSIC CONJG

CONTAINS
  FUNCTION dotc(v1, v2) RESULT(res)
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: v1, v2
    COMPLEX (KIND=dp) :: res

  !  res = SUM(CONJG(v1)*v2)
    res = SUM(v1*v2)
  END FUNCTION dotc

  PURE FUNCTION dotr(v1, v2) RESULT(res)
    REAL (KIND=dp), DIMENSION(:), INTENT(IN) :: v1, v2
    REAL (KIND=dp) :: res

    res = SUM(v1*v2)
  END FUNCTION dotr

  FUNCTION normc(v) RESULT(res)
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: v
    REAL (KIND=dp) :: res

   ! res = SQRT(REAL(dotc(v,v), KIND=dp))
    res = SQRT(REAL(dotc(CONJG(v),v), KIND=dp))
  END FUNCTION normc

  PURE FUNCTION normr(v) RESULT(res)
    REAL (KIND=dp), DIMENSION(:), INTENT(IN) :: v
    REAL (KIND=dp) :: res

    !res = SQRT(dotr(v,v))
    res = SQRT(SUM(v*v))
  END FUNCTION normr

  FUNCTION crossr(v1, v2) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: v1, v2
    REAL (KIND=dp), DIMENSION(3) :: res

    res(1) = v1(2)*v2(3) - v1(3)*v2(2)
    res(2) = v1(3)*v2(1) - v1(1)*v2(3)
    res(3) = v1(1)*v2(2) - v1(2)*v2(1)
  END FUNCTION crossr

  FUNCTION crossc(v1, v2) RESULT(res)
    COMPLEX (KIND=dp), DIMENSION(3), INTENT(IN) :: v1, v2
    COMPLEX (KIND=dp), DIMENSION(3) :: res

    res(1) = v1(2)*v2(3) - v1(3)*v2(2)
    res(2) = v1(3)*v2(1) - v1(1)*v2(3)
    res(3) = v1(1)*v2(2) - v1(2)*v2(1)
  END FUNCTION crossc

  SUBROUTINE vcrosscr(v1, v2, N, res)
    COMPLEX (KIND=dp), DIMENSION(3,N), INTENT(IN) :: v1
    REAL (KIND=dp), DIMENSION(3,N), INTENT(IN) :: v2
    INTEGER, INTENT(IN) :: N
    COMPLEX (KIND=dp), DIMENSION(3,N), INTENT(INOUT) :: res

    res(1,:) = v1(2,:)*v2(3,:) - v1(3,:)*v2(2,:)
    res(2,:) = v1(3,:)*v2(1,:) - v1(1,:)*v2(3,:)
    res(3,:) = v1(1,:)*v2(2,:) - v1(2,:)*v2(1,:)
  END SUBROUTINE vcrosscr

  FUNCTION invmat33r(m) RESULT(mi)
    REAL (KIND=dp), DIMENSION(3,3), INTENT(IN) :: m 
    REAL (KIND=dp), DIMENSION(3,3) :: mi
    REAL (KIND=dp) :: det

    det = dotr(m(:,1), crossr(m(:,2), m(:,3)))

    mi(1,:) = crossr(m(:,2), m(:,3))
    mi(2,:) = crossr(m(:,3), m(:,1))
    mi(3,:) = crossr(m(:,1), m(:,2))

    mi = mi/det
  END FUNCTION invmat33r

  SUBROUTINE matrix_inverse(A, Ainv)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: A
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: Ainv
    INTEGER :: dim, INFO, LWORK
    INTEGER, DIMENSION(1:SIZE(A,1)) :: IPIV
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: WORK

    dim = SIZE(A,1)

    Ainv = A

    CALL ZGETRF(dim, dim, Ainv, dim, IPIV, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Matrix factorization failed!'
       STOP
    END IF

    ALLOCATE(WORK(dim))

    CALL ZGETRI(dim, Ainv, dim, IPIV, WORK, -1, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Solving of system failed!'
       STOP
    END IF

    LWORK = WORK(1)

    DEALLOCATE(WORK)
    ALLOCATE(WORK(LWORK))

    CALL ZGETRI(dim, Ainv, dim, IPIV, WORK, LWORK, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Solving of system failed!'
       STOP
    END IF

    DEALLOCATE(WORK)

  END SUBROUTINE matrix_inverse

  SUBROUTINE solve_linsys(mat, src)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: mat
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: src

    INTEGER, DIMENSION(SIZE(src)) :: IPIV
    INTEGER :: dim, INFO

    dim = SIZE(src)

    CALL ZGETRF(dim, dim, mat, dim, IPIV, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Matrix factorization failed!'
       STOP
    END IF

    CALL ZGETRS('N', dim, 1, mat, dim, IPIV, src, dim, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Solving of system failed!'
       STOP
    END IF
  END SUBROUTINE solve_linsys

  SUBROUTINE matrix_eigenvalues(A, eigval, eigvec)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: A
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: eigval
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: eigvec
    COMPLEX (KIND=dp), DIMENSION(SIZE(eigvec,1),SIZE(eigvec,2)) :: VL
    INTEGER :: dim
    INTEGER :: INFO, LWORK
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: WORK
    REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: RWORK
   
    dim = SIZE(A,1)

    LWORK = 2*dim
    ALLOCATE(WORK(1:LWORK), RWORK(1:(2*dim)))

    ! Query workspace.
    CALL ZGEEV('N', 'V', dim, A, dim, eigval, VL, 1, eigvec, dim, WORK, -1, RWORK, INFO)
    
    IF(INFO/=0) THEN
       WRITE(*,*) 'Matrix eigenvalue solution was unsuccessful!'
       IF(INFO<0) THEN
          WRITE(*,*) 'Argument ', -INFO, ' was illegal.'
       END IF

       STOP
    END IF
    
    LWORK = WORK(1)
    DEALLOCATE(WORK)
    ALLOCATE(WORK(1:LWORK))
    
    ! Solve eigenvalue problem.
    CALL ZGEEV('N', 'V', dim, A, dim, eigval, VL, 1, eigvec, dim, WORK, LWORK, RWORK, INFO)
    
    IF(INFO/=0) THEN
       WRITE(*,*) 'Matrix eigenvalue solution was unsuccessful!'
       IF(INFO<0) THEN
          WRITE(*,*) 'Argument ', -INFO, ' was illegal.'
       END IF

       STOP
    END IF

    DEALLOCATE(WORK, RWORK)

  END SUBROUTINE matrix_eigenvalues

  SUBROUTINE matmul_blocklt(B, A, dim)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: B
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: A
    INTEGER, INTENT(IN) :: dim

    CALL zgemm('N', 'N', dim, dim, dim, 1.0_dp, B, dim, A(1:dim,1:dim),&
         dim, 0.0_dp, A(1:dim,1:dim), dim)

    CALL zgemm('N', 'N', dim, dim, dim, 1.0_dp, B, dim, A((dim+1):(2*dim),1:dim),&
         dim, 0.0_dp, A((dim+1):(2*dim),1:dim), dim)

    CALL zgemm('N', 'N', dim, dim, dim, 1.0_dp, B, dim, A(1:dim,(dim+1):(2*dim)),&
         dim, 0.0_dp, A(1:dim,(dim+1):(2*dim)), dim)

    CALL zgemm('N', 'N', dim, dim, dim, 1.0_dp, B, dim, A((dim+1):(2*dim),(dim+1):(2*dim)),&
         dim, 0.0_dp, A((dim+1):(2*dim),(dim+1):(2*dim)), dim)

    !A(1:dim,1:dim) = MATMUL(B, A(1:dim,1:dim))
    !A((dim+1):(2*dim),1:dim) = MATMUL(B, A((dim+1):(2*dim),1:dim))
    !A(1:dim,(dim+1):(2*dim)) = MATMUL(B, A(1:dim,(dim+1):(2*dim)))
    !A((dim+1):(2*dim),(dim+1):(2*dim)) = MATMUL(B, A((dim+1):(2*dim),(dim+1):(2*dim)))
  END SUBROUTINE matmul_blocklt

  SUBROUTINE matmul_blockrt(A, B, dim)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: A
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: B
    INTEGER, INTENT(IN) :: dim

    CALL zgemm('N', 'N', dim, dim, dim, 1.0_dp, A(1:dim,1:dim), dim, B,&
         dim, 0.0_dp, A(1:dim,1:dim), dim)

    CALL zgemm('N', 'N', dim, dim, dim, 1.0_dp, A((dim+1):(2*dim),1:dim), dim, B,&
         dim, 0.0_dp, A((dim+1):(2*dim),1:dim), dim)

    CALL zgemm('N', 'N', dim, dim, dim, 1.0_dp, A(1:dim,(dim+1):(2*dim)), dim, B,&
         dim, 0.0_dp, A(1:dim,(dim+1):(2*dim)), dim)

    CALL zgemm('N', 'N', dim, dim, dim, 1.0_dp, A((dim+1):(2*dim),(dim+1):(2*dim)), dim, B,&
         dim, 0.0_dp, A((dim+1):(2*dim),(dim+1):(2*dim)), dim)

    !A(1:dim,1:dim) = MATMUL(A(1:dim,1:dim), B)
    !A((dim+1):(2*dim),1:dim) = MATMUL(A((dim+1):(2*dim),1:dim), B)
    !A(1:dim,(dim+1):(2*dim)) = MATMUL(A(1:dim,(dim+1):(2*dim)), B)
    !A((dim+1):(2*dim),(dim+1):(2*dim)) = MATMUL(A((dim+1):(2*dim),(dim+1):(2*dim)), B)
  END SUBROUTINE matmul_blockrt

  ! Matrix describing rotation around z-axis by given angle
  ! in radians.
  FUNCTION matrix_rz(angle) RESULT(rm)
    REAL (KIND=dp), INTENT(IN) :: angle
    REAL (KIND=dp), DIMENSION(3,3) :: rm
    REAL (KIND=dp) :: st, ct

    st = SIN(angle)
    ct = COS(angle)

    rm(:,1) = (/ct, st, 0.0_dp/)
    rm(:,2) = (/-st, ct, 0.0_dp/)
    rm(:,3) = (/0.0_dp, 0.0_dp, 1.0_dp/)
  END FUNCTION matrix_rz
END MODULE linalg
