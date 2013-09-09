! MODULE: symmetry
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implements the type group_action and routines for
! manipulating arrays of this type.
MODULE symmetry
  USE linalg

  IMPLICIT NONE

  REAL (KIND=dp), DIMENSION(3,3), PARAMETER :: id33r = (/1.0_dp,0.0_dp,0.0_dp,&
       0.0_dp,1.0_dp,0.0_dp,&
       0.0_dp,0.0_dp,1.0_dp/)

  INTEGER, PARAMETER :: gid_identity = 1,&
       gid_mxp = 2,&
       gid_myp = 3,&
       gid_mzp = 4,&
       gid_rz = 5

  TYPE group_action
     ! Action for electric field. Action for H-field is ef*detj.
     ! There can be multiple sub-problems with different field actions.
     COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: ef

     ! Jacobian the group action for points in R^3.
     REAL (KIND=dp), DIMENSION(3,3) :: j

     ! Determinant of the jacobian.
     REAL (KIND=dp) :: detj

     ! Identifier of a special action.
     INTEGER :: id
  END type group_action

CONTAINS
  SUBROUTINE group_id(ga)
    TYPE(group_action), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ga

    ALLOCATE(ga(1))
    ALLOCATE(ga(1)%ef(1))

    ga(1)%j = id33r
    ga(1)%detj = 1.0_dp
    ga(1)%ef(1) = 1.0_dp
    ga(1)%id = gid_identity
  END SUBROUTINE group_id

  SUBROUTINE group_mp(np, ga)
    TYPE(group_action), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ga
    INTEGER, INTENT(IN) :: np
    INTEGER :: n, na

    na = 2

    ALLOCATE(ga(na))
    DO n=1,na
       ALLOCATE(ga(n)%ef(na))
    END DO

    ga(1)%ef = (/1.0_dp, 1.0_dp/)
    ga(1)%j = id33r
    ga(1)%detj = 1.0_dp
    ga(1)%id = gid_identity

    ga(2)%ef = (/1.0_dp, -1.0_dp/)
    ga(2)%j = id33r
    ga(2)%j(np,np) = -1.0_dp
    ga(2)%detj = -1.0_dp

    IF(np==1) THEN
       ga(2)%id = gid_mxp
    ELSE IF(np==2) THEN
       ga(2)%id = gid_myp
    ELSE IF(np==3) THEN
       ga(2)%id = gid_mzp
    END IF
  END SUBROUTINE group_mp

  ! Group of n rotation with respect to the z-axis.
  SUBROUTINE group_rz(n, ga)
    INTEGER, INTENT(IN) :: n
    TYPE(group_action), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ga
    INTEGER :: k, l, na
    REAL (KIND=dp) :: angle
    COMPLEX (KIND=dp) :: phase

    na = n

    ALLOCATE(ga(na))
    DO k=1,na
       ALLOCATE(ga(k)%ef(na))
    END DO

    ! Set point actions.
    DO k=1,n
       angle = 2.0_dp*pi*(k-1)/n

       ga(k)%j = matrix_rz(angle)
       ga(k)%detj = 1.0_dp
       ga(k)%id = 0
    END DO

    ga(1)%id = gid_identity
    ga(n)%id = gid_rz

    ! Declare field actions to unity.
    DO k=1,n
       ga(k)%ef(:) = 1.0_dp
    END DO

    ! Generator of the 1-dimensional complex group.
    angle = 2.0_dp*pi/n
    phase = EXP(-(0,1)*angle)

    ! Field actions of all representations.
    DO k=2,n
       DO l=2,n
          ga(k)%ef(l) = phase**(k+l-2)
          IF(MOD(k+l-2,n)==0) THEN
             ga(k)%ef(l) = phase
          END IF
       END DO
    END DO
  END SUBROUTINE group_rz

  ! Given groups ga1 and ga2, returns the group ga,
  ! whose elements are products of members in the two groups.
  SUBROUTINE product_group(ga, ga2)
    TYPE(group_action), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ga
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga2
    TYPE(group_action), DIMENSION(:), ALLOCATABLE :: ga1
    INTEGER :: n, m, k, na1, na2, na, n2, m2, k2

    na1 = SIZE(ga)
    na2 = SIZE(ga2)

    na = na1*na2

    ! Copy ga to ga1.
    ALLOCATE(ga1(na1))
    DO n=1,na1
       ALLOCATE(ga1(n)%ef(na1))
       ga1(n)%ef(:) = ga(n)%ef(:)
       ga1(n)%j(:,:) = ga(n)%j(:,:)
       ga1(n)%detj = ga(n)%detj
       ga1(n)%id = ga(n)%id
    END DO

    ! Destroy old ga.
    DO n=1,na1
       DEALLOCATE(ga(n)%ef)
    END DO
    DEALLOCATE(ga)

    ! Reallocate ga to hold the product group.
    ALLOCATE(ga(na))
    DO n=1,na
       ALLOCATE(ga(n)%ef(na))
    END DO

    DO n=1,na1
       DO m=1,na2
          k = (n-1)*na2 + m

          ga(k)%j = MATMUL(ga1(n)%j, ga2(m)%j)
          ga(k)%detj = ga1(n)%detj*ga2(m)%detj
          ga(k)%id = 0

          IF(ga1(n)%id==gid_identity .AND. ga2(m)%id==gid_identity) THEN
             ga(k)%id = gid_identity
          ELSE IF(ga1(n)%id/=0 .AND. ga1(n)%id/=gid_identity .AND. ga2(m)%id==gid_identity) THEN
             ga(k)%id = ga1(n)%id
          ELSE IF(ga2(m)%id/=0 .AND. ga2(m)%id/=gid_identity .AND. ga1(n)%id==gid_identity) THEN
             ga(k)%id = ga2(m)%id
          END IF

          DO n2=1,na1
             DO m2=1,na2
                k2 = (n2-1)*na2 + m2

                ga(k)%ef(k2) = ga1(n)%ef(n2)*ga2(m)%ef(m2)
             END DO
          END DO
       END DO
    END DO

    DEALLOCATE(ga1)
  END SUBROUTINE product_group

END MODULE symmetry