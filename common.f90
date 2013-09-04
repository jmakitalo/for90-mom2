! MODULE: common
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implements various types for the description of scattering
! problems, including media and domains.
MODULE common
  USE source
  USE greenprd
  USE nlsurf

  IMPLICIT NONE

  TYPE medium_nlb
     COMPLEX (KIND=dp) :: delta_prime
     COMPLEX (KIND=dp) :: gamma
  END type medium_nlb

  TYPE medium_prop
     COMPLEX (KIND=dp) :: ri
     COMPLEX (KIND=dp) :: shri
     TYPE(medium_nls) :: nls
     TYPE(medium_nlb) :: nlb
  END type medium_prop

  TYPE medium
     INTEGER :: type

     ! Range of wavelengths.
     TYPE(medium_prop), DIMENSION(:), ALLOCATABLE :: prop
  END type medium

  ! Medium types.
  INTEGER, PARAMETER :: mtype_linear = 1,&
       mtype_nls = 2,&
       mtype_nlb = 3

  ! Solution data.
  TYPE solution
     ! Solution vector. First dimension denotes basis coefficients of (J,M).
     ! Second dimension denotes sub-problems related to group representations.
     COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: x, nlx

     ! Eigenvectors of a spectral problem. First dimension denotes basis coefficients.
     ! Second dimension denotes eigenvalue index.
     COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec

     ! Eigenvalues of a spectral problem.
     COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval

     ! Wavelength.
     REAL (KIND=dp) :: wl
  END TYPE solution

  ! Data for computational domain or a sub-domain.
  TYPE domain
     ! Index of precomputed Green's function data. Index -1 denotes non-periodic GF.
     INTEGER :: gf_index

     ! Index for medium associated with the domain.
     INTEGER :: medium_index

     ! Mesh of the domain boundary. May contain also a volume mesh.
     TYPE(mesh_container) :: mesh
  END TYPE domain

  ! Data for a computation batch, which way involve a set of wavelengths.
  TYPE batch
     ! Name of the problem (used for output file naming).
     CHARACTER (LEN=256) :: name

     ! Name of the mesh file.
     CHARACTER (LEN=256) :: mesh_file

     ! Mesh data (may contain multiple domains).
     TYPE(mesh_container) :: mesh

     ! Scale factor for the mesh.
     REAL (KIND=dp) :: scale

     ! nwl: Number of wavelengths.
     INTEGER :: nwl

     ! Source data.
     TYPE(srcdata) :: src

     ! Solution data for each wavelength.
     TYPE(solution), DIMENSION(:), ALLOCATABLE :: sols

     ! Periodic Green's function data.
     TYPE(prdnfo), DIMENSION(:), POINTER :: prd

     ! Array of group actions.
     TYPE(group_action), DIMENSION(:), ALLOCATABLE :: ga

     ! Domains of the problem.
     TYPE(domain), DIMENSION(:), ALLOCATABLE :: domains

     ! Media.
     TYPE(medium), DIMENSION(:), ALLOCATABLE :: media
  END TYPE batch

CONTAINS
  SUBROUTINE batch_defaults(b)
    TYPE(batch), INTENT(INOUT) :: b

    b%name = 'unnamed'
    b%mesh_file = 'default.msh'
    b%scale = 1d-9
    b%nwl = 0
    b%src%type = 0
  END SUBROUTINE batch_defaults

  SUBROUTINE delete_batch(b)
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: i

    CALL delete_mesh(b%mesh)

    IF(ALLOCATED(b%sols)) THEN
       DO i=1,SIZE(b%sols)
          CALL delete_solution(b%sols(i))
       END DO

       DEALLOCATE(b%sols)
    END IF

    IF(ALLOCATED(b%ga)) THEN
       DO i=1,SIZE(b%ga)
          DEALLOCATE(b%ga(i)%ef)
       END DO

       DEALLOCATE(b%ga)
    END IF

    IF(ALLOCATED(b%domains)) THEN
       DO i=1,SIZE(b%domains)
          CALL delete_domain(b%domains(i))
       END DO

       DEALLOCATE(b%domains)
    END IF

    IF(ALLOCATED(b%media)) THEN
       DEALLOCATE(b%media)
    END IF

    IF(ASSOCIATED(b%prd)) THEN
       DO i=1,SIZE(b%prd)
          CALL clear_prd(b%prd(i))
       END DO

       DEALLOCATE(b%prd)
    END IF

    IF(ALLOCATED(b%media)) THEN
       DO i=1,SIZE(b%media)
          DEALLOCATE(b%media(i)%prop)
       END DO

       DEALLOCATE(b%media)
    END IF
  END SUBROUTINE delete_batch

  SUBROUTINE delete_solution(s)
    TYPE(solution), INTENT(INOUT) :: s

    IF(ALLOCATED(s%x)) THEN
       DEALLOCATE(s%x)
    END IF

    IF(ALLOCATED(s%nlx)) THEN
       DEALLOCATE(s%nlx)
    END IF

    IF(ALLOCATED(s%eigvec)) THEN
       DEALLOCATE(s%eigvec)
    END IF

    IF(ALLOCATED(s%eigval)) THEN
       DEALLOCATE(s%eigval)
    END IF
  END SUBROUTINE delete_solution

  SUBROUTINE delete_domain(d)
    TYPE(domain), INTENT(INOUT) :: d

    CALL delete_mesh(d%mesh)
  END SUBROUTINE delete_domain
END MODULE common
