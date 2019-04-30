! MODULE: common
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implements various data aggregation types for the description of scattering
! problems, including media and domains. These form the basis for the high-level
! code that is used to define certain types of problems via the user interface.
MODULE common
  USE source
  USE greenprd
  USE nlsurf
  USE nlbulk

  IMPLICIT NONE

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
       mtype_nlb_nonlocal = 3,&
       mtype_nlb_dipole = 4

  ! Solution data.
  TYPE solution
     ! Solution vector. First dimension denotes basis coefficients of (J,M).
     ! Second dimension denotes sub-problems related to group representations.
     ! Third dimension denotes excitation source.
     COMPLEX (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: x, nlx

     ! Basis coefficients for jumps in M and J due to surface sources.
     ! 1st dim: coefficients of M and J jump expansions.
     ! 2nd dim: domain, 3rd dim: group representation, 4th dim: excitation source.
     COMPLEX (KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: src_coef

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
     TYPE(srcdata), DIMENSION(:), ALLOCATABLE :: src

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

     ! Quadrature data.
     TYPE(quad_data) :: qd_tri
     TYPE(quad_data) :: qd_tetra
  END TYPE batch

CONTAINS
  SUBROUTINE batch_defaults(b)
    TYPE(batch), INTENT(INOUT) :: b

    b%name = 'unnamed'
    b%mesh_file = 'default.msh'
    b%scale = 1d-9
    b%nwl = 0
    ALLOCATE(b%src(1))
    b%src(1)%type = 0

    b%qd_tri = tri_quad_data('tri_gl13')
    b%qd_tetra = tetra_quad_data('tetra_gl4')
  END SUBROUTINE batch_defaults

  SUBROUTINE delete_batch(b)
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: i

    CALL delete_quad_data(b%qd_tri)
    CALL delete_quad_data(b%qd_tetra)

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

    IF(ALLOCATED(b%src)) THEN
       DEALLOCATE(b%src)
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

    IF(ALLOCATED(s%src_coef)) THEN
       DEALLOCATE(s%src_coef)
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
