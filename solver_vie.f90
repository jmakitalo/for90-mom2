MODULE solver_vie
  USE time
  USE vie
  USE common

  IMPLICIT NONE

CONTAINS
  SUBROUTINE solve_batch_vie(b)
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: n, m, l, nbasis, nga, nsrc
    REAL (KIND=dp) :: wl, omega
    COMPLEX (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: A

    COMPLEX (KIND=dp) :: phdx, phdy, ri, k
    TYPE(prdnfo), POINTER :: prd
    TYPE(medium_prop) :: mprop

    WRITE(*,*) '--- Begin wavelength batch ---'

    ! Check that necessary input data exists.
    IF(ALLOCATED(b%sols)==.FALSE.) THEN
       WRITE(*,*) 'Setup wavelengths prior to solving!'
       RETURN
    END IF

    IF(ALLOCATED(b%mesh%nodes)==.FALSE.) THEN
       WRITE(*,*) 'Load mesh prior to solving!'
       RETURN
    END IF

    IF(b%mesh%nsolid_faces==0) THEN
       CALL build_solid_faces(b%mesh)
       CALL compute_basis_data(b%mesh)
    END IF

    IF(ALLOCATED(b%domains)==.FALSE.) THEN
       WRITE(*,*) 'Set up domains prior to solving!'
       RETURN
    END IF

    IF(ALLOCATED(b%media)==.FALSE.) THEN
       WRITE(*,*) 'Set up media prior to solving!'
       RETURN
    END IF

    IF(ALLOCATED(b%src)==.FALSE.) THEN
       WRITE(*,*) 'Set up source prior to solving!'
       RETURN
    END IF

    nbasis = b%mesh%nsolid_faces
    nga = SIZE(b%ga)
    nsrc = SIZE(b%src)

    ALLOCATE(A(nbasis,nbasis,nga))

    WRITE(*,*) 'Name: ', TRIM(b%name)

    WRITE(*,*) 'Triangle quadrature: ', TRIM(b%qd_tri%description)
    WRITE(*,*) 'Tetrahedron quadrature: ', TRIM(b%qd_tetra%description)

    CALL print_source_info(b%src(1))

    prd => NULL()

    ! Go through all given wavelengths.
    DO n=1, b%nwl
       wl = b%sols(n)%wl
       omega = 2.0_dp*pi*c0/wl

       ri = b%media(b%domains(1)%medium_index)%prop(n)%ri
       k = ri*omega/c0

       ! Setup periodic GF for current wavelength.
       DO m=1,SIZE(b%domains)
          IF(b%domains(m)%gf_index/=-1) THEN
             b%prd(b%domains(m)%gf_index)%cwl = find_closest(wl,&
                  b%prd(b%domains(m)%gf_index)%coef(:)%wl)
          END IF
       END DO

       ! Set lattice cell phase shifts according to domain 1.
       IF(b%domains(1)%gf_index/=-1) THEN
          prd => b%prd(b%domains(1)%gf_index)
          phdx = EXP((0,1)*prd%dx*prd%coef(prd%cwl)%k0x)
          phdy = EXP((0,1)*prd%dy*prd%coef(prd%cwl)%k0y)
       END IF

       ! Print some information.
       WRITE(*,'(A,F6.1,A,I0,A,I0,A)') ' Wavelength: ', wl*1d9, ' nm (', n, ' of ', b%nwl, ')'

       DO m=1,SIZE(b%domains)
          ri = b%media(b%domains(m)%medium_index)%prop(n)%ri
          WRITE(*,'(A,I0,A,"(",F6.3,",",F6.3,")")') ' Refractive index of domain ', m, ': ', ri
       END DO

       WRITE(*,*) 'Building matrices'

       ! Compute the system matrix for all group representations.
       ! This is O(N^2), but in practice the most time consuming part.
       CALL timer_start()
       CALL cpu_timer_start()
       CALL vie_matrix(b%mesh, k, b%ga(1), prd, b%qd_tri, b%qd_tetra, xi_hom, A(:,:,1))
       WRITE(*,*) 'Wall-clock time:'
       WRITE(*,*) sec_to_str(timer_end())
       WRITE(*,*) 'CPU time:'
       WRITE(*,*) sec_to_str(cpu_timer_end())

       ! Allocate memory for source/solution vector.
       ALLOCATE(b%sols(n)%x(nbasis, nga, nsrc))

       ! Compute excitation source vectors for each representation.
       ! Excitation is assumed to be in domain 1.
       WRITE(*,*) 'Computing source vectors'
       CALL timer_start()

       ri = b%media(b%domains(1)%medium_index)%prop(n)%ri
       DO l=1,nsrc
          CALL vie_srcvec(b%mesh, omega, ri, b%ga(1), b%qd_tetra, b%src(l), b%sols(n)%x(:,1,l))
       END DO
       WRITE(*,*) sec_to_str(timer_end())

       ! Solve the linear system of equations for each representation,
       ! enforcing possible boundary conditions.
       ! This is O(N^3) but often in practice very fast.
       WRITE(*,*) 'Solving system'
       CALL timer_start()
       CALL cpu_timer_start()

       ! Solve the remaining linear system of dimension dim.
       CALL solve_multi_linsys(A(:,:,1), b%sols(n)%x(:,1,:))

       WRITE(*,*) 'Wall-clock time:'
       WRITE(*,*) sec_to_str(timer_end())
       WRITE(*,*) 'CPU time:'
       WRITE(*,*) sec_to_str(cpu_timer_end())

    END DO

    DEALLOCATE(A)

    WRITE(*,*) '--- End wavelength batch ---'

  CONTAINS
    FUNCTION xi_hom(pos, s) RESULT(xires)
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: pos
      INTEGER, INTENT(IN) :: s
      COMPLEX (KIND=dp) :: eps, diag
      COMPLEX, DIMENSION(3,3) :: xires

      eps = b%media(b%domains(2)%medium_index)%prop(n)%ri**2
      diag = 1.0_dp - 1.0_dp/eps

      xires(:,:) = 0.0_dp

      xires(1,1) = diag
      xires(2,2) = diag
      xires(3,3) = diag
    END FUNCTION xi_hom

  END SUBROUTINE solve_batch_vie
END MODULE solver_vie
