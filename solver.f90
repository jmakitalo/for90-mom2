! MODULE: solver
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Contains high-level routines for solving linear and nonlinear
! scattering problems, e.g., for a range of wavelengths.
MODULE solver
  USE sysmat
  USE source
  USE bc
  USE time

  IMPLICIT NONE

CONTAINS
  ! Solves linear and optionally second-order scattering problems
  ! for a range of wavelengths as specified in input structure b.
  SUBROUTINE solve_batch(b)
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: n, m, l, nbasis, dim, nf, nga, nind, nsrc
    REAL (KIND=dp) :: wl, omega
    COMPLEX (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: A

    COMPLEX (KIND=dp) :: phdx, phdy, ri
    TYPE(prdnfo), POINTER :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: src_coef
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: src_vec
    COMPLEX (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: src_vec2
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: epsp
    INTEGER, DIMENSION(b%mesh%nedges) :: ind
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

    nbasis = b%mesh%nedges
    nga = SIZE(b%ga)
    nsrc = SIZE(b%src)

    ALLOCATE(A(nbasis*2,nbasis*2,nga))

    WRITE(*,*) 'Name: ', TRIM(b%name)

    WRITE(*,*) 'Triangle quadrature: ', TRIM(b%qd_tri%description)
    WRITE(*,*) 'Tetrahedron quadrature: ', TRIM(b%qd_tetra%description)

    CALL print_source_info(b%src(1))

    prd => NULL()

    ! Go through all given wavelengths.
    DO n=1, b%nwl
       wl = b%sols(n)%wl
       omega = 2.0_dp*pi*c0/wl

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
       CALL sysmat_pmchwt(b, n, A)
       WRITE(*,*) 'Wall-clock time:'
       WRITE(*,*) sec_to_str(timer_end())
       WRITE(*,*) 'CPU time:'
       WRITE(*,*) sec_to_str(cpu_timer_end())

       ! Allocate memory for source/solution vector.
       ALLOCATE(b%sols(n)%x(nbasis*2, nga, nsrc))

       ! Compute excitation source vectors for each representation.
       ! Excitation is assumed to be in domain 1.
       WRITE(*,*) 'Computing source vectors'
       CALL timer_start()

       ri = b%media(b%domains(1)%medium_index)%prop(n)%ri
       DO l=1,nsrc
          CALL srcvec(b%domains(1)%mesh, b%mesh%nedges, omega, ri, b%ga, b%src(l),&
               b%qd_tri, b%sols(n)%x(:,:,l))
       END DO
       WRITE(*,*) sec_to_str(timer_end())

       ! Solve the linear system of equations for each representation,
       ! enforcing possible boundary conditions.
       ! This is O(N^3) but often in practice very fast.
       WRITE(*,*) 'Solving system'
       CALL timer_start()
       CALL cpu_timer_start()
       CALL solve_systems(b%mesh, b%ga, phdx, phdy, A, b%sols(n)%x)
       WRITE(*,*) 'Wall-clock time:'
       WRITE(*,*) sec_to_str(timer_end())
       WRITE(*,*) 'CPU time:'
       WRITE(*,*) sec_to_str(cpu_timer_end())

       ! If nonlinear media have been specified, compute second-order response.
       IF(is_nl_centrosym(b)) THEN
          WRITE(*,*) 'Solving nonlinear scattering from centrosymmetric material'

          ! Allocate memory for SH solution and auxiliary arrays.
          ALLOCATE(b%sols(n)%nlx(nbasis*2, nga, nsrc))
          ALLOCATE(src_vec(b%mesh%nedges*2))
          ALLOCATE(b%sols(n)%src_coef(b%mesh%nedges*2,SIZE(b%domains),nga,nsrc))
          ALLOCATE(epsp(b%mesh%nfaces))

          b%sols(n)%nlx(:,:,:) = 0.0_dp
          b%sols(n)%src_coef(:,:,:,:) = 0.0_dp
          src_vec(:) = 0.0_dp

          ! epsp is the selvedge-region permittivity, which should be that external
          ! to the nonlinear medium.
          ! Assumes external medium vacuum.
          !epsp(:) = eps0

          ! Setup periodic GF for current wavelength.
          ! Periodic GF must be setup again as it depends on wavelength,
          ! which now corresponds to the SH.
          DO m=1,SIZE(b%domains)
             IF(b%domains(m)%gf_index/=-1) THEN
                b%prd(b%domains(m)%gf_index)%cwl = find_closest(wl*0.5_dp,&
                     b%prd(b%domains(m)%gf_index)%coef(:)%wl)
             END IF
          END DO
          
          ! Set lattice cell phase shifts according to domain 1.
          IF(b%domains(1)%gf_index/=-1) THEN
             prd => b%prd(b%domains(1)%gf_index)
             phdx = EXP((0,1)*prd%dx*prd%coef(prd%cwl)%k0x)
             phdy = EXP((0,1)*prd%dy*prd%coef(prd%cwl)%k0y)
          END IF

          ! Go through each domain that is associated with nonlinear medium.
          DO m=1,SIZE(b%domains)
             IF(b%media(b%domains(m)%medium_index)%type/=mtype_nls) THEN
                CYCLE
             END IF

             CALL determine_epsp(b, m, n, epsp(1:b%domains(m)%mesh%nfaces))

             ! Get medium properties for current wavelength.
             mprop = b%media(b%domains(m)%medium_index)%prop(n)

             ! Get mapping from local edge indices of domain m to global
             ! edge indices, which correspond to solution coefficient indices.
             nind = b%domains(m)%mesh%nedges
             ind(1:nind) = b%domains(m)%mesh%edges(:)%parent_index

             ! Compute part of the excitation source vector and the expansion
             ! coefficients to the surface polarization in RWG basis.
             DO nf=1,nga
                DO l=1,nsrc
                   CALL nlsurf_coef(b%domains(m)%mesh, b%mesh%nedges, omega, mprop%ri, mprop%shri,&
                        epsp(1:b%domains(m)%mesh%nfaces), b%sols(n)%x(:,:,l), b%ga, nf, mprop%nls,&
                        b%qd_tri, phdx, phdy, b%sols(n)%src_coef(1:(2*nind),m,nf,l),&
                        src_vec(1:(2*nind)))

                   !CALL nlsurf_srcvec(b%domains(m)%mesh, b%mesh%nedges, omega, mprop%ri, mprop%shri,&
                   !     epsp(1:b%domains(m)%mesh%nfaces), b%sols(n)%x(:,:,l), b%ga, nf, mprop%nls,&
                   !     src_vec(1:(2*nind)))

                   ! Place the source vector elements to proper places by the use of
                   ! the edge index mappings.
                   b%sols(n)%nlx(ind(1:nind),nf,l) = b%sols(n)%nlx(ind(1:nind),nf,l) + src_vec(1:nind)
                
                   b%sols(n)%nlx(ind(1:nind)+nbasis,nf,l) = b%sols(n)%nlx(ind(1:nind)+nbasis,nf,l) +&
                        src_vec((nind+1):(2*nind))
                END DO
             END DO
          END DO

          ! Compute the system matrix for the nonlinear problem and add the remaining
          ! part to the source vector (this part depends on the submatrices).
          CALL timer_start()
          CALL sysmat_pmchwt_nls(b, n, A, b%sols(n)%src_coef, b%sols(n)%nlx)
          WRITE(*,*) sec_to_str(timer_end())

          ! Solve the linear system of equations, enforcing boundary conditions.
          WRITE(*,*) 'Solving system'
          CALL timer_start()
          CALL solve_systems(b%mesh, b%ga, phdx, phdy, A, b%sols(n)%nlx)
          WRITE(*,*) sec_to_str(timer_end())

          ! Free up memory reserved for the auxiliary arrays.
          DEALLOCATE(src_vec, epsp)
       END IF

       ! If nonlinear media have been specified, compute second-order response.
       IF(is_nl_noncentrosym(b)) THEN
          WRITE(*,*) 'Solving nonlinear scattering from non-centrosymmetric material'

          ! Allocate memory for SH solution and auxiliary arrays.
          ALLOCATE(b%sols(n)%nlx(nbasis*2, nga, nsrc))
          ALLOCATE(src_coef(b%mesh%nedges*2,SIZE(b%domains),nga,nsrc))
          ALLOCATE(src_vec2(b%mesh%nedges*2, nga, nsrc))

          b%sols(n)%nlx(:,:,:) = 0.0_dp
          src_coef(:,:,:,:) = 0.0_dp
          src_vec2(:,:,:) = 0.0_dp

          ! Setup periodic GF for current wavelength.
          ! Periodic GF must be setup again as it depends on wavelength,
          ! which now corresponds to the SH.
          DO m=1,SIZE(b%domains)
             IF(b%domains(m)%gf_index/=-1) THEN
                b%prd(b%domains(m)%gf_index)%cwl = find_closest(wl*0.5_dp,&
                     b%prd(b%domains(m)%gf_index)%coef(:)%wl)
             END IF
          END DO
          
          ! Set lattice cell phase shifts according to domain 1.
          IF(b%domains(1)%gf_index/=-1) THEN
             prd => b%prd(b%domains(1)%gf_index)
             phdx = EXP((0,1)*prd%dx*prd%coef(prd%cwl)%k0x)
             phdy = EXP((0,1)*prd%dy*prd%coef(prd%cwl)%k0y)
          END IF

          WRITE(*,*) 'Computing sources'
          CALL timer_start()

          ! Go through each domain that is associated with nonlinear medium.
          DO m=1,SIZE(b%domains)
             IF(b%media(b%domains(m)%medium_index)%type/=mtype_nlb_dipole) THEN
                CYCLE
             END IF

             ! Get medium properties for current wavelength.
             mprop = b%media(b%domains(m)%medium_index)%prop(n)

             ! Get mapping from local edge indices of domain m to global
             ! edge indices, which correspond to solution coefficient indices.
             nind = b%domains(m)%mesh%nedges
             ind(1:nind) = b%domains(m)%mesh%edges(:)%parent_index

             ! Compute the excitation source vectors for each source and representation.
             CALL srcvec_nlbulk_dipole(b%domains(m)%mesh, b%mesh%nedges, omega, mprop%ri,&
                  mprop%shri, b%sols(n)%x, b%ga, mprop%nlb, b%qd_tri, b%qd_tetra,&
                  src_vec2(1:(2*nind),:,:))
                
             DO l=1,nsrc
                ! Place the source vector elements to proper places by the use of
                ! the edge index mappings.
                b%sols(n)%nlx(ind(1:nind),:,l) = b%sols(n)%nlx(ind(1:nind),:,l)&
                     + src_vec2(1:nind,:,l)
                
                b%sols(n)%nlx(ind(1:nind)+nbasis,:,l) = b%sols(n)%nlx(ind(1:nind)+nbasis,:,l) +&
                     src_vec2((nind+1):(2*nind),:,l)
             END DO
          END DO

          WRITE(*,*) sec_to_str(timer_end())

          ! Compute the system matrix for the nonlinear problem.
          WRITE(*,*) 'Generating system matrix'
          CALL timer_start()
          CALL sysmat_pmchwt_nls(b, n, A, src_coef, b%sols(n)%nlx)
          WRITE(*,*) sec_to_str(timer_end())

          ! Solve the linear system of equations, enforcing boundary conditions.
          WRITE(*,*) 'Solving system'
          CALL timer_start()
          CALL solve_systems(b%mesh, b%ga, phdx, phdy, A, b%sols(n)%nlx)
          WRITE(*,*) sec_to_str(timer_end())

          ! Free up memory reserved for the auxiliary arrays.
          DEALLOCATE(src_vec2, src_coef)
       END IF
    END DO

    DEALLOCATE(A)

    WRITE(*,*) '--- End wavelength batch ---'

  END SUBROUTINE solve_batch

  ! Solves a linear system of equations for all group representations
  ! as described by matrices A and source vectors x.
  ! Removes and copies elements of A and x in order to enforce
  ! boundary conditions.
  SUBROUTINE solve_systems(mesh, ga, phdx, phdy, A, x)
    TYPE(mesh_container), INTENT(IN) :: mesh
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), INTENT(IN) :: phdx, phdy
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: A
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: x

    INTEGER :: dim, r, nga
    INTEGER, DIMENSION(mesh%nedges*2) :: id
    COMPLEX (KIND=dp), DIMENSION(mesh%nedges*2) :: phase

    nga = SIZE(ga)

    ! Go through all representations.
    DO r=1,nga
       ! Determine boundary conditions arising from symmetry and
       ! periodicity. BC information is stored in id and phase.
       CALL edge_bc(mesh, ga, phdx, phdy, r, id, phase)
       
       ! Enforce the boundary conditions on matrix A and source x.
       CALL resolve_system_dependencies(A(:,:,r), x(:,r,:), id, phase)
       CALL reduce_system(A(:,:,r), x(:,r,:), dim, id)
       
       ! Solve the remaining linear system of dimension dim.
       CALL solve_multi_linsys(A(1:dim,1:dim,r), x(1:dim,r,:))
       
       ! Expand the solution vector.
       ! If some elements of x were deemed zero or lin. dep.
       ! put these values into the solution vectors and make
       ! its size 2*nedges.
       CALL expand_solution(dim, id, phase, x(:,r,:))
    END DO
  END SUBROUTINE solve_systems

  ! Checks whether the batch contains nonlinear centrosymmetric materials.
  FUNCTION is_nl_centrosym(b) RESULT(res)
    TYPE(batch), INTENT(IN) :: b
    LOGICAL :: res
    INTEGER :: n

    res = .FALSE.

    DO n=1,SIZE(b%media)
       IF(b%media(n)%type==mtype_nls .OR. b%media(n)%type==mtype_nlb_nonlocal) THEN
          res = .TRUE.
          RETURN
       END IF
    END DO
  END FUNCTION is_nl_centrosym

  ! Checks whether the batch contains nonlinear non-centrosymmetric materials.
  FUNCTION is_nl_noncentrosym(b) RESULT(res)
    TYPE(batch), INTENT(IN) :: b
    LOGICAL :: res
    INTEGER :: n

    res = .FALSE.

    DO n=1,SIZE(b%media)
       IF(b%media(n)%type==mtype_nlb_dipole) THEN
          res = .TRUE.
          RETURN
       END IF
    END DO
  END FUNCTION is_nl_noncentrosym

  ! Determines the 'primed' permittivity for each face of the surface mesh.
  ! This permittivity if that of the complement domain just outside the boundary.
  SUBROUTINE determine_epsp(b, dindex, wlindex, epsp)
    TYPE(batch), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: dindex, wlindex
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: epsp
    COMPLEX (KIND=dp) :: ri
    LOGICAL :: found_face

    INTEGER :: n, m, l, parent

    DO n=1,b%domains(dindex)%mesh%nfaces
       parent = b%domains(dindex)%mesh%faces(n)%parent_index

       found_face = .FALSE.

       DO m=1,SIZE(b%domains)
          IF(m==dindex) THEN
             CYCLE
          END IF

          DO l=1,b%domains(m)%mesh%nfaces
             IF(b%domains(m)%mesh%faces(l)%parent_index==parent) THEN
                ri = b%media( b%domains(m)%medium_index )%prop(wlindex)%ri
                epsp(n) = (ri**2)*eps0

                found_face = .TRUE.

                EXIT
             END IF
          END DO

          IF(found_face) THEN
             EXIT
          END IF
       END DO

       IF(found_face==.FALSE.) THEN
          WRITE(*,*) 'Could not find epsp!'
          STOP
       END IF
    END DO
  END SUBROUTINE determine_epsp

END MODULE solver
