! MODULE: interface
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implementation of text based interface for the program.
! The interface basically help in filling the data in the batch type
! (see common.f90) in the most frequent uses.
MODULE interface
  USE solver
  USE diffr
  USE ffields
  USE cs
  USE nfpost
  USE solver_vie

  IMPLICIT NONE

CONTAINS
  SUBROUTINE read_quad(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    CHARACTER (LEN=256) :: tri_name, tetra_name

    READ(line,*) tri_name, tetra_name

    CALL delete_quad_data(b%qd_tri)
    CALL delete_quad_data(b%qd_tetra)

    b%qd_tri = tri_quad_data(TRIM(tri_name))
    b%qd_tetra = tetra_quad_data(TRIM(tetra_name))
  END SUBROUTINE read_quad

  SUBROUTINE read_ndom(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: ndom

    READ(line,*) ndom

    ALLOCATE(b%domains(1:ndom))

    WRITE(*,'(A,I0,A)') ' Allocated ', ndom, ' domains.'
  END SUBROUTINE read_ndom

  SUBROUTINE read_sdom(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: dom_index, nsurf, nvol, n
    INTEGER, DIMENSION(:), ALLOCATABLE :: surf_ids
    INTEGER, DIMENSION(:), POINTER :: vol_ids
    CHARACTER (LEN=256) :: numstr, oname

    IF(ALLOCATED(b%domains)==.FALSE.) THEN
       WRITE(*,*) 'Error: no domains allocated!'
       STOP
    END IF

    IF(ALLOCATED(b%mesh%faces)==.FALSE.) THEN
       WRITE(*,*) 'Error: master mesh has not been loaded!'
       STOP
    END IF

    READ(line,*) dom_index, nsurf

    ALLOCATE(surf_ids(1:nsurf))

    READ(line,*) dom_index, nsurf, surf_ids(1:nsurf), nvol

    vol_ids => NULL()

    IF(nvol>0) THEN
       ALLOCATE(vol_ids(1:nvol))
       
       READ(line,*) dom_index, nsurf, surf_ids(1:nsurf), nvol, vol_ids(1:nvol),&
            b%domains(dom_index)%medium_index, b%domains(dom_index)%gf_index
    ELSE
       READ(line,*) dom_index, nsurf, surf_ids(1:nsurf), nvol,&
            b%domains(dom_index)%medium_index, b%domains(dom_index)%gf_index
    END IF

    b%domains(dom_index)%mesh = extract_submesh(b%mesh, ABS(surf_ids), vol_ids)

    DO n=1,nsurf
       IF(surf_ids(n)==0) THEN
          WRITE(*,*) 'Zero is not a valid surface id!'
          STOP
       END IF

       IF(surf_ids(n)<0) THEN
          CALL invert_faces(b%domains(dom_index)%mesh, ABS(surf_ids(n)))
       END IF
    END DO

    !WRITE(numstr, '(I0)') dom_index
    !oname = TRIM(b%name) // '-' // TRIM(ADJUSTL(numstr)) // '.msh'
    !CALL save_msh(b%domains(dom_index)%mesh, b%scale, oname)

    CALL build_mesh(b%domains(dom_index)%mesh, 1.0_dp)

    DEALLOCATE(surf_ids)

    IF(ASSOCIATED(vol_ids)) THEN
       DEALLOCATE(vol_ids)
    END IF

  END SUBROUTINE read_sdom

  SUBROUTINE read_fdom(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    CHARACTER (LEN=256) :: numstr, oname
    INTEGER :: i

    CALL determine_edge_couples(b%mesh, 1D-12)
    CALL submesh_edge_connectivity(b%mesh, b%domains(:)%mesh)
    CALL orient_basis(b%mesh, b%domains(:)%mesh)

    !CALL export_mesh(TRIM(b%name) // '.pmf', b%mesh)
    !DO i=1,SIZE(b%domains)
    !   WRITE(numstr, '(I0)') i
    !   oname = TRIM(b%name) // '-' // TRIM(ADJUSTL(numstr)) // '.pmf'
    !   CALL export_mesh(oname, b%domains(i)%mesh)
    !END DO
  END SUBROUTINE read_fdom

  SUBROUTINE read_sbnd(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: id, bnd
    CHARACTER (LEN=256) :: bndname

    READ(line,*) id, bndname

    IF(TRIM(bndname)=='prdx1') THEN
       bnd = mesh_bnd_prdx1
    ELSE IF(TRIM(bndname)=='prdx2') THEN
       bnd = mesh_bnd_prdx2
    ELSE IF(TRIM(bndname)=='prdy1') THEN
       bnd = mesh_bnd_prdy1
    ELSE IF(TRIM(bndname)=='prdy2') THEN
       bnd = mesh_bnd_prdy2
    ELSE IF(TRIM(bndname)=='xplane') THEN
       bnd = mesh_bnd_xplane
    ELSE IF(TRIM(bndname)=='yplane') THEN
       bnd = mesh_bnd_yplane
    ELSE IF(TRIM(bndname)=='zplane') THEN
       bnd = mesh_bnd_zplane
    ELSE IF(TRIM(bndname)=='rz1') THEN
       bnd = mesh_bnd_rz1
    ELSE IF(TRIM(bndname)=='rz2') THEN
       bnd = mesh_bnd_rz2
    ELSE
       WRITE(*,*) 'Unrecognized boundary type!'
       STOP
    END IF

    CALL classify_edges(b%mesh, id, bnd)
  END SUBROUTINE read_sbnd

  ! Define the groups that describe the symmetry of the geometry.
  SUBROUTINE read_symm(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: nsubgroups
    TYPE(group_action), DIMENSION(:), ALLOCATABLE :: ga
    CHARACTER (LEN=3), DIMENSION(4) :: group_names
    CHARACTER (LEN=3) :: name
    INTEGER :: n, m, nr

    READ(line,*) nsubgroups
    READ(line,*) nsubgroups, group_names(1:nsubgroups)

    CALL group_id(b%ga)

    DO n=1,nsubgroups
       name = group_names(n)

       IF(TRIM(name)=='id') THEN
          CALL group_id(ga)
       ELSE IF(TRIM(name)=='mxp') THEN
          CALL group_mp(1, ga)
       ELSE IF(TRIM(name)=='myp') THEN
          CALL group_mp(2, ga)
       ELSE IF(TRIM(name)=='mzp') THEN
          CALL group_mp(3, ga)
       ELSE IF(name(1:1)=='r') THEN
          READ(name(2:3),*) nr
          CALL group_rz(nr, ga)
       ELSE
          WRITE(*,*) 'Error: unrecognized group!'
          IF(ALLOCATED(b%ga)) THEN
             DO m=1,SIZE(b%ga)
                DEALLOCATE(b%ga(m)%ef)
             END DO
             DEALLOCATE(b%ga)
          END IF
          RETURN
       END IF

       CALL product_group(b%ga, ga)

       DO m=1,SIZE(ga)
          DEALLOCATE(ga(m)%ef)
       END DO
       DEALLOCATE(ga)
    END DO

    WRITE(*,'(A,I0,A)') ' Created a symmetry group with ', SIZE(b%ga), ' elements.'

    !CALL print_group(b%ga)
  END SUBROUTINE read_symm

  ! Allocate media.
  SUBROUTINE read_nmed(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: n, i

    IF(ALLOCATED(b%sols)==.FALSE.) THEN
       WRITE(*,*) 'Error: no wavelengths set up!'
       RETURN
    END IF

    READ(line,*) n

    ALLOCATE(b%media(n))

    DO i=1,n
       ALLOCATE(b%media(i)%prop(b%nwl))
    END DO

    WRITE(*,'(A,I0,A)') ' Allocated ', n, ' media.'
  END SUBROUTINE read_nmed

  ! Set properties of medium.
  SUBROUTINE read_smed(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: mindex, n
    CHARACTER (LEN=256) :: type, method, file, file2, file3
    COMPLEX (KIND=dp) :: val
    REAL (KIND=dp), DIMENSION(3,6) :: chi2

    IF(ALLOCATED(b%media)==.FALSE.) THEN
       WRITE(*,*) 'Error: no media allocated!'
       RETURN
    END IF

    READ(line,*) mindex, type, method

    ! Set default medium type as linear.
    b%media(mindex)%type = mtype_linear

    IF(TRIM(type)=='linear') THEN
       IF(method=='file') THEN
          READ(line,*) mindex, type, method, file

          ! Fetch the refractive indices for the domain.
          DO n=1,b%nwl
             b%media(mindex)%prop(n)%ri = get_refind(TRIM(ADJUSTL(file)), b%sols(n)%wl)
             b%media(mindex)%prop(n)%shri = get_refind(TRIM(ADJUSTL(file)), b%sols(n)%wl*0.5_dp)
          END DO

          WRITE(*,'(A,I0,A,A)') 'Set medium ', mindex, ' refractive index to ', TRIM(file)
       ELSE IF(method=='value') THEN
          READ(line,*) mindex, type, method, val

          ! Set constant value for refractive index.
          DO n=1,b%nwl
             b%media(mindex)%prop(n)%ri = val
             b%media(mindex)%prop(n)%shri = val
          END DO

          WRITE(*,'(A,I0,A,"(",F6.3,",",F6.3,")")') 'Set medium ', mindex,&
               ' refractive index to ', val
       ELSE
          WRITE(*,*) 'Unrecongnized method for retrieveing medium properties!'
          RETURN
       END IF
    ELSE IF(TRIM(type)=='nlsurf') THEN
       READ(line,*) mindex, type, method

       b%media(mindex)%type = mtype_nls

       IF(method=='value') THEN
          READ(line,*) mindex, type, method, b%media(mindex)%prop(1)%nls%chi2_nnn,&
               b%media(mindex)%prop(1)%nls%chi2_ntt, b%media(mindex)%prop(1)%nls%chi2_ttn

          DO n=1,b%nwl
             b%media(mindex)%prop(n)%nls%chi2_nnn = b%media(mindex)%prop(1)%nls%chi2_nnn
             b%media(mindex)%prop(n)%nls%chi2_ntt = b%media(mindex)%prop(1)%nls%chi2_ntt
             b%media(mindex)%prop(n)%nls%chi2_ttn = b%media(mindex)%prop(1)%nls%chi2_ttn
          END DO
       ELSE
          WRITE(*,*) 'Unrecongnized method for retrieveing medium properties!'
          RETURN
       END IF
    ELSE IF(TRIM(type)=='nlbulk-nonlocal') THEN
       READ(line,*) mindex, type, method

       b%media(mindex)%type = mtype_nlb_nonlocal

       IF(method=='value') THEN
          READ(line,*) mindex, type, method, b%media(mindex)%prop(1)%nlb%delta_prime,&
               b%media(mindex)%prop(1)%nlb%gamma

          DO n=1,b%nwl
             b%media(mindex)%prop(n)%nlb%delta_prime = b%media(mindex)%prop(1)%nlb%delta_prime
             b%media(mindex)%prop(n)%nlb%gamma = b%media(mindex)%prop(1)%nlb%gamma
          END DO
       ELSE
          WRITE(*,*) 'Unrecongnized method for retrieveing medium properties!'
          RETURN
       END IF
    ELSE IF(TRIM(type)=='nlbulk-dipole') THEN
       READ(line,*) mindex, type, method

       b%media(mindex)%type = mtype_nlb_dipole

       IF(method=='file') THEN
          READ(line,*) mindex, type, method, file, file2, file3

          CALL get_matrix(TRIM(ADJUSTL(file)), chi2, 3, 6)
          b%media(mindex)%prop(1)%nlb%chi2 = chi2

          CALL get_matrix(TRIM(ADJUSTL(file2)), b%media(mindex)%prop(1)%nlb%T, 3, 3)
          CALL get_matrix(TRIM(ADJUSTL(file3)), b%media(mindex)%prop(1)%nlb%invT, 3, 3)

          WRITE(*,*) 'chi2'
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%chi2(1,:)
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%chi2(2,:)
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%chi2(3,:)

          WRITE(*,*) 'crystal to lab'
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%T(1,:)
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%T(2,:)
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%T(3,:)

          WRITE(*,*) 'lab to crystal'
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%invT(1,:)
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%invT(2,:)
          WRITE(*,*) b%media(mindex)%prop(1)%nlb%invT(3,:)

          DO n=1,b%nwl
             b%media(mindex)%prop(n)%nlb%chi2 = b%media(mindex)%prop(1)%nlb%chi2
             b%media(mindex)%prop(n)%nlb%T = b%media(mindex)%prop(1)%nlb%T
             b%media(mindex)%prop(n)%nlb%invT = b%media(mindex)%prop(1)%nlb%invT
          END DO
       ELSE
          WRITE(*,*) 'Unrecongnized method for retrieveing medium properties!'
          RETURN
       END IF
    ELSE
       WRITE(*,*) 'Unrecognized material type!'
    END IF
  END SUBROUTINE read_smed

  SUBROUTINE read_mesh(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b

    READ(line,*) b%mesh_file, b%scale

    b%mesh = load_mesh(b%mesh_file)

    CALL build_mesh(b%mesh, b%scale)
  END SUBROUTINE read_mesh

  SUBROUTINE read_wlrange(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    REAL (KIND=dp) :: wl1, wl2
    INTEGER :: n

    IF(ALLOCATED(b%sols)) THEN
       DEALLOCATE(b%sols)
    END IF

    READ(line,*) b%nwl, wl1, wl2

    ALLOCATE(b%sols(b%nwl))

    IF(b%nwl==1) THEN
       b%sols(1)%wl = wl1
    ELSE
       DO n=1,b%nwl
          b%sols(n)%wl = linterp(wl1, wl2, REAL(n-1, KIND=dp)/REAL(b%nwl-1, KIND=dp))
       END DO
    END IF
  END SUBROUTINE read_wlrange

  SUBROUTINE read_wllist(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    REAL (KIND=dp) :: wl1, wl2
    INTEGER :: n

    IF(ALLOCATED(b%sols)) THEN
       DEALLOCATE(b%sols)
    END IF

    READ(line,*) b%nwl

    ALLOCATE(b%sols(b%nwl))

    READ(line,*) b%nwl, (b%sols(n)%wl, n=1,b%nwl)
  END SUBROUTINE read_wllist

  SUBROUTINE read_scan_source(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: npt, n, m, index
    REAL (KIND=dp) :: sz, d
    TYPE(srcdata) :: src

    IF(ALLOCATED(b%src)==.FALSE.) THEN
       WRITE(*,*) 'Source must be setup before specifying scanning!'
       STOP
    END IF

    READ(line,*) sz, d, npt

    src = b%src(1)

    DEALLOCATE(b%src)

    ALLOCATE(b%src(1:(npt*npt)))

    ! Construct sources so that in column major matrix R
    ! of results R(1,1) corresponds to top left position.
    DO n=1,npt
       DO m=1,npt
          index = (n-1)*npt + m
          b%src(index) = src
          b%src(index)%pos = (/-d/2+d*REAL(n-1,KIND=dp)/(npt-1), d/2-d*REAL(m-1,KIND=dp)/(npt-1), sz/)
       END DO
    END DO

    WRITE(*,*) 'Source scanning was setup'

  END SUBROUTINE read_scan_source

  SUBROUTINE read_move_source(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: index
    REAL (KIND=dp), DIMENSION(3) :: pos
    TYPE(srcdata) :: src

    IF(ALLOCATED(b%src)==.FALSE.) THEN
       WRITE(*,*) 'Source must be setup before specifying move transform!'
       STOP
    END IF

    READ(line,*) index, pos(1:3)

    b%src(index)%pos = pos
  END SUBROUTINE read_move_source

  ! Allocate source definitions.
  SUBROUTINE read_nsrc(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: n, i

    READ(line,*) n

    ALLOCATE(b%src(n))

    WRITE(*,'(A,I0,A)') ' Allocated ', n, ' sources.'
  END SUBROUTINE read_nsrc

  SUBROUTINE read_osrc(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: n
    REAL (KIND=dp), DIMENSION(3) :: pos

    READ(line,*) n, pos(1:3)

    b%src(n)%pos = pos

    WRITE(*,*) 'Offset source ', n, ' by ', pos
    
  END SUBROUTINE read_osrc

  SUBROUTINE read_ssrc(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    CHARACTER (LEN=32) :: token
    CHARACTER (LEN=256) :: solname, src_meshname, nfocus
    CHARACTER (LEN=3) :: ext
    INTEGER :: nga, nfrags, index

    IF(ALLOCATED(b%src)==.FALSE.) THEN
       WRITE(*,*) 'No sources allocated!'
       STOP
    END IF

    READ(line,*) index

    b%src(index)%pos(:) = 0.0_dp

    READ(line,*) index, token

    IF(token=='pw') THEN
       b%src(index)%type = src_pw
       READ(line,*) index, token, b%src(index)%theta, b%src(index)%phi, b%src(index)%psi
       b%src(index)%theta = b%src(index)%theta*degtorad
       b%src(index)%phi = b%src(index)%phi*degtorad
       b%src(index)%psi = b%src(index)%psi*degtorad
    ELSE IF(token=='pw_elliptic') THEN
       b%src(index)%type = src_pw_elliptic
       READ(line,*) index, token, b%src(index)%theta, b%src(index)%phi, b%src(index)%phase
       b%src(index)%theta = b%src(index)%theta*degtorad
       b%src(index)%phi = b%src(index)%phi*degtorad
    ELSE IF(token=='focus_rad') THEN
       b%src(index)%type = src_focus_rad
       READ(line,*) index, token, b%src(index)%focal, b%src(index)%waist, b%src(index)%napr, nfocus
    ELSE IF(token=='focus_x') THEN
       b%src(index)%type = src_focus_x
       READ(line,*) index, token, b%src(index)%focal, b%src(index)%waist, b%src(index)%napr, nfocus
    ELSE IF(token=='focus_y') THEN
       b%src(index)%type = src_focus_y
       READ(line,*) index, token, b%src(index)%focal, b%src(index)%waist, b%src(index)%napr, nfocus
    ELSE IF(token=='focus_hg01') THEN
       b%src(index)%type = src_focus_hg01
       READ(line,*) index, token, b%src(index)%focal, b%src(index)%waist, b%src(index)%napr, nfocus
    ELSE IF(token=='focus_azimut') THEN
       b%src(index)%type = src_focus_azim
       READ(line,*) index, token, b%src(index)%focal, b%src(index)%waist, b%src(index)%napr, nfocus
    ELSE IF(token=='dipole') THEN
       b%src(index)%type = src_dipole
       READ(line,*) index, token, b%src(index)%pos, b%src(index)%dmom
    ELSE
       WRITE(*,*) 'Invalid source type!'
    END IF

    IF(b%src(index)%type==src_focus_rad .OR. b%src(index)%type==src_focus_x .OR.&
         b%src(index)%type==src_focus_y .OR. b%src(index)%type==src_focus_hg01 .OR.&
         b%src(index)%type==src_focus_azim) THEN
       IF(nfocus=='true') THEN
          b%src(index)%nfocus = .TRUE.
       ELSE IF(nfocus=='false') THEN
          b%src(index)%nfocus = .FALSE.
       ELSE
          WRITE(*,*) 'Unrecognized source argument value!'
          STOP
       END IF
    END IF
  END SUBROUTINE read_ssrc

  SUBROUTINE read_npgf(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: n

    READ(line,*) n

    ALLOCATE(b%prd(1:n))
  END SUBROUTINE read_npgf

  SUBROUTINE read_ipgw(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    CHARACTER (LEN=16) :: tag
    CHARACTER (LEN=256) :: filename
    INTEGER :: index

    READ(line,*) index, filename

    CALL import_pgfw_interp_1d(TRIM(filename), b%prd(index))
  END SUBROUTINE read_ipgw

  SUBROUTINE read_fres(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: fid=10, iovar, n, srcindex, drindex, dtindex, nf
    REAL (KIND=dp) :: wl, omega, zr, zt
    COMPLEX (KIND=dp) :: ri, r, t
    COMPLEX (KIND=dp), DIMENSION(3) :: er, hr, et, ht
    REAL (KIND=dp), DIMENSION(3) :: pol
    TYPE(prdnfo), POINTER :: prd
    CHARACTER (LEN=256) :: oname, numstr
    TYPE(srcdata) :: src

    ! Symmetry utilization not supported. Use the first fragment.
    nf = 1

    READ(line,*) srcindex, drindex, dtindex, zr, zt

    src = b%src(srcindex)
    pol = get_pol(src%theta, src%phi, src%psi)

    WRITE(numstr, '(A,I0)') '-s', srcindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-fresnel.dat'

    OPEN(fid, FILE=oname, ACTION='WRITE', IOSTAT=iovar)
    IF(iovar>0) THEN
       WRITE(*,*) 'Could not open output file for Fresnel data!'
       STOP
    END IF

    DO n=1,b%nwl
       wl = b%sols(n)%wl
       omega = 2.0_dp*pi*c0/wl

       ! Reflected field.
       ri = b%media(b%domains(drindex)%medium_index)%prop(n)%ri
       prd => b%prd(b%domains(drindex)%gf_index)
       prd%cwl = find_closest(wl, prd%coef(:)%wl)

       CALL diff_fields(b%domains(drindex)%mesh, b%ga, 1, b%sols(n)%x(:,nf,srcindex), b%mesh%nedges, omega, ri,&
            prd, (/0.0_dp,0.0_dp,zr/), 0, 0, b%qd_tri, er, hr)

       ! Transmitted field.
       ri = b%media(b%domains(dtindex)%medium_index)%prop(n)%ri
       prd => b%prd(b%domains(dtindex)%gf_index)
       prd%cwl = find_closest(wl, prd%coef(:)%wl)

       CALL diff_fields(b%domains(dtindex)%mesh, b%ga, 1, b%sols(n)%x(:,nf,srcindex), b%mesh%nedges, omega, ri,&
            prd, (/0.0_dp,0.0_dp,zt/), 0, 0, b%qd_tri, et, ht)

       ! Fresnel coefficients: project field vectors to incident wave polarization vector.
       r = dotc(CMPLX(pol,KIND=dp), er)
       t = dotc(CMPLX(pol,KIND=dp), et)

       WRITE(fid,*) wl, REAL(r), AIMAG(r), REAL(t), AIMAG(t)
    END DO

    CLOSE(fid)
  END SUBROUTINE read_fres

  SUBROUTINE read_diff(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: fid=10, iovar, n, dindex, orderx, ordery, srcindex
    REAL (KIND=dp) :: wl, omega, polangle
    REAL (KIND=dp), DIMENSION(b%nwl) :: irr
    COMPLEX (KIND=dp) :: ri, ri_inc
    TYPE(prdnfo), POINTER :: prd
    CHARACTER (LEN=256) :: oname, numstr, polstr
    LOGICAL :: polarize

    READ(line,*) srcindex, dindex, orderx, ordery, polstr, polangle

    WRITE(*,'(A,I0,A,I0,A)') 'Computing diffracted irradiance to order (', orderx, ',', ordery, ')'

    IF(TRIM(polstr)=='true') THEN
       polarize = .TRUE.
       WRITE(*,*) 'Using a polarizer at angle', polangle, ' deg'
    ELSE
       polarize = .FALSE.
    END IF

    ! Transform from degrees to radians.
    polangle = polangle*pi/180

    DO n=1,b%nwl
       wl = b%sols(n)%wl

       omega = 2.0_dp*pi*c0/wl
       ri = b%media(b%domains(dindex)%medium_index)%prop(n)%ri
       ri_inc = b%media(b%domains(1)%medium_index)%prop(n)%ri

       prd => b%prd(b%domains(dindex)%gf_index)

       prd%cwl = find_closest(wl, prd%coef(:)%wl)

       irr(n) = diff_irradiance(b%domains(dindex)%mesh, b%ga, dindex==1, b%src(srcindex),&
            b%sols(n)%x(:,:,srcindex), b%mesh%nedges, omega, ri, ri_inc, prd,&
            orderx, ordery, b%qd_tri, polarize, polangle)

       WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
    END DO

    IF(orderx==0 .AND. ordery==0) THEN
       WRITE(numstr, '(A,I0)') '-s', srcindex
    ELSE
       WRITE(numstr, '(A,I0,A,I0,I0)') '-s', srcindex, '-o', orderx, ordery
    END IF

    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.dif'

    OPEN(fid, FILE=oname, ACTION='WRITE', IOSTAT=iovar)
    IF(iovar>0) THEN
       WRITE(*,*) 'Could not open output file for diffraction data!'
       STOP
    END IF

    DO n=1,b%nwl
       wl = b%sols(n)%wl

       WRITE(fid,*) wl, irr(n)
    END DO

    CLOSE(fid)

    IF(ALLOCATED(b%sols(1)%nlx)) THEN
       DO n=1,b%nwl
          ! SHG
          wl = b%sols(n)%wl*0.5_dp
          
          omega = 2.0_dp*pi*c0/wl
          ri = b%media(b%domains(dindex)%medium_index)%prop(n)%shri
          ri_inc = b%media(b%domains(1)%medium_index)%prop(n)%ri

          prd => b%prd(b%domains(dindex)%gf_index)
          
          prd%cwl = find_closest(wl, prd%coef(:)%wl)
          
          irr(n) = diff_irradiance(b%domains(dindex)%mesh, b%ga, .FALSE., b%src(srcindex),&
               b%sols(n)%nlx(:,:,srcindex), b%mesh%nedges, omega, ri, ri_inc, prd,&
               orderx, ordery, b%qd_tri, polarize, polangle)
          
          WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
       END DO

       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-sh.dif'

       OPEN(fid, FILE=oname, ACTION='WRITE', IOSTAT=iovar)
       IF(iovar>0) THEN
          WRITE(*,*) 'Could not open output file for diffraction data!'
          STOP
       END IF
       
       DO n=1,b%nwl
          wl = b%sols(n)%wl
          
          WRITE(fid,*) wl, irr(n)
       END DO
       
       CLOSE(fid)
    END IF

  END SUBROUTINE read_diff

  SUBROUTINE read_trns(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: fid=10, iovar, n, dindexT, dindexR, srcindex
    REAL (KIND=dp) :: wl, omega, z0T, z0R
    REAL (KIND=dp), DIMENSION(b%nwl,3) :: data
    REAL (KIND=dp), DIMENSION(b%nwl,SIZE(b%src)) :: Tpower, Rpower
    COMPLEX (KIND=dp) :: ri, ri_inc
    TYPE(prdnfo), POINTER :: prd
    CHARACTER (LEN=64) :: addsrc_str
    CHARACTER (LEN=256) :: oname, numstr
    LOGICAL :: addsrc

    READ(line,*) dindexT, z0T, addsrc_str, dindexR, z0R

    IF(TRIM(addsrc_str)=='true') THEN
       addsrc = .TRUE.
    ELSE
       addsrc = .FALSE.
    END IF

    WRITE(*,*) 'Computing transmittance'

    DO n=1,b%nwl
       wl = b%sols(n)%wl
       omega = 2.0_dp*pi*c0/wl
       ri_inc = b%media(b%domains(1)%medium_index)%prop(n)%ri

       ri = b%media(b%domains(dindexT)%medium_index)%prop(n)%ri
       prd => b%prd(b%domains(dindexT)%gf_index)
       prd%cwl = find_closest(wl, prd%coef(:)%wl)

       Tpower(n,:) = transmittance(b%domains(dindexT)%mesh, b%ga, addsrc, b%src,&
            b%sols(n)%x, b%mesh%nedges, omega, ri, ri_inc, prd, z0T, -1.0_dp, b%qd_tri)

       ri = b%media(b%domains(dindexR)%medium_index)%prop(n)%ri
       prd => b%prd(b%domains(dindexR)%gf_index)
       prd%cwl = find_closest(wl, prd%coef(:)%wl)

       Rpower(n,:) = transmittance(b%domains(dindexR)%mesh, b%ga, .FALSE., b%src,&
            b%sols(n)%x, b%mesh%nedges, omega, ri, ri_inc, prd, z0R, 1.0_dp, b%qd_tri)

       WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
    END DO

    data(:,1) = b%sols(:)%wl

    DO n=1,SIZE(b%src)
       data(:,2) = Tpower(:,n)
       data(:,3) = Rpower(:,n)

       WRITE(numstr, '(A,I0)') '-s', n
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-trns.dat'

       CALL write_data(TRIM(oname), data)
    END DO

    IF(ALLOCATED(b%sols(1)%nlx)) THEN
       DO n=1,b%nwl
          ! SH wavelength.
          wl = b%sols(n)%wl/2
          omega = 2.0_dp*pi*c0/wl
          ri_inc = b%media(b%domains(1)%medium_index)%prop(n)%ri
                    
          ri = b%media(b%domains(dindexT)%medium_index)%prop(n)%shri
          prd => b%prd(b%domains(dindexT)%gf_index)
          prd%cwl = find_closest(wl, prd%coef(:)%wl)
          
          Tpower(n,:) = transmittance(b%domains(dindexT)%mesh, b%ga, .FALSE., b%src,&
               b%sols(n)%nlx, b%mesh%nedges, omega, ri, ri_inc, prd, z0T, -1.0_dp, b%qd_tri)
          
          ri = b%media(b%domains(dindexR)%medium_index)%prop(n)%shri
          prd => b%prd(b%domains(dindexR)%gf_index)
          prd%cwl = find_closest(wl, prd%coef(:)%wl)
          
          Rpower(n,:) = transmittance(b%domains(dindexR)%mesh, b%ga, .FALSE., b%src,&
               b%sols(n)%nlx, b%mesh%nedges, omega, ri, ri_inc, prd, z0R, 1.0_dp, b%qd_tri)
          
          WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
       END DO

       DO n=1,SIZE(b%src)
          data(:,2) = Tpower(:,n)
          data(:,3) = Rpower(:,n)
          
          WRITE(numstr, '(A,I0)') '-s', n
          oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-trns-sh.dat'
          
          CALL write_data(TRIM(oname), data)
       END DO
       
    END IF
  END SUBROUTINE read_trns

  SUBROUTINE read_nfms(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: wlindex, srcindex, dindex
    CHARACTER (LEN=256) :: oname, numstr
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri

    READ(line,*) wlindex, srcindex, dindex

    WRITE(numstr, '(A,I0,A,I0,A,I0)') '-wl', wlindex, '-s', srcindex, '-d', dindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.msh'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri
       
    CALL field_mesh(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
         b%sols(wlindex)%x(:,:,srcindex), b%ga, omega, ri)

    !CALL gradPnls_mesh('gradPnls.msh', b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
    !     b%sols(wlindex)%x(:,:,srcindex), b%ga, omega, ri)

    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-sh.msh'
       
       ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%shri
       
       CALL field_mesh(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
            b%sols(wlindex)%nlx(:,:,srcindex), b%ga, 2.0_dp*omega, ri)

    END IF
  END SUBROUTINE read_nfms

  SUBROUTINE read_nfpl(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    TYPE(nfield_plane) :: nfplane
    INTEGER :: wlindex, srcindex, dindex, n1, n2, ns
    CHARACTER (LEN=256) :: oname, numstr, tag
    REAL (KIND=dp), DIMENSION(3) :: v1, v2, origin
    REAL (KIND=dp) :: omega, wl
    COMPLEX (KIND=dp) :: ri
    COMPLEX (KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: ef, hf
    TYPE(prdnfo), POINTER :: prd

    READ(line,*) wlindex, dindex, origin, v1, v2, n1, n2, tag

    ALLOCATE(ef(3,SIZE(b%src),n2,n1), hf(3,SIZE(b%src),n2,n1))

    wl = b%sols(wlindex)%wl
    omega = 2.0_dp*pi*c0/wl
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri

    IF(b%domains(dindex)%gf_index/=-1) THEN
       prd => b%prd(b%domains(dindex)%gf_index)
       prd%cwl = find_closest(wl, prd%coef(:)%wl)
    ELSE
       prd => NULL()
    END IF

    CALL field_plane(b%domains(dindex)%mesh, b%mesh%nedges, b%sols(wlindex)%x, b%ga,&
         omega, ri, prd, dindex==1, b%src, b%qd_tri, origin, v1, v2, n1, n2, ef, hf)

    DO ns=1,SIZE(b%src)
       WRITE(numstr, '(A,I0,A,I0,A,I0,A,A)') '-wl', wlindex, '-s', ns, '-d', dindex, '-', TRIM(ADJUSTL(tag))

       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ex-re.dat', REAL(ef(1,ns,:,:),KIND=dp))
       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ey-re.dat', REAL(ef(2,ns,:,:),KIND=dp))
       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ez-re.dat', REAL(ef(3,ns,:,:),KIND=dp))

       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ex-im.dat', AIMAG(ef(1,ns,:,:)))
       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ey-im.dat', AIMAG(ef(2,ns,:,:)))
       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ez-im.dat', AIMAG(ef(3,ns,:,:)))

       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hx-re.dat', REAL(hf(1,ns,:,:),KIND=dp))
       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hy-re.dat', REAL(hf(2,ns,:,:),KIND=dp))
       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hz-re.dat', REAL(hf(3,ns,:,:),KIND=dp))

       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hx-im.dat', AIMAG(hf(1,ns,:,:)))
       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hy-im.dat', AIMAG(hf(2,ns,:,:)))
       CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hz-im.dat', AIMAG(hf(3,ns,:,:)))
    END DO

    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       wl = b%sols(wlindex)%wl/2 ! SHG wavelength
       !omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl   mistake found 9.9.2014
       omega = 2.0_dp*pi*c0/wl
       ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%shri
       
       IF(b%domains(dindex)%gf_index/=-1) THEN
          prd => b%prd(b%domains(dindex)%gf_index)
          prd%cwl = find_closest(wl, prd%coef(:)%wl)
       ELSE
          prd => NULL()
       END IF
       
       CALL field_plane(b%domains(dindex)%mesh, b%mesh%nedges, b%sols(wlindex)%nlx, b%ga,&
            omega, ri, prd, .FALSE., b%src, b%qd_tri, origin, v1, v2, n1, n2, ef, hf)
       
       DO ns=1,SIZE(b%src)
          WRITE(numstr, '(A,I0,A,I0,A,I0,A,A,A)') '-wl', wlindex, '-s', ns, '-d', dindex, '-', TRIM(ADJUSTL(tag)), '-sh'
          
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ex-re.dat', REAL(ef(1,ns,:,:),KIND=dp))
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ey-re.dat', REAL(ef(2,ns,:,:),KIND=dp))
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ez-re.dat', REAL(ef(3,ns,:,:),KIND=dp))
          
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ex-im.dat', AIMAG(ef(1,ns,:,:)))
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ey-im.dat', AIMAG(ef(2,ns,:,:)))
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-ez-im.dat', AIMAG(ef(3,ns,:,:)))
          
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hx-re.dat', REAL(hf(1,ns,:,:),KIND=dp))
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hy-re.dat', REAL(hf(2,ns,:,:),KIND=dp))
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hz-re.dat', REAL(hf(3,ns,:,:),KIND=dp))
          
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hx-im.dat', AIMAG(hf(1,ns,:,:)))
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hy-im.dat', AIMAG(hf(2,ns,:,:)))
          CALL write_data(TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-hz-im.dat', AIMAG(hf(3,ns,:,:)))
       END DO
    END IF

    DEALLOCATE(ef, hf)

  END SUBROUTINE read_nfpl

  SUBROUTINE read_test(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    TYPE(nfield_plane) :: nfplane
    INTEGER :: wlindex, srcindex, dindex, nind
    CHARACTER (LEN=256) :: oname, numstr, addsrcstr, meshname
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    TYPE(mesh_container) :: extmesh
    LOGICAL :: addsrc
    TYPE(prdnfo), POINTER :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: nlx
    INTEGER, DIMENSION(b%mesh%nedges) :: ind
    INTEGER :: nd, na

    DO nd=1,SIZE(b%domains)
       WRITE(*,*) 'domain ', nd

       DO na=1,SIZE(b%ga)
          WRITE(*,'(A,I0,A)',ADVANCE='NO') 'action ', na, ': '
          IF(BTEST(b%ga(na)%genbits,gid_mxp)) THEN
             WRITE(*,'(A)',ADVANCE='NO') 'mxp '
          END IF
          IF(BTEST(b%ga(na)%genbits,gid_myp)) THEN
             WRITE(*,'(A)',ADVANCE='NO') 'myp '
          END IF
          IF(BTEST(b%ga(na)%genbits,gid_mzp)) THEN
             WRITE(*,'(A)',ADVANCE='NO') 'mzp '
          END IF
          IF(BTEST(b%ga(na)%genbits,gid_rz)) THEN
             WRITE(*,'(A)',ADVANCE='NO') 'rz '
          END IF

          IF(admissible_ga(b%domains(nd)%mesh, b%ga(na), nd==1)) THEN
             WRITE(*,'(A)') ', admissible'
          ELSE
             WRITE(*,'(A)') ', not admissible'
          END IF
       END DO
    END DO

    RETURN

    ! Extinction theorem

    READ(line,*) wlindex, srcindex, dindex, meshname

    extmesh = load_mesh(meshname)

    CALL scale_mesh(extmesh, b%scale)
    CALL compute_mesh_essentials(extmesh)

    WRITE(numstr, '(A,I0,A,I0,A,I0)') '-wl', wlindex, '-s', srcindex, '-d', dindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.msh'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl

    ! Internal FF

    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-fftest.msh'
    
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri
    
    IF(b%domains(dindex)%gf_index>0) THEN
       prd => b%prd(b%domains(dindex)%gf_index)
       prd%cwl = find_closest(b%sols(wlindex)%wl, prd%coef(:)%wl)
    ELSE
       prd => NULL()
    END IF

    CALL field_external_mesh(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
         b%sols(wlindex)%x(:,:,srcindex), b%ga, omega, ri, prd,&
         .FALSE., b%src(srcindex), extmesh, b%qd_tri)

    ! SH polarization
    
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-shtest.msh'
    
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%shri
    
    IF(b%domains(dindex)%gf_index>0) THEN
       prd => b%prd(b%domains(dindex)%gf_index)
       prd%cwl = find_closest(b%sols(wlindex)%wl*0.5_dp, prd%coef(:)%wl)
    ELSE
       prd => NULL()
    END IF
    
    ALLOCATE(nlx(SIZE(b%sols(wlindex)%nlx,1), SIZE(b%sols(wlindex)%nlx,2),&
         SIZE(b%sols(wlindex)%nlx,3)))
    
    nlx(:,:,:) = 0.0_dp
    
    IF(b%media(b%domains(dindex)%medium_index)%type==mtype_nls) THEN
       ! Get mapping from local edge indices of domain dindex to global
       ! edge indices, which correspond to solution coefficient indices.
       nind = b%domains(dindex)%mesh%nedges
       ind(1:nind) = b%domains(dindex)%mesh%edges(:)%parent_index
       
       nlx(ind(1:nind),:,:) = b%sols(wlindex)%src_coef((nind+1):(2*nind),dindex,:,:)
       
       nlx(ind(1:nind)+b%mesh%nedges,:,:) = b%sols(wlindex)%src_coef(1:nind,dindex,:,:)
    ELSE
       nlx(:,:,:) = b%sols(wlindex)%nlx(:,:,:)
    END IF
    
    CALL field_external_mesh(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
         nlx(:,:,srcindex), b%ga, 2.0_dp*omega, ri, prd,&
         .FALSE., b%src(srcindex), extmesh, b%qd_tri)
    
    DEALLOCATE(nlx)

    CALL delete_mesh(extmesh)
  END SUBROUTINE read_test

  SUBROUTINE read_nfem(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    TYPE(nfield_plane) :: nfplane
    INTEGER :: wlindex, srcindex, dindex, nind
    CHARACTER (LEN=256) :: oname, numstr, addsrcstr, meshname
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    TYPE(mesh_container) :: extmesh
    LOGICAL :: addsrc
    TYPE(prdnfo), POINTER :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: nlx
    INTEGER, DIMENSION(b%mesh%nedges) :: ind

    READ(line,*) wlindex, srcindex, dindex, meshname, addsrcstr

    extmesh = load_mesh(meshname)

    CALL scale_mesh(extmesh, b%scale)
    CALL compute_mesh_essentials(extmesh)

    addsrc = (TRIM(addsrcstr)=='true')

    WRITE(numstr, '(A,I0,A,I0,A,I0)') '-wl', wlindex, '-s', srcindex, '-d', dindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.msh'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri

    IF(b%domains(dindex)%gf_index>0) THEN
       prd => b%prd(b%domains(dindex)%gf_index)
       prd%cwl = find_closest(b%sols(wlindex)%wl, prd%coef(:)%wl)
    ELSE
       prd => NULL()
    END IF

    CALL field_external_mesh(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
         b%sols(wlindex)%x(:,:,srcindex), b%ga, omega, ri, prd,&
         addsrc, b%src(srcindex), extmesh, b%qd_tri)

    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-sh.msh'
       
       ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%shri

       IF(b%domains(dindex)%gf_index>0) THEN
          prd => b%prd(b%domains(dindex)%gf_index)
          prd%cwl = find_closest(b%sols(wlindex)%wl*0.5_dp, prd%coef(:)%wl)
       ELSE
          prd => NULL()
       END IF

       ALLOCATE(nlx(SIZE(b%sols(wlindex)%nlx,1), SIZE(b%sols(wlindex)%nlx,2),&
            SIZE(b%sols(wlindex)%nlx,3)))

       nlx(:,:,:) = 0.0_dp

       IF(b%media(b%domains(dindex)%medium_index)%type==mtype_nls) THEN
          ! Get mapping from local edge indices of domain dindex to global
          ! edge indices, which correspond to solution coefficient indices.
          nind = b%domains(dindex)%mesh%nedges
          ind(1:nind) = b%domains(dindex)%mesh%edges(:)%parent_index

          nlx(ind(1:nind),:,:) = b%sols(wlindex)%src_coef((nind+1):(2*nind),dindex,:,:)&
               + b%sols(wlindex)%nlx(ind(1:nind),:,:)

          nlx(ind(1:nind)+b%mesh%nedges,:,:) = b%sols(wlindex)%src_coef(1:nind,dindex,:,:)&
               + b%sols(wlindex)%nlx(ind(1:nind)+b%mesh%nedges,:,:)
       ELSE
          nlx(:,:,:) = b%sols(wlindex)%nlx(:,:,:)
       END IF
       
       CALL field_external_mesh(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
            nlx(:,:,srcindex), b%ga, 2.0_dp*omega, ri, prd,&
            .FALSE., b%src(srcindex), extmesh, b%qd_tri)

       DEALLOCATE(nlx)
    END IF

    CALL delete_mesh(extmesh)
  END SUBROUTINE read_nfem

  SUBROUTINE read_nfdm(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    TYPE(nfield_plane) :: nfplane
    INTEGER :: wlindex, srcindex, dindex
    CHARACTER (LEN=256) :: oname, numstr, addsrcstr
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    LOGICAL :: addsrc
    TYPE(prdnfo), POINTER :: prd
    REAL (KIND=dp), DIMENSION(3) :: origin, dsize
    INTEGER, DIMENSION(3) :: npoints

    READ(line,*) wlindex, srcindex, dindex, addsrcstr, origin, dsize, npoints

    addsrc = (TRIM(addsrcstr)=='true')

    WRITE(numstr, '(A,I0,A,I0,A,I0)') '-wl', wlindex, '-s', srcindex, '-d', dindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.msh'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri

    IF(b%domains(dindex)%gf_index>0) THEN
       prd => b%prd(b%domains(dindex)%gf_index)
       prd%cwl = find_closest(b%sols(wlindex)%wl, prd%coef(:)%wl)
    ELSE
       prd => NULL()
    END IF
    
    CALL field_domain(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
         b%sols(wlindex)%x(:,:,srcindex), b%ga, omega, ri, prd, b%src(srcindex), addsrc,&
         origin, dsize, npoints, b%qd_tri)

    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-sh.msh'
       
       ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%shri

       IF(b%domains(dindex)%gf_index>0) THEN
          prd%cwl = find_closest(b%sols(wlindex)%wl/2, prd%coef(:)%wl)
       END IF
       
       CALL field_domain(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
            b%sols(wlindex)%nlx(:,:,srcindex), b%ga, 2.0_dp*omega, ri, prd, b%src(srcindex), .FALSE.,&
            origin, dsize, npoints, b%qd_tri)
    END IF
  END SUBROUTINE read_nfdm

  SUBROUTINE read_nfst(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    TYPE(nfield_plane) :: nfplane
    INTEGER :: wlindex, srcindex, dindex
    CHARACTER (LEN=256) :: oname, numstr, addsrcstr
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    LOGICAL :: addsrc
    TYPE(prdnfo), POINTER :: prd

    READ(line,*) wlindex, srcindex, dindex, addsrcstr

    addsrc = (TRIM(addsrcstr)=='true')

    IF(addsrc) THEN
       WRITE(*,*) 'Adding source'
    END IF

    WRITE(numstr, '(A,I0,A,I0,A,I0)') '-wl', wlindex, '-s', srcindex, '-d', dindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-stream.msh'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri

    IF(b%domains(dindex)%gf_index>0) THEN
       prd => b%prd(b%domains(dindex)%gf_index)
       prd%cwl = find_closest(b%sols(wlindex)%wl, prd%coef(:)%wl)
    ELSE
       prd => NULL()
    END IF

    CALL field_stream(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
         b%sols(wlindex)%x(:,:,srcindex), b%ga, omega, ri, prd, b%src(srcindex), addsrc, b%qd_tri)
  END SUBROUTINE read_nfst

  ! Computes a 2D image, where pixels correspond to source positions.
  ! Data may be plotted direcly with MATLAB's scimage, whence x-axis
  ! points right and y-axis points up.
  SUBROUTINE read_rcs2(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: iovar, npt, wlindex, n
    REAL (KIND=dp) :: omega, theta1, theta2
    COMPLEX (KIND=dp) :: ri
    REAL (KIND=dp), DIMENSION(SIZE(b%src)) :: scatp
    CHARACTER (LEN=256) :: oname, numstr

    READ(line,*) wlindex, theta1, theta2

    WRITE(*,*) 'Computing scattered power'

    CALL timer_start()

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(1)%medium_index)%prop(wlindex)%ri

    npt = NINT(SQRT(REAL(SIZE(b%src))))

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(b,omega,ri,scatp,wlindex,theta1,theta2)&
    !$OMP PRIVATE(n)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,SIZE(b%src)
       !theta_max = ASIN(b%src(n)%napr/REAL(ri,KIND=dp))

       CALL rcs_solangle(b%domains(1)%mesh, b%mesh%nedges, omega, ri, b%ga,&
            b%sols(wlindex)%x(:,:,n), theta1, theta2, b%qd_tri, scatp(n))
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    WRITE(numstr, '(I0)') wlindex
    oname = TRIM(b%name) // '-wl' // TRIM(ADJUSTL(numstr)) // '-scan.dat'

    CALL write_data(oname, RESHAPE(scatp,(/npt,npt/)))

    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       oname = TRIM(b%name) // '-wl' // TRIM(ADJUSTL(numstr)) // '-scan-sh.dat'
       
       ! Second-harmonic frequency.
       ri = b%media(b%domains(1)%medium_index)%prop(wlindex)%shri

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(b,omega,ri,scatp,wlindex,theta1,theta2)&
       !$OMP PRIVATE(n)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,SIZE(b%src)
          !theta_max = ASIN(b%src(n)%napr/REAL(ri,KIND=dp))

          CALL rcs_solangle(b%domains(1)%mesh, b%mesh%nedges, 2.0_dp*omega, ri, b%ga,&
               b%sols(wlindex)%nlx(:,:,n), theta1, theta2, b%qd_tri, scatp(n))
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       CALL write_data(oname, RESHAPE(scatp,(/npt,npt/)))

    END IF

    WRITE(*,*) sec_to_str(timer_end())

  END SUBROUTINE read_rcs2

  ! Compute RCS by integrating over boundaries of a selected domain.
  SUBROUTINE read_rcs3(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: iovar, ntheta_rcs, nphi_rcs, wlindex, srcindex, dindex
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: rcsdata
    CHARACTER (LEN=256) :: oname, numstr

    READ(line,*) wlindex, srcindex, dindex, ntheta_rcs, nphi_rcs

    WRITE(*,*) 'Computing radar cross-sections'
    CALL timer_start()

    ALLOCATE(rcsdata(1:ntheta_rcs,1:nphi_rcs))

    WRITE(numstr, '(A,I0,A,I0)') '-wl', wlindex, '-s', srcindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.rcs'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri

    CALL rcs(b%domains(dindex)%mesh, b%mesh%nedges, omega, ri, b%ga, b%sols(wlindex)%x(:,:,srcindex),&
         ntheta_rcs, nphi_rcs, b%qd_tri, rcsdata)
    CALL write_data(oname, rcsdata)

    ! Nonlinear radar cross-sections.
    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-sh.rcs'
       
       ! Second-harmonic frequency.
       ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%shri
       
       CALL rcs(b%domains(dindex)%mesh, b%mesh%nedges, 2.0_dp*omega, ri, b%ga,&
            b%sols(wlindex)%nlx(:,:,srcindex), ntheta_rcs, nphi_rcs, b%qd_tri, rcsdata)
       CALL write_data(oname, rcsdata)

    END IF

    DEALLOCATE(rcsdata)

    WRITE(*,*) sec_to_str(timer_end())

  END SUBROUTINE read_rcs3

  SUBROUTINE read_rcst(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: iovar, ntheta_rcs, nphi_rcs, wlindex, srcindex
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: rcsdata
    CHARACTER (LEN=256) :: oname, numstr

    READ(line,*) wlindex, srcindex, ntheta_rcs, nphi_rcs

    WRITE(*,*) 'Computing radar cross-sections'
    CALL timer_start()

    ALLOCATE(rcsdata(1:ntheta_rcs,1:nphi_rcs))

    WRITE(numstr, '(A,I0,A,I0)') '-wl', wlindex, '-s', srcindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.rcs'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(1)%medium_index)%prop(wlindex)%ri

    CALL rcs(b%domains(1)%mesh, b%mesh%nedges, omega, ri, b%ga, b%sols(wlindex)%x(:,:,srcindex),&
         ntheta_rcs, nphi_rcs, b%qd_tri, rcsdata)
    CALL write_data(oname, rcsdata)

    ! Nonlinear radar cross-sections.
    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-sh.rcs'
       
       ! Second-harmonic frequency.
       ri = b%media(b%domains(1)%medium_index)%prop(wlindex)%shri
       
       CALL rcs(b%domains(1)%mesh, b%mesh%nedges, 2.0_dp*omega, ri, b%ga,&
            b%sols(wlindex)%nlx(:,:,srcindex), ntheta_rcs, nphi_rcs, b%qd_tri, rcsdata)
       CALL write_data(oname, rcsdata)

    END IF

    DEALLOCATE(rcsdata)

    WRITE(*,*) sec_to_str(timer_end())

  END SUBROUTINE read_rcst

  SUBROUTINE read_plrz(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: n
    REAL (KIND=dp) :: omega, a
    COMPLEX (KIND=dp) :: ri
    COMPLEX (KIND=dp), DIMENSION(3) :: alpha
    REAL (KIND=dp), DIMENSION(b%nwl,7) :: data
    TYPE(prdnfo), POINTER :: prd

    READ(line,*) a

    WRITE(*,*) 'Computing polarizability'

    prd => NULL()

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(b,data,prd,a)&
    !$OMP PRIVATE(n,omega,ri,alpha)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,b%nwl
       ri = b%media( b%domains(1)%medium_index )%prop(n)%ri
       omega = 2*pi*c0/b%sols(n)%wl

       !IF(b%domains(1)%gf_index>0) THEN
       !   prd => b%prd(b%domains(1)%gf_index)
       !ELSE
       !   prd => NULL()
       !END IF

       alpha = polarizability(b%domains(1)%mesh, b%ga, b%sols(n)%x,&
            b%mesh%nedges, omega, ri, prd, a, b%qd_tri)

       data(n,1) = b%sols(n)%wl
       data(n,2) = REAL(alpha(1),KIND=dp)
       data(n,3) = IMAG(alpha(1))
       data(n,4) = REAL(alpha(2),KIND=dp)
       data(n,5) = IMAG(alpha(2))
       data(n,6) = REAL(alpha(3),KIND=dp)
       data(n,7) = IMAG(alpha(3))
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    CALL write_data(TRIM(b%name) // '-plrz.dat', data)

  END SUBROUTINE read_plrz

  SUBROUTINE read_crst(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: fid=10, iovar, n, ns, medium_index, domain_index
    REAL (KIND=dp) :: csca, cabs, wl
    REAL (KIND=dp), DIMENSION(b%nwl,3) :: data
    REAL (KIND=dp) :: omega, factor
    COMPLEX (KIND=dp) :: ri
    TYPE(srcdata) :: nlsrc
    CHARACTER (LEN=256) :: oname, numstr

    ! Can be set to -1 if boundary orientation is to be flipped.
    factor = 1.0_dp

    medium_index = 1

    domain_index = 1

    IF(LEN_TRIM(line)/=0) THEN
       READ(line,*) medium_index, domain_index, factor
    END IF

    WRITE(*,*) 'Computing extinction cross-section'

    ! Consider all excitation sources.
    DO ns=1,SIZE(b%src)

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(b,data,ns,medium_index,domain_index,factor)&
       !$OMP PRIVATE(n,omega,ri,wl,csca,cabs)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,b%nwl
          wl = b%sols(n)%wl
          
          omega = 2.0_dp*pi*c0/wl
          !ri = b%media(b%domains(1)%medium_index)%prop(n)%ri
          ri = b%media(medium_index)%prop(n)%ri
          
          CALL cs_prtsrf(b%domains(domain_index)%mesh, b%mesh%nedges, omega, ri, b%ga,&
               b%sols(n)%x(:,:,ns), b%src(ns), b%qd_tri, csca, cabs)
          
          data(n,1) = wl
          data(n,2) = csca*factor
          data(n,3) = cabs*factor
          
          WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       WRITE(numstr, '(A,I0)') '-s', ns
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.crs'
       
       CALL write_data(oname, data)
    END DO

    ! Second-harmonic.
    IF(ALLOCATED(b%sols(1)%nlx)) THEN
       nlsrc%type = 0

       ! Consider all excitation sources.
       DO ns=1,SIZE(b%src)
          
          !$OMP PARALLEL DEFAULT(NONE)&
          !$OMP SHARED(b,data,nlsrc,ns,medium_index,domain_index,factor)&
          !$OMP PRIVATE(n,omega,ri,wl,csca,cabs)
          !$OMP DO SCHEDULE(STATIC)
          DO n=1,b%nwl
             wl = b%sols(n)%wl
             
             omega = 2.0_dp*pi*c0/wl
             ri = b%media(medium_index)%prop(n)%shri
             
             CALL cs_prtsrf(b%domains(domain_index)%mesh, b%mesh%nedges, 2.0_dp*omega, ri, b%ga,&
                  b%sols(n)%nlx(:,:,ns), nlsrc, b%qd_tri, csca, cabs)
             
             data(n,1) = wl
             data(n,2) = csca*factor
             data(n,3) = cabs*factor
             
             WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
          END DO
          !$OMP END DO
          !$OMP END PARALLEL

          WRITE(numstr, '(A,I0)') '-s', ns
          oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-sh.crs'
          
          CALL write_data(oname, data)
       END DO
    END IF
  END SUBROUTINE read_crst

  SUBROUTINE read_csca(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: fid=10, iovar, n
    REAL (KIND=dp) :: csca, cabs, wl
    REAL (KIND=dp), DIMENSION(b%nwl,2) :: data
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    TYPE(srcdata) :: nlsrc

    WRITE(*,*) 'Computing scattering cross-section'

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(b,data)&
    !$OMP PRIVATE(n,omega,ri,wl,csca)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,b%nwl
       wl = b%sols(n)%wl

       omega = 2.0_dp*pi*c0/wl
       ri = b%media(b%domains(1)%medium_index)%prop(n)%ri

       CALL csca_ff(b%domains(1)%mesh, b%mesh%nedges, omega, ri, b%ga,&
            b%sols(n)%x(:,:,1), b%qd_tri, csca)

       data(n,1) = wl
       data(n,2) = csca

       WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    CALL write_data(TRIM(b%name) // '.crs', data)

    ! Second-harmonic.
    IF(ALLOCATED(b%sols(1)%nlx)) THEN
       nlsrc%type = 0

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(b,data)&
       !$OMP PRIVATE(n,omega,ri,wl,csca)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,b%nwl
          wl = b%sols(n)%wl
          
          omega = 2.0_dp*pi*c0/wl
          ri = b%media(b%domains(1)%medium_index)%prop(n)%shri
          
          CALL csca_ff(b%domains(1)%mesh, b%mesh%nedges, 2.0_dp*omega, ri, b%ga,&
               b%sols(n)%nlx(:,:,1), b%qd_tri, csca)
          
          data(n,1) = wl
          data(n,2) = csca
          
          WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
       
       CALL write_data(TRIM(b%name) // '-sh.crs', data)
    END IF
  END SUBROUTINE read_csca

  SUBROUTINE read_vext(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: fid=10, iovar, n, ns
    REAL (KIND=dp) :: cext, wl
    REAL (KIND=dp), DIMENSION(b%nwl,2) :: data
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri, k
    TYPE(srcdata) :: nlsrc
    CHARACTER (LEN=256) :: oname, numstr
    TYPE(prdnfo), POINTER :: prd

    WRITE(*,*) 'Computing extinction cross-section from VIE solution'

    prd => NULL()

    ! Consider all excitation sources.
    DO ns=1,SIZE(b%src)

!       !$OMP PARALLEL DEFAULT(NONE)&
!       !$OMP SHARED(b,data,ns,prd)&
!       !$OMP PRIVATE(n,omega,ri,k,wl,cext)
!       !$OMP DO SCHEDULE(STATIC)
       DO n=1,b%nwl
          wl = b%sols(n)%wl
          
          omega = 2.0_dp*pi*c0/wl
          ri = b%media(b%domains(1)%medium_index)%prop(n)%ri
          k = ri*omega/c0

          cext = cext_vie(b%mesh, k, b%ga(1), prd, b%qd_tetra, xi_hom,&
               b%sols(n)%x(:,1,ns), b%src(ns))
          
          data(n,1) = wl
          data(n,2) = cext
          
          WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
       END DO
!       !$OMP END DO
!       !$OMP END PARALLEL

       WRITE(numstr, '(A,I0)') '-s', ns
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.crs'
       
       CALL write_data(oname, data)
    END DO

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

  END SUBROUTINE read_vext

  SUBROUTINE test()
    COMPLEX (KIND=dp) :: res

    CALL asqz2(fun, 0.0_dp, pi, 0.0_dp, 3.0_dp, 1d-6, 10, res)

    WRITE(*,*) REAL(res)

  CONTAINS
    FUNCTION fun(x,y) RESULT(z)
      REAL (KIND=dp), INTENT(IN) :: x, y
      COMPLEX (KIND=dp) :: z

      z = EXP(-y)*(SIN(4*x)**2)

    END FUNCTION fun
  END SUBROUTINE test

  SUBROUTINE msgloop()
    CHARACTER (LEN=256) :: line
    CHARACTER (LEN=128) :: scmd
    CHARACTER (LEN=256) :: filename
    TYPE(batch) :: b
    INTEGER :: stat, n, m

    CALL batch_defaults(b)

    DO
       stat = 0
       WRITE(*,'(A2)', ADVANCE='NO') '> '
       READ(*,'(A)',IOSTAT=stat) line
       ! READ(*,'(A4)',IOSTAT=stat, ADVANCE='NO') scmd

       IF(LEN_TRIM(line)==0) THEN
          CYCLE
       END IF

       READ(line, *) scmd

       IF(LEN_TRIM(line)>LEN_TRIM(scmd)) THEN
          line = line((LEN_TRIM(scmd)+1):LEN_TRIM(line))
       ELSE
          line = ''
       END IF

       IF(stat<0 .OR. scmd=='exit') THEN
          EXIT
       ELSE IF(scmd=='name') THEN
          READ(line,*) b%name
       ELSE IF(scmd=='mesh') THEN
          CALL read_mesh(line, b)
       ELSE IF(scmd=='wlrg') THEN
          CALL read_wlrange(line, b)
       ELSE IF(scmd=='wlls') THEN
          CALL read_wllist(line, b)
       ELSE IF(scmd=='symm') THEN
          CALL read_symm(line, b)
       ELSE IF(scmd=='nmed') THEN
          CALL read_nmed(line, b)
       ELSE IF(scmd=='smed') THEN
          CALL read_smed(line, b)
       ELSE IF(scmd=='nsrc') THEN
          CALL read_nsrc(line, b)
       ELSE IF(scmd=='ssrc') THEN
          CALL read_ssrc(line, b)
       ELSE IF(scmd=='osrc') THEN
          CALL read_osrc(line, b)
       ELSE IF(scmd=='solv') THEN
          CALL solve_batch(b)
       ELSE IF(scmd=='npgf') THEN
          CALL read_npgf(line,b)
       ELSE IF(scmd=='ipgw') THEN
          CALL read_ipgw(line, b)
       ELSE IF(scmd=='ndom') THEN
          CALL read_ndom(line, b)
       ELSE IF(scmd=='sdom') THEN
          CALL read_sdom(line, b)
       ELSE IF(scmd=='fdom') THEN
          CALL read_fdom(line, b)
       ELSE IF(scmd=='sbnd') THEN
          CALL read_sbnd(line, b)
       ELSE IF(scmd=='nfms') THEN
          CALL read_nfms(line, b)
       ELSE IF(scmd=='nfdm') THEN
          CALL read_nfdm(line, b)
       ELSE IF(scmd=='nfem') THEN
          CALL read_nfem(line, b)
       ELSE IF(scmd=='nfst') THEN
          CALL read_nfst(line, b)
       ELSE IF(scmd=='nfpl') THEN
          CALL read_nfpl(line, b)
       ELSE IF(scmd=='rcst') THEN
          CALL read_rcst(line, b)
       ELSE IF(scmd=='rcs2') THEN
          CALL read_rcs2(line, b)
       ELSE IF(scmd=='rcs3') THEN
          CALL read_rcs3(line, b)
       ELSE IF(scmd=='crst') THEN
          CALL read_crst(line, b)
       ELSE IF(scmd=='diff') THEN
          CALL read_diff(line, b)
       ELSE IF(scmd=='plrz') THEN
          CALL read_plrz(line, b)
       ELSE IF(scmd=='scan') THEN
          CALL read_scan_source(line, b)
       ELSE IF(scmd=='trns') THEN
          CALL read_trns(line, b)
       ELSE IF(scmd=='quad') THEN
          CALL read_quad(line, b)
       ELSE IF(scmd=='csca') THEN
          CALL read_csca(line, b)
       ELSE IF(scmd=='test') THEN
          CALL read_test(line, b)
       ELSE IF(scmd=='vsol') THEN
          CALL solve_batch_vie(b)
       ELSE IF(scmd=='vext') THEN
          CALL read_vext(line, b)
       ELSE IF(scmd=='viem') THEN
          CALL vie_eigen_spec()
       ELSe IF(scmd=='msrc') THEN
          CALL read_move_source(line, b)
       ELSe IF(scmd=='fres') THEN
          CALL read_fres(line, b)
       ELSE
          WRITE(*,*) 'Unrecognized command ', scmd, '!'
       END IF
    END DO

    CALL delete_batch(b)

  END SUBROUTINE msgloop

END MODULE interface
