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

  IMPLICIT NONE

CONTAINS
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
    CHARACTER (LEN=256) :: type, method, file
    COMPLEX (KIND=dp) :: val

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

       IF(method=='value') THEN
          READ(line,*) mindex, type, method, b%media(mindex)%prop(1)%nlb%chi2zzz

          DO n=1,b%nwl
             b%media(mindex)%prop(n)%nlb%chi2zzz = b%media(mindex)%prop(1)%nlb%chi2zzz
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
    CHARACTER (LEN=3) :: ext

    READ(line,*) b%mesh_file, b%scale

    ext = getext(b%mesh_file)
    IF(ext=='msh') THEN
       b%mesh = load_mesh(b%mesh_file)
       CALL build_mesh(b%mesh, b%scale)
    ELSE
       WRITE(*,*) 'Invalid mesh extension!'
       STOP
    END IF
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

  SUBROUTINE read_source(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    CHARACTER (LEN=32) :: token
    CHARACTER (LEN=256) :: solname, src_meshname, nfocus
    CHARACTER (LEN=3) :: ext
    INTEGER :: nga, nfrags

    ALLOCATE(b%src(1))

    b%src(1)%pos(:) = 0.0_dp

    READ(line,*) token
    IF(token=='pw') THEN
       b%src(1)%type = src_pw
       READ(line,*) token, b%src(1)%theta, b%src(1)%phi, b%src(1)%psi
       b%src(1)%theta = b%src(1)%theta*degtorad
       b%src(1)%phi = b%src(1)%phi*degtorad
       b%src(1)%psi = b%src(1)%psi*degtorad
    ELSE IF(token=='focus_rad') THEN
       b%src(1)%type = src_focus_rad
       READ(line,*) token, b%src(1)%focal, b%src(1)%waist, b%src(1)%napr, nfocus
    ELSE IF(token=='focus_x') THEN
       b%src(1)%type = src_focus_x
       READ(line,*) token, b%src(1)%focal, b%src(1)%waist, b%src(1)%napr, nfocus
    ELSE IF(token=='focus_y') THEN
       b%src(1)%type = src_focus_y
       READ(line,*) token, b%src(1)%focal, b%src(1)%waist, b%src(1)%napr, nfocus
    ELSE IF(token=='focus_hg01') THEN
       b%src(1)%type = src_focus_hg01
       READ(line,*) token, b%src(1)%focal, b%src(1)%waist, b%src(1)%napr, nfocus
    ELSE IF(token=='focus_azimut') THEN
       b%src(1)%type = src_focus_azim
       READ(line,*) token, b%src(1)%focal, b%src(1)%waist, b%src(1)%napr, nfocus
    ELSE IF(token=='dipole') THEN
       b%src(1)%type = src_dipole
       READ(line,*) token, b%src(1)%pos, b%src(1)%dmom
    ELSE
       WRITE(*,*) 'Invalid source type!'
    END IF

    IF(b%src(1)%type==src_focus_rad .OR. b%src(1)%type==src_focus_x .OR.&
         b%src(1)%type==src_focus_y .OR. b%src(1)%type==src_focus_hg01 .OR.&
         b%src(1)%type==src_focus_azim) THEN
       IF(nfocus=='true') THEN
          b%src(1)%nfocus = .TRUE.
       ELSE IF(nfocus=='false') THEN
          b%src(1)%nfocus = .FALSE.
       ELSE
          WRITE(*,*) 'Unrecognized source argument value!'
          STOP
       END IF
    END IF
  END SUBROUTINE read_source

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

  SUBROUTINE read_diff(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: fid=10, iovar, n, dindex
    REAL (KIND=dp) :: wl, omega
    REAL (KIND=dp), DIMENSION(b%nwl) :: irr
    COMPLEX (KIND=dp) :: ri, ri_inc
    TYPE(prdnfo), POINTER :: prd

    READ(line,*) dindex

    WRITE(*,*) 'Computing diffracted irradiance'

    DO n=1,b%nwl
       wl = b%sols(n)%wl

       omega = 2.0_dp*pi*c0/wl
       ri = b%media(b%domains(dindex)%medium_index)%prop(n)%ri
       ri_inc = b%media(b%domains(1)%medium_index)%prop(n)%ri

       prd => b%prd(b%domains(dindex)%gf_index)

       prd%cwl = find_closest(wl, prd%coef(:)%wl)

       irr(n) = diff_irradiance(b%domains(dindex)%mesh, b%ga, dindex==1, b%src(1), b%sols(n)%x(:,:,1),&
            b%mesh%nedges, omega, ri, ri_inc, prd)

       WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
    END DO

    OPEN(fid, FILE=(TRIM(b%name) // '.dif'), ACTION='WRITE', IOSTAT=iovar)
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
          
          irr(n) = diff_irradiance(b%domains(dindex)%mesh, b%ga, .FALSE., b%src(1), b%sols(n)%nlx(:,:,1),&
               b%mesh%nedges, omega, ri, ri_inc, prd)
          
          WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
       END DO

       OPEN(fid, FILE=(TRIM(b%name) // '-sh.dif'), ACTION='WRITE', IOSTAT=iovar)
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

  SUBROUTINE read_nfms(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    TYPE(nfield_plane) :: nfplane
    INTEGER :: wlindex, dindex
    CHARACTER (LEN=256) :: oname, numstr
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri

    READ(line,*) wlindex, dindex

    WRITE(numstr, '(A,I0,A,I0)') '-wl', wlindex, '-d', dindex
    oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '.msh'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%ri
       
    CALL field_mesh(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
         b%sols(wlindex)%x(:,:,1), b%ga, omega, ri)

    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       oname = TRIM(b%name) // TRIM(ADJUSTL(numstr)) // '-sh.msh'
       
       ri = b%media(b%domains(dindex)%medium_index)%prop(wlindex)%shri
       
       CALL field_mesh(oname, b%domains(dindex)%mesh, b%scale, b%mesh%nedges,&
            b%sols(wlindex)%nlx(:,:,1), b%ga, 2.0_dp*omega, ri)

    END IF
  END SUBROUTINE read_nfms

  ! Computes a 2D image, where pixels correspond to source positions.
  ! Data may be plotted direcly with MATLAB's scimage, whence x-axis
  ! points right and y-axis points up.
  SUBROUTINE read_rcs2(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: iovar, npt, wlindex, n
    REAL (KIND=dp) :: omega, theta_max
    COMPLEX (KIND=dp) :: ri
    REAL (KIND=dp), DIMENSION(SIZE(b%src)) :: scatp
    CHARACTER (LEN=256) :: oname, numstr

    READ(line,*) wlindex

    WRITE(*,*) 'Computing scattered power'

    CALL timer_start()

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(1)%medium_index)%prop(wlindex)%ri

    npt = NINT(SQRT(REAL(SIZE(b%src))))

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(b,omega,ri,scatp,wlindex)&
    !$OMP PRIVATE(n,theta_max)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,SIZE(b%src)
       theta_max = ASIN(b%src(n)%napr)

       CALL rcs_solangle(b%domains(1)%mesh, b%mesh%nedges, omega, ri, b%ga,&
            b%sols(wlindex)%x(:,:,n), theta_max, scatp(n))
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
       !$OMP SHARED(b,omega,ri,scatp,wlindex)&
       !$OMP PRIVATE(n,theta_max)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,SIZE(b%src)
          theta_max = ASIN(b%src(n)%napr)

          CALL rcs_solangle(b%domains(1)%mesh, b%mesh%nedges, 2.0_dp*omega, ri, b%ga,&
               b%sols(wlindex)%nlx(:,:,n), theta_max, scatp(n))
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       CALL write_data(oname, RESHAPE(scatp,(/npt,npt/)))

    END IF

    WRITE(*,*) sec_to_str(timer_end())

  END SUBROUTINE read_rcs2

  SUBROUTINE read_rcst(line, b)
    CHARACTER (LEN=*), INTENT(IN) :: line
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: iovar, ntheta_rcs, nphi_rcs, wlindex
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: rcsdata
    CHARACTER (LEN=256) :: oname, numstr

    READ(line,*) wlindex, ntheta_rcs, nphi_rcs

    WRITE(*,*) 'Computing radar cross-sections'
    CALL timer_start()

    ALLOCATE(rcsdata(1:ntheta_rcs,1:nphi_rcs))

    WRITE(numstr, '(I0)') wlindex
    oname = TRIM(b%name) // '-wl' // TRIM(ADJUSTL(numstr)) // '.rcs'

    omega = 2.0_dp*pi*c0/b%sols(wlindex)%wl
    ri = b%media(b%domains(1)%medium_index)%prop(wlindex)%ri

    CALL rcs(b%domains(1)%mesh, b%mesh%nedges, omega, ri, b%ga, b%sols(wlindex)%x(:,:,1),&
         ntheta_rcs, nphi_rcs, rcsdata)
    CALL write_data(oname, rcsdata)

    ! Nonlinear radar cross-sections.
    IF(ALLOCATED(b%sols(wlindex)%nlx)) THEN
       oname = TRIM(b%name) // '-wl' // TRIM(ADJUSTL(numstr)) // '-sh.rcs'
       
       ! Second-harmonic frequency.
       ri = b%media(b%domains(1)%medium_index)%prop(wlindex)%shri
       
       CALL rcs(b%domains(1)%mesh, b%mesh%nedges, 2.0_dp*omega, ri, b%ga, b%sols(wlindex)%nlx(:,:,1),&
            ntheta_rcs, nphi_rcs, rcsdata)
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

       alpha = polarizability(b%domains(1)%mesh, b%ga, b%sols(n)%x(:,:,1),&
            b%mesh%nedges, omega, ri, prd, a)

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
    INTEGER :: fid=10, iovar, n
    REAL (KIND=dp) :: csca, cabs, wl
    REAL (KIND=dp), DIMENSION(b%nwl,3) :: data
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    TYPE(srcdata) :: nlsrc

    WRITE(*,*) 'Computing extinction cross-section'

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(b,data)&
    !$OMP PRIVATE(n,omega,ri,wl,csca,cabs)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,b%nwl
       wl = b%sols(n)%wl

       omega = 2.0_dp*pi*c0/wl
       ri = b%media(b%domains(1)%medium_index)%prop(n)%ri

       CALL cs_prtsrf(b%domains(1)%mesh, b%mesh%nedges, omega, ri, b%ga,&
            b%sols(n)%x(:,:,1), b%src(1), csca, cabs)

       data(n,1) = wl
       data(n,2) = csca
       data(n,3) = cabs

       WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    CALL write_data(TRIM(b%name) // '.crs', data)

    ! Second-harmonic.
    IF(ALLOCATED(b%sols(1)%nlx)) THEN
       nlsrc%type = 0

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(b,data,nlsrc)&
       !$OMP PRIVATE(n,omega,ri,wl,csca,cabs)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,b%nwl
          wl = b%sols(n)%wl
          
          omega = 2.0_dp*pi*c0/wl
          ri = b%media(b%domains(1)%medium_index)%prop(n)%shri
          
          CALL cs_prtsrf(b%domains(1)%mesh, b%mesh%nedges, 2.0_dp*omega, ri, b%ga,&
               b%sols(n)%nlx(:,:,1), nlsrc, csca, cabs)
          
          data(n,1) = wl
          data(n,2) = csca
          data(n,3) = cabs
          
          WRITE(*,'(A,I0,A,I0,:)') ' Wavelength ',  n, ' of ', b%nwl
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
       
       CALL write_data(TRIM(b%name) // '-sh.crs', data)
    END IF
  END SUBROUTINE read_crst

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
       ELSE IF(scmd=='srce') THEN
          CALL read_source(line, b)
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
       ELSE IF(scmd=='rcst') THEN
          CALL read_rcst(line, b)
       ELSE IF(scmd=='rcs2') THEN
          CALL read_rcs2(line, b)
       ELSE IF(scmd=='crst') THEN
          CALL read_crst(line, b)
       ELSE IF(scmd=='diff') THEN
          CALL read_diff(line, b)
       ELSE IF(scmd=='plrz') THEN
          CALL read_plrz(line, b)
       ELSE IF(scmd=='scan') THEN
          CALL read_scan_source(line, b)
       ELSE
          WRITE(*,*) 'Unrecognized command ', scmd, '!'
       END IF
    END DO

    CALL delete_batch(b)

  END SUBROUTINE msgloop

END MODULE interface
