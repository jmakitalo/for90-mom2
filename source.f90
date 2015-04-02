! MODULE: source
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for evaluating various source excitations, such as a plane-wave
! and focused Gaussian beams of various forms. Contains also routines for
! evaluating source vectors of the form <f_m,s>, where f_m is an RWG
! function and s is a source excitation.
MODULE source
  USE bessel
  USE aux
  USE mesh
  USE symmetry
  USE quad
  USE rwgf
  USE dipole

  IMPLICIT NONE

  INTEGER, PARAMETER :: src_pw = 1,&
       src_focus_rad = 2,&
       src_focus_x = 3,&
       src_focus_y = 4,&
       src_nlsurf = 5,&
       src_focus_hg01 = 6,&
       src_focus_azim = 7,&
       src_dipole = 8,&
       src_nlbulk = 9,&
       src_pw_elliptic = 10

  INTEGER, PARAMETER :: focustype_radial = 1,&
       focustype_x = 2,&
       focustype_y = 3,&
       focustype_hg01 = 4,&
       focustype_azimut = 5

  TYPE srcdata
     ! src_type: Identifier for the excitation source.
     INTEGER :: type

     ! Angles for a plane-wave excitation.
     REAL (KIND=dp) :: theta, phi, psi

     ! Phase factor between Ex and Ey, where (Ex,Ey,k) form a right-handed system.
     COMPLEX (KIND=dp) :: phase

     ! Focused beam excitation parameters.
     REAL (KIND=dp) :: focal, waist, napr

     ! Normalize beam fields to maximum at focal plane?
     LOGICAL :: nfocus

     ! Excitation source position for, e.g., beams and a dipole.
     REAL (KIND=dp), DIMENSION(3) :: pos

     ! Dipole moment of a dipole source.
     COMPLEX (KIND=dp), DIMENSION(3) :: dmom

  END type srcdata

CONTAINS
  ! Electric field of plane-wave with polarization pv,
  ! wave-vector dir*k evaluated at point r.
  FUNCTION pw_efield(pv, dir, k, r) RESULT(efield)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pv, dir, r
    COMPLEX (KIND=dp), INTENT(IN) :: k
    COMPLEX (KIND=dp), DIMENSION(3) :: efield

    efield = pv*EXP((0,1)*k*dotr(dir, r))
  END FUNCTION pw_efield

  ! Magnetic field of plane-wave with polarization pv,
  ! wave-vector dir*k evaluated at point r.
  FUNCTION pw_hfield(pv, dir, k, omega, r) RESULT(hfield)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pv, dir, r
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), DIMENSION(3) :: hfield

    hfield = k/(omega*mu0)*crossc(CMPLX(dir,KIND=dp), pw_efield(pv, dir, k, r))
  END FUNCTION pw_hfield

  ! Electric and magnetic field of a plane-wave with k-vector and
  ! polarization determined by the angles theta, phi and psi.
  SUBROUTINE pw_fields(theta, phi, psi, omega, ri, pt, ef, hf)
    REAL (KIND=dp), INTENT(IN) :: theta, phi, psi, omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(INOUT) :: ef, hf

    REAL (KIND=dp), DIMENSION(3) :: pol, dir
    COMPLEX (KIND=dp) :: k

    dir = get_dir(theta, phi)
    pol = get_pol(theta, phi, psi)

    k = ri*omega/c0

    ef = pw_efield(pol, dir, k, pt)
    hf = pw_hfield(pol, dir, k, omega, pt)
  END SUBROUTINE pw_fields

  SUBROUTINE pw_fields_elliptic(theta, phi, phase, omega, ri, pt, ef, hf)
    REAL (KIND=dp), INTENT(IN) :: theta, phi, omega
    COMPLEX (KIND=dp), INTENT(IN) :: phase, ri
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(INOUT) :: ef, hf

    COMPLEX (KIND=dp), DIMENSION(3) :: efx, efy
    REAL (KIND=dp), DIMENSION(3) :: pol, dir
    COMPLEX (KIND=dp) :: k

    dir = get_dir(theta, phi)
    pol = get_pol(theta, phi, 0.0_dp)

    k = ri*omega/c0

    efx = pw_efield(pol, dir, k, pt)
    efy = crossc(CMPLX(dir,KIND=dp), efx)
    ef = efx + phase*efy
    ef = ef/normc(ef)

    hf = k/(omega*mu0)*crossc(CMPLX(dir,KIND=dp), ef)
  END SUBROUTINE pw_fields_elliptic

  ! Source vector <f_m,s> with f_m RWG function and s a specified
  ! source excitation function.
  ! mesh: boundary mesh of the domain where the source is defined in.
  ! nedgestot: total number of edges in the whole system.
  ! omega: angular frequency
  ! ri: refractive index in the source domain.
  ! ga: group actions.
  ! src: source parameters
  ! q: the source vector for all group representations
  SUBROUTINE srcvec(mesh, nedgestot, omega, ri, ga, src, qd, q)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(srcdata), INTENT(IN) :: src
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: q

    INTEGER :: nweights, m, p, r, index, focustype, nr
    COMPLEX (KIND=dp), DIMENSION(3,SIZE(ga),qd%num_nodes,mesh%nfaces) :: ef, hf
    COMPLEX (KIND=dp), DIMENSION(3) :: emax
    COMPLEX (KIND=dp), DIMENSION(SIZE(ga)) :: Ie, Ih
    REAL (KIND=dp) :: Am, pm
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpm
    REAL (KIND=dp), DIMENSION(3) :: fm, pos
    LOGICAL :: print_info

    print_info = .FALSE.

    nweights = qd%num_nodes

    q(:,:) = 0.0_dp

    ! Do we have a dipole source? Requires special treatment.
    IF(src%type==src_dipole) THEN
       CALL srcvec_dipole(mesh, nedgestot, omega, ri, ga, src%pos, src%dmom, qd, q)
       RETURN
    END IF

    ! If we have a focused beam, determine its type and compute
    ! maximum electric field value.
    IF(print_info .AND. src%type/=src_pw .AND. src%type/=src_pw_elliptic) THEN
       IF(src%type==src_focus_rad) THEN
          focustype = focustype_radial
       ELSE IF(src%type==src_focus_x) THEN
          focustype = focustype_x
       ELSE IF(src%type==src_focus_y) THEN
          focustype = focustype_y
       ELSE IF(src%type==src_focus_hg01) THEN
          focustype = focustype_hg01
       ELSE IF(src%type==src_focus_azim) THEN
          focustype = focustype_azimut
       END IF

       emax = compute_focus_maximum(src%focal, src%napr, src%waist, ri, omega, focustype)

       WRITE(*,*) '|max(E)| at focus: ', normc(emax)

       IF(focustype==focustype_hg01) THEN
          pm = beam_max_pt(omega/c0, src%focal, src%waist, ASIN(src%napr/REAL(ri,KIND=dp)), focustype)
          WRITE(*,*) 'Beam max. location: ', (/0.0_dp,pm,0.0_dp/)
       ELSE IF(focustype==focustype_azimut) THEN
          pm = beam_max_pt(omega/c0, src%focal, src%waist, ASIN(src%napr/REAL(ri,KIND=dp)), focustype)
          WRITE(*,*) 'Beam max. location: ', (/pm,0.0_dp,0.0_dp/)
       END IF
    END IF

    ! Pre-compute excitation fields in parallel.

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(mesh,nweights,qd,src,omega,ri,ga,ef,hf)&
    !$OMP PRIVATE(m,qpm,r)
    !$OMP DO SCHEDULE(STATIC)
    DO m=1,mesh%nfaces
       qpm = quad_tri_points(qd, m, mesh)

       DO r=1,nweights
          ! Compute the source excitation for all representations.
          CALL src_field_frags(src, omega, ri, ga, qpm(:,r), ef(:,:,r,m), hf(:,:,r,m))
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    ! Go through faces of the mesh.
    DO m=1,mesh%nfaces
       qpm = quad_tri_points(qd, m, mesh)
       Am = mesh%faces(m)%area
    
       ! Go through edges of the current face.
       DO p=1,3
          Ie = 0.0_dp
          Ih = 0.0_dp
        
          ! Go through integration points on the face.
          DO r=1,nweights
             fm = rwg(qpm(:,r), m, p, mesh)

             ! Compute inner-product integrands for each representation.
             DO nr=1,SIZE(ga)
                Ie(nr) = Ie(nr) + qd%weights(r)*dotc(CMPLX(fm,KIND=dp), ef(:,nr,r,m))
                Ih(nr) = Ih(nr) + qd%weights(r)*dotc(CMPLX(fm,KIND=dp), hf(:,nr,r,m))
             END DO
          END DO
        
          ! Get local index.
          index = mesh%faces(m)%edge_indices(p)

          ! Map to global index.
          index = mesh%edges(index)%parent_index

          DO nr=1,SIZE(ga)
             q(index,nr) = q(index,nr) + Ie(nr)*Am
             
             q(index+nedgestot,nr) = q(index+nedgestot,nr) + Ih(nr)*Am
          END DO
       END DO
    END DO
  END SUBROUTINE srcvec

  ! Computes source excitation fields in all given group representations.
  SUBROUTINE src_field_frags(src, omega, ri, ga, ptin, einc, hinc)
    TYPE(srcdata), INTENT(IN) :: src
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: ptin

    COMPLEX (KIND=dp), DIMENSION(3,SIZE(ga)), INTENT(INOUT) :: einc, hinc
    COMPLEX (KIND=dp), DIMENSION(3) :: ef, hf
    REAL (KIND=dp), DIMENSION(3) :: pt
    INTEGER :: nr, na, focustype
    COMPLEX (KIND=dp) :: gae
    REAL (KIND=dp) :: detj

    einc(:,:) = 0.0_dp
    hinc(:,:) = 0.0_dp

    IF(src%type==src_focus_rad) THEN
       focustype = focustype_radial
    ELSE IF(src%type==src_focus_x) THEN
       focustype = focustype_x
    ELSE IF(src%type==src_focus_y) THEN
       focustype = focustype_y
    ELSE IF(src%type==src_focus_hg01) THEN
       focustype = focustype_hg01
    ELSE IF(src%type==src_focus_azim) THEN
       focustype = focustype_azimut
    END IF

    ! Loop through group actions.
    DO na=1,SIZE(ga)
       
       ! Evaluate the source at p_g^-1(ptin).
       IF(src%type==src_pw) THEN
          pt = MATMUL(TRANSPOSE(ga(na)%j), ptin)
          
          CALL pw_fields(src%theta, src%phi, src%psi, omega, ri, pt, ef, hf)
       ELSE IF(src%type==src_pw_elliptic) THEN
          pt = MATMUL(TRANSPOSE(ga(na)%j), ptin)
          
          CALL pw_fields_elliptic(src%theta, src%phi, src%phase, omega, ri, pt, ef, hf)
       ELSE IF(src%type==src_focus_rad .OR. src%type==src_focus_azim .OR.&
            src%type==src_focus_x .OR. src%type==src_focus_y .OR.&
            src%type==src_focus_hg01) THEN
          pt = MATMUL(TRANSPOSE(ga(na)%j), ptin) - src%pos
          
          CALL compute_focus(src%focal, src%napr, src%waist, ri, omega, pt,&
               ef, hf, focustype, src%nfocus)
       END IF

       ! Loop through group representations.
       DO nr=1,SIZE(ga)

          detj = ga(na)%detj
          gae = ga(na)%ef(nr)

          ! Compute the projections of the fields in the group representations.
          einc(:,nr) = einc(:,nr) + MATMUL(ga(na)%j, ef)*gae
          hinc(:,nr) = hinc(:,nr) + MATMUL(ga(na)%j, hf)*gae*detj
       END DO
    END DO

    ! Divide by the number of group elements.
    einc(:,:) = einc(:,:)/SIZE(ga)
    hinc(:,:) = hinc(:,:)/SIZE(ga)
  END SUBROUTINE src_field_frags

  ! Computes the source at point pt. This is a simplified version of
  ! src_field_frags and does not consider symmetry representations.
  SUBROUTINE src_fields(src, omega, ri, pt, einc, hinc)
    TYPE(srcdata), INTENT(IN) :: src
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(INOUT) :: einc, hinc
    INTEGER :: focustype

    einc(:) = 0.0_dp
    hinc(:) = 0.0_dp

    IF(src%type==src_focus_rad) THEN
       focustype = focustype_radial
    ELSE IF(src%type==src_focus_x) THEN
       focustype = focustype_x
    ELSE IF(src%type==src_focus_y) THEN
       focustype = focustype_y
    ELSE IF(src%type==src_focus_hg01) THEN
       focustype = focustype_hg01
    ELSE IF(src%type==src_focus_azim) THEN
       focustype = focustype_azimut
    END IF

    IF(src%type==src_pw) THEN
       CALL pw_fields(src%theta, src%phi, src%psi, omega, ri, pt, einc, hinc)
    ELSE IF(src%type==src_pw_elliptic) THEN
       CALL pw_fields_elliptic(src%theta, src%phi, src%phase, omega, ri, pt, einc, hinc)
    ELSE IF(src%type==src_focus_rad .OR. src%type==src_focus_azim .OR.&
         src%type==src_focus_x .OR. src%type==src_focus_y .OR.&
         src%type==src_focus_hg01) THEN
       CALL compute_focus(src%focal, src%napr, src%waist, ri, omega, pt - src%pos,&
            einc, hinc, focustype, src%nfocus)
    END IF

  END SUBROUTINE src_fields

  FUNCTION compute_focus_maximum(f, na, w0, ri, omega, focustype) RESULT(emax)
    REAL (KIND=dp), INTENT(IN) :: f, na, w0, omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    INTEGER, INTENT(IN) :: focustype

    COMPLEX (KIND=dp), DIMENSION(3) :: emax, hmax
    REAL (KIND=dp) :: pm

    IF(focustype==focustype_hg01) THEN
       pm = beam_max_pt(omega/c0, f, w0, ASIN(na/REAL(ri,KIND=dp)), focustype)
       CALL compute_focus_NH(f, na, w0, ri, omega, (/0.0_dp,pm,0.0_dp/), emax, hmax, focustype)
    ELSE IF(focustype==focustype_azimut) THEN
       pm = beam_max_pt(omega/c0, f, w0, ASIN(na/REAL(ri,KIND=dp)), focustype)
       CALL compute_focus_NH(f, na, w0, ri, omega, (/pm,0.0_dp,0.0_dp/), emax, hmax, focustype)
    ELSE IF(focustype==focustype_y) THEN
       CALL compute_focus_NH(f, na, w0, ri, omega, (/0.0_dp,0.0_dp,0.0_dp/), emax, hmax, focustype_x)
    ELSE
       CALL compute_focus_NH(f, na, w0, ri, omega, (/0.0_dp,0.0_dp,0.0_dp/), emax, hmax, focustype)
    END IF
  END FUNCTION compute_focus_maximum

  ! The expressions for the fields are based on the cartesian coordinate system (x',y',z') in
  ! which the beam propagates in the positive z'-direction, as presented by Novotny.
  ! We wish the beam to propagate in the (-z)-direction in out cartesian coordinate
  ! system (x,y,z). For this reason we apply a contravariant coordinate transformation for
  ! the results obtained in the frame (x',y',z').
  SUBROUTINE compute_focus(f, na, w0, ri, omega, posLocal, e, h, focustype, nfocus)
    REAL (KIND=dp), INTENT(IN) :: f, na, w0, omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: posLocal
    INTEGER, INTENT(IN) :: focustype
    LOGICAL, INTENT(IN) :: nfocus

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(OUT) :: e, h
    COMPLEX (KIND=dp), DIMENSION(3) :: e0, h0, emax, hmax
    REAL (KIND=dp) :: scale, pm
    REAL (KIND=dp), DIMENSION(3) :: pos

    ! Map the posLocal from (x,y,z) to pos in (x',y',z')
    pos = (/posLocal(1), -posLocal(2), -posLocal(3)/)

    IF(focustype==focustype_y) THEN
       CALL compute_focus_NH(f, na, w0, ri, omega, (/pos(2), -pos(1), pos(3)/), e, h, focustype_x)
       e = (/-e(2), e(1), e(3)/)
       h = (/-h(2), h(1), h(3)/)
    ELSE
       CALL compute_focus_NH(f, na, w0, ri, omega, pos, e, h, focustype)
    END IF

    ! Normalize fields to maximum value at focal plane?
    IF(nfocus) THEN
       emax = compute_focus_maximum(f, na, w0, ri, omega, focustype)
       scale = 1.0_dp/normc(emax)
       
       e = e*scale
       h = h*scale
    END IF

    ! The obtained (e,h) is now in the frame (x',y',z'), so we transform it to (x,y,z).
    ! e and h transform similarly although h is a pseudovector. This is because the rotation
    ! is a proper transformation.
    e = (/e(1), -e(2), -e(3)/)
    h = (/h(1), -h(2), -h(3)/)
  END SUBROUTINE compute_focus

  ! Always call the subroutine compute_focus to evaluate the fields.
  ! Returns fields of a beam propagating along +z direction.
  ! The fields are normalized to incident plane-wave amplitude E0.
  ! Based on expressions by Novotny and Hecht.
  SUBROUTINE compute_focus_NH(f, na, w0, ri, omega, pos, e, h, focustype)
    REAL (KIND=dp), INTENT(IN) :: f, na, w0, omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pos
    INTEGER, INTENT(IN) :: focustype

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(OUT) :: e, h

    REAL (KIND=dp) :: r, z, phi, theta_max, k, maxerr
    COMPLEX (KIND=dp) :: I10, I11, I12, I00, I01, I02, I13, I14, prefix
    INTEGER :: maxDepth

    maxerr = 1D-4
    maxDepth = 10

    r = SQRT(pos(1)**2 + pos(2)**2)
    phi = ATAN2(pos(2), pos(1))
    z = pos(3)

    k = REAL(ri,KIND=dp)*omega/c0

    theta_max = ASIN(na/REAL(ri,KIND=dp))

    IF(focustype==focustype_radial) THEN
       prefix = (1.0_dp/SQRT(ri))*(0,1)*k*(f**2)*EXP(-(0,1)*k*f)/(2*w0)

       I10 = asqz(integI10, 0.0_dp, theta_max, maxerr, maxDepth)
       I11 = asqz(integI11, 0.0_dp, theta_max, maxerr, maxDepth)
       I12 = asqz(integI12, 0.0_dp, theta_max, maxerr, maxDepth)

       !I10 = focus_int(r,z,k,f,w0,theta_max,0)
       !I11 = focus_int(r,z,k,f,w0,theta_max,1)
       !I12 = focus_int(r,z,k,f,w0,theta_max,2)

       e = prefix*(/(0,1)*(I11 - I12)*COS(phi), (0,1)*(I11 - I12)*SIN(phi), -4.0_dp*I10/)
       h = ri*prefix/eta0*(/-(0,1)*(I11 + 3.0_dp*I12)*&
            SIN(phi), (0,1)*(I11 + 3.0_dp*I12)*COS(phi), (0.0_dp,0.0_dp)/)
    ELSE IF(focustype==focustype_x) THEN
       prefix = (1.0_dp/SQRT(ri))*(0,1)*k*f*EXP(-(0,1)*k*f)/2.0_dp

       I00 = asqz(integI00, 0.0_dp, theta_max, maxerr, maxDepth)
       I01 = asqz(integI01, 0.0_dp, theta_max, maxerr, maxDepth)
       I02 = asqz(integI02, 0.0_dp, theta_max, maxerr, maxDepth)

       !I00 = focus_int(r,z,k,f,w0,theta_max,3)
       !I01 = focus_int(r,z,k,f,w0,theta_max,4)
       !I02 = focus_int(r,z,k,f,w0,theta_max,5)

       e = prefix*(/I00 + I02*COS(2.0_dp*phi), I02*SIN(2.0_dp*phi), -2.0_dp*(0,1)*I01*COS(phi)/)
       h = ri*prefix/eta0*(/I02*SIN(2.0_dp*phi), I00 - I02*COS(2.0_dp*phi), -2.0_dp*(0,1)*I01*&
            SIN(phi)/)
    ELSE IF(focustype==focustype_y) THEN
       WRITE(*,*) 'Focused y-polarized beam not yet implemented.'
       STOP
    ELSE IF(focustype==focustype_hg01) THEN
       prefix = (1.0_dp/SQRT(ri))*(0,1)*k*(f**2)*EXP(-(0,1)*k*f)/(2.0_dp*w0)

       I10 = asqz(integI10, 0.0_dp, theta_max, maxerr, maxDepth)
       I11 = asqz(integI11, 0.0_dp, theta_max, maxerr, maxDepth)
       I12 = asqz(integI12, 0.0_dp, theta_max, maxerr, maxDepth)
       I13 = asqz(integI13, 0.0_dp, theta_max, maxerr, maxDepth)
       I14 = asqz(integI14, 0.0_dp, theta_max, maxerr, maxDepth)

       !I10 = focus_int(r,z,k,f,w0,theta_max,0)
       !I11 = focus_int(r,z,k,f,w0,theta_max,1)
       !I12 = focus_int(r,z,k,f,w0,theta_max,2)
       !I13 = focus_int(r,z,k,f,w0,theta_max,6)
       !I14 = focus_int(r,z,k,f,w0,theta_max,7)

       e = prefix*(/(0,1)*(I11 + 2.0_dp*I12)*SIN(phi) + (0,1)*I14*SIN(3.0_dp*phi),&
            -(0,1)*I12*COS(phi) - (0,1)*I14*COS(3.0_dp*phi), 2.0_dp*I13*SIN(2.0_dp*phi)/)
       h = ri*prefix/eta0*(/-(0,1)*I12*COS(phi) - (0,1)*I14*COS(3.0_dp*phi),&
            (0,1)*I11*SIN(phi) - (0,1)*I14*SIN(3.0_dp*phi), -2.0_dp*I10 - 2.0_dp*I13*&
            COS(2.0_dp*phi)/)
    ELSE IF(focustype==focustype_azimut) THEN
       prefix = (1.0_dp/SQRT(ri))*(0,1)*k*(f**2)*EXP(-(0,1)*k*f)/(2*w0)

       I10 = asqz(integI10, 0.0_dp, theta_max, maxerr, maxDepth)
       I11 = asqz(integI11, 0.0_dp, theta_max, maxerr, maxDepth)
       I12 = asqz(integI12, 0.0_dp, theta_max, maxerr, maxDepth)

       !I10 = focus_int(r,z,k,f,w0,theta_max,0)
       !I11 = focus_int(r,z,k,f,w0,theta_max,1)
       !I12 = focus_int(r,z,k,f,w0,theta_max,2)

       e = prefix*(/(0,1)*(I11 + 3.0_dp*I12)*SIN(phi), -(0,1)*(I11 + 3.0_dp*I12)*&
            COS(phi), (0.0_dp,0.0_dp)/)
       h = ri*prefix/eta0*(/(0,1)*(I11 - I12)*COS(phi), (0,1)*(I11 - I12)*SIN(phi), -4.0_dp*I10/)
    END IF

    ! Define integrands of variable theta.
  CONTAINS
    FUNCTION integCommon(theta) RESULT(common)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common

      common = EXP(-(f**2)*(SIN(theta)**2)/(w0**2))*SQRT(COS(theta))*EXP((0,1)*k*z*COS(theta))
    END FUNCTION integCommon

    FUNCTION integI10(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j0

      common = integCommon(theta)
      CALL besselj(besm,0,maxerr,k*r*SIN(theta),j0)
      integ = common*(SIN(theta)**3)*j0
    END FUNCTION integI10

    FUNCTION integI11(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j1

      common = integCommon(theta)
      CALL besselj(besm,1,maxerr,k*r*SIN(theta),j1)
      integ = common*(SIN(theta)**2)*(1.0_dp + 3.0_dp*COS(theta))*j1
    END FUNCTION integI11

    FUNCTION integI12(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j1

      common = integCommon(theta)
      CALL besselj(besm,1,maxerr,k*r*SIN(theta),j1)
      integ = common*(SIN(theta)**2)*(1.0_dp - COS(theta))*j1
    END FUNCTION integI12

    FUNCTION integI00(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j0

      common = integCommon(theta)
      CALL besselj(besm,0,maxerr,k*r*SIN(theta),j0)
      integ = common*SIN(theta)*(1.0_dp + COS(theta))*j0
    END FUNCTION integI00

    FUNCTION integI01(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j1

      common = integCommon(theta)
      CALL besselj(besm,1,maxerr,k*r*SIN(theta),j1)
      integ = common*(SIN(theta)**2)*j1
    END FUNCTION integI01

    FUNCTION integI02(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j2

      common = integCommon(theta)
      CALL besselj(besm,2,maxerr,k*r*SIN(theta),j2)
      integ = common*SIN(theta)*(1.0_dp - COS(theta))*j2
    END FUNCTION integI02

    FUNCTION integI13(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j2

      common = integCommon(theta)
      CALL besselj(besm,2,maxerr,k*r*SIN(theta),j2)
      integ = common*SQRT(CMPLX(COS(theta)))*(SIN(THETA)**3)*j2
    END FUNCTION integI13

    FUNCTION integI14(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j3

      common = integCommon(theta)
      CALL besselj(besm,3,maxerr,k*r*SIN(theta),j3)
      integ = common*SQRT(CMPLX(COS(theta)))*(SIN(THETA)**2)*(1.0_dp-COS(theta))*j3
    END FUNCTION integI14
  END SUBROUTINE compute_focus_NH

  SUBROUTINE test_focus(focustype)
    INTEGER, INTENT(IN) :: focustype
    INTEGER :: n, m, k, npoints
    INTEGER, DIMENSION(6) :: fids
    REAL (KIND=dp) :: width, height, tn, tm, f, na, w0, omega, wl
    REAL (KIND=dp), DIMENSION(3) :: pt
    COMPLEX (KIND=dp), DIMENSION(3) :: e, h
    COMPLEX (KIND=dp) :: ri

    npoints = 100
    width = 1500e-9
    height = 1500e-9

    wl = 1060.0_dp*1D-9
    omega = 2.0_dp*pi*c0/wl
    f = 1.0_dp*1D-3
    w0 = 0.8_dp*1D-3
    na = 0.8_dp
    ri = 1.0_dp

    fids = (/10,11,12,13,14,15/)

    OPEN(fids(1), FILE='focus-ex.amp', ACTION='WRITE')
    OPEN(fids(2), FILE='focus-ex.pha', ACTION='WRITE')
    OPEN(fids(3), FILE='focus-ey.amp', ACTION='WRITE')
    OPEN(fids(4), FILE='focus-ey.pha', ACTION='WRITE')
    OPEN(fids(5), FILE='focus-ez.amp', ACTION='WRITE')
    OPEN(fids(6), FILE='focus-ez.pha', ACTION='WRITE')

    DO n=1,npoints
       DO m=1,npoints
          tn = REAL(n-1,KIND=dp)/REAL(npoints-1,KIND=dp)
          tm = REAL(m-1,KIND=dp)/REAL(npoints-1,KIND=dp)
          pt = (/-width/2+width*tn, 0.0_dp, -height/2+height*tm/)

          CALL compute_focus(f, na, w0, ri, omega, pt, e, h, focustype, .FALSE.)

          WRITE(fids(1),'(EN15.3)', ADVANCE='NO') ABS(e(1))
          WRITE(fids(2),'(EN15.3)', ADVANCE='NO') ATAN2(IMAG(e(1)), REAL(e(1)))
          WRITE(fids(3),'(EN15.3)', ADVANCE='NO') ABS(e(2))
          WRITE(fids(4),'(EN15.3)', ADVANCE='NO') ATAN2(IMAG(e(2)), REAL(e(2)))
          WRITE(fids(5),'(EN15.3)', ADVANCE='NO') ABS(e(3))
          WRITE(fids(6),'(EN15.3)', ADVANCE='NO') ATAN2(IMAG(e(3)), REAL(e(3)))
       END DO

       DO k=1,6
          WRITE(fids(k), '(/)')
       END DO

    END DO

    DO k=1,6
       CLOSE(fids(k))
    END DO

    !pt(:) = 0.0_dp
    !pt(2) = beam_max_pt(2*pi/wl, f, w0, ASIN(na/REAL(ri,KIND=dp)), focustype)
    !CALL compute_focus_denorm(f, na, w0, ri, omega, pt, e, h, focustype)
    !WRITE(*,*) pt(2), normc(e)
  END SUBROUTINE test_focus

  FUNCTION besselj1d(a,x) RESULT(y)
    REAL (KIND=dp), INTENT(IN) :: a, x
    REAL (KIND=dp) :: y, j0, j1, maxerr = 1e-6
    INTEGER :: besm

    CALL besselj(besm,0,maxerr,a*x,j0)
    CALL besselj(besm,1,maxerr,a*x,j1)
    y = a*j0 - 1.0_dp/x*j1
  END FUNCTION besselj1d

  FUNCTION besselj1dd(a,x) RESULT(y)
    REAL (KIND=dp), INTENT(IN) :: a, x
    REAL (KIND=dp) :: y, j0, j1, maxerr = 1e-6
    INTEGER :: besm

    CALL besselj(besm,0,maxerr,a*x,j0)
    CALL besselj(besm,1,maxerr,a*x,j1)
    y = -(a**2)*j1 - a/x*j0 + 1.0_dp/(x**2)*j1
  END FUNCTION besselj1dd

  FUNCTION besselj3d(a,x) RESULT(y)
    REAL (KIND=dp), INTENT(IN) :: a, x
    REAL (KIND=dp) :: y, j2, j3, maxerr = 1e-6
    INTEGER :: besm

    CALL besselj(besm,2,maxerr,a*x,j2)
    CALL besselj(besm,3,maxerr,a*x,j3)
    y = a*j2 - 3.0_dp/x*j3
  END FUNCTION besselj3d

  FUNCTION besselj3dd(a,x) RESULT(y)
    REAL (KIND=dp), INTENT(IN) :: a, x
    REAL (KIND=dp) :: y, j1, j2, j3, maxerr = 1e-6
    INTEGER :: besm

    CALL besselj(besm,1,maxerr,a*x,j1)
    CALL besselj(besm,2,maxerr,a*x,j2)
    CALL besselj(besm,3,maxerr,a*x,j3)
    y = (a**2)*j1 - 2.0_dp*a/x*j2 - 3.0_dp*a/x*j2 + 9.0_dp/(x**2)*j3
  END FUNCTION besselj3dd

  FUNCTION beam_max_pt(k, f, w0, theta_max, focustype) RESULT(r)
    INTEGER, INTENT(IN) :: focustype
    REAL (KIND=dp), INTENT(IN) :: k, f, w0, theta_max

    REAL (KIND=dp) :: r, z, fd, fdd, err, r_prev, maxerr
    INTEGER :: n, maxDepth
    COMPLEX (KIND=dp) :: I11, I12, I11d, I12d, I11dd, I12dd, I14, I14d, I14dd, c, cd, cdd

    maxerr = 1D-6
    maxDepth = 10

    z = 0.0_dp

    ! Initial guess.
    IF(focustype==focustype_azimut) THEN
       r = 0.5_dp*3.8317_dp/(k*SIN(theta_max))
    ELSE IF(focustype==focustype_hg01) THEN
       r = 0.5_dp*3.8317_dp/(k*SIN(theta_max))
    END IF

    DO n=1,20
       I11 = asqz(integI11, 0.0_dp, theta_max, maxerr, maxDepth)
       I12 = asqz(integI12, 0.0_dp, theta_max, maxerr, maxDepth)
       I11d = asqz(integDrhoI11, 0.0_dp, theta_max, maxerr, maxDepth)
       I12d = asqz(integDrhoI12, 0.0_dp, theta_max, maxerr, maxDepth)
       I11dd = asqz(integDDrhoI11, 0.0_dp, theta_max, maxerr, maxDepth)
       I12dd = asqz(integDDrhoI12, 0.0_dp, theta_max, maxerr, maxDepth)

       IF(focustype==focustype_azimut) THEN
          !I11 = focus_int(rho,0.0_dp,k,f,w0,theta_max,1)
          !I12 = focus_int(rho,0.0_dp,k,f,w0,theta_max,2)
          !I11d = focus_int(rho,0.0_dp,k,f,w0,theta_max,8)
          !I12d = focus_int(rho,0.0_dp,k,f,w0,theta_max,9)
          !I11dd = focus_int(rho,0.0_dp,k,f,w0,theta_max,10)
          !I12dd = focus_int(rho,0.0_dp,k,f,w0,theta_max,11)

          fd = REAL((I11d + 3*I12d)*CONJG(I11 + 3*I12) + (I11 + 3*I12)*CONJG(I11d + 3*I12d))
          fdd = REAL((I11dd + 3*I12dd)*CONJG(I11 + 3*I12) + (I11d + 3*I12d)*CONJG(I11d + 3*I12d) +&
               (I11d + 3*I12d)*CONJG(I11d + 3*I12d) + (I11 + 3*I12)*CONJG(I11dd + 3*I12dd))
       ELSE IF(focustype==focustype_hg01) THEN
          I14 = asqz(integI14, 0.0_dp, theta_max, maxerr, maxDepth)
          I14d = asqz(integDrhoI14, 0.0_dp, theta_max, maxerr, maxDepth)
          I14dd = asqz(integDDrhoI14, 0.0_dp, theta_max, maxerr, maxDepth)

          !I14 = focus_int(rho,0.0_dp,k,f,w0,theta_max,7)
          !I14d = focus_int(rho,0.0_dp,k,f,w0,theta_max,12)
          !I14dd = focus_int(rho,0.0_dp,k,f,w0,theta_max,13)

          c = I11 + 2.0_dp*I12 - I14
          cd = I11d + 2.0_dp*I12d - I14d
          cdd = I11dd + 2.0_dp*I12dd - I14dd
          fd = REAL(cd*CONJG(c) + c*CONJG(cd))
          fdd = REAL(cdd*CONJG(c) + 2.0_dp*cd*CONJG(cd) + c*CONJG(cdd))
       END IF

       r_prev = r
       r = r - fd/fdd

       err = ABS(r_prev - r)
       IF((err < 1e-4) .AND. n>4) THEN
          RETURN
       END IF
       
    END DO

    WRITE(*,*) 'Could not obtain error bounds for source normalization!'
    WRITE(*,*) 'Obtained error was ', err

  CONTAINS
    FUNCTION integCommon(theta) RESULT(common)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common

      common = EXP(-(f**2)*(SIN(theta)**2)/(w0**2))*SQRT(COS(theta))*EXP((0,1)*k*z*COS(theta))
    END FUNCTION integCommon

    FUNCTION integI11(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j1

      common = integCommon(theta)
      CALL besselj(besm,1,maxerr,k*r*SIN(theta),j1)
      integ = common*(SIN(theta)**2)*(1.0_dp + 3.0_dp*COS(theta))*j1
    END FUNCTION integI11

    FUNCTION integI12(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j1

      common = integCommon(theta)
      CALL besselj(besm,1,maxerr,k*r*SIN(theta),j1)
      integ = common*(SIN(theta)**2)*(1.0_dp - COS(theta))*j1
    END FUNCTION integI12

    FUNCTION integI14(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      REAL (KIND=dp) :: j3

      common = integCommon(theta)
      CALL besselj(besm,3,maxerr,k*r*SIN(theta),j3)
      integ = common*SQRT(CMPLX(COS(theta)))*(SIN(THETA)**2)*(1.0_dp-COS(theta))*j3
    END FUNCTION integI14

    FUNCTION integDrhoI11(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm

      common = integCommon(theta)
      integ = common*(SIN(theta)**2)*(1.0_dp + 3.0_dp*COS(theta))*besselj1d(k*SIN(theta),r)
    END FUNCTION integDrhoI11

    FUNCTION integDrhoI12(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm

      common = integCommon(theta)
      integ = common*(SIN(theta)**2)*(1.0_dp - COS(theta))*besselj1d(k*SIN(theta),r)
    END FUNCTION integDrhoI12

    FUNCTION integDDrhoI11(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm

      common = integCommon(theta)
      integ = common*(SIN(theta)**2)*(1.0_dp + 3.0_dp*COS(theta))*besselj1dd(k*SIN(theta),r)
    END FUNCTION integDDrhoI11

    FUNCTION integDDrhoI12(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm

      common = integCommon(theta)
      integ = common*(SIN(theta)**2)*(1.0_dp - COS(theta))*besselj1dd(k*SIN(theta),r)
    END FUNCTION integDDrhoI12

    FUNCTION integDrhoI14(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      COMPLEX (KIND=dp) :: jd

      common = integCommon(theta)
      jd = besselj3d(k*SIN(theta),r)
      integ = common*SQRT(CMPLX(COS(theta)))*(SIN(THETA)**2)*(1.0_dp-COS(theta))*jd
    END FUNCTION integDrhoI14

    FUNCTION integDDrhoI14(theta) RESULT(integ)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: common, integ
      INTEGER :: besm
      COMPLEX (KIND=dp) :: jd

      common = integCommon(theta)
      jd = besselj3dd(k*SIN(theta),r)
      integ = common*SQRT(CMPLX(COS(theta)))*(SIN(THETA)**2)*(1.0_dp-COS(theta))*jd
    END FUNCTION integDDrhoI14

  END FUNCTION beam_max_pt

  SUBROUTINE print_source_info(src)
    TYPE(srcdata), INTENT(IN) :: src

    IF(src%type==src_pw) THEN
       WRITE(*,*) 'source type: ', 'pw'
       WRITE(*,'(A,T9,F6.2,A)') ' theta: ', src%theta*180/pi, ' deg'
       WRITE(*,'(A,T7,F6.2,A)') ' phi: ', src%phi*180/pi, ' deg'
       WRITE(*,'(A,T7,F6.2,A)') ' psi: ', src%psi*180/pi, ' deg'
    ELSE IF(src%type==src_pw_elliptic) THEN
       WRITE(*,*) 'source type: ', 'pw_elliptic'
       WRITE(*,'(A,T9,F6.2,A)') ' theta: ', src%theta*180/pi, ' deg'
       WRITE(*,'(A,T7,F6.2,A)') ' phi: ', src%phi*180/pi, ' deg'
       WRITE(*,'(A,T7,F6.2,A,F6.2,A)') ' phase: (', REAL(src%phase), ', ', AIMAG(src%phase), ')'
    ELSE IF(src%type==src_focus_rad) THEN
       WRITE(*,*) 'source type: ', 'focus_rad'
       WRITE(*,'(A,T16,E10.3)') ' focal length: ', src%focal
       WRITE(*,'(A,T14,E10.3)') ' beam waist: ', src%waist
       WRITE(*,'(A,T22,F6.2)') ' numerical aperture: ', src%napr
       IF(src%nfocus) THEN
          WRITE(*,*) 'normalized to maximum at focal plane'
       ELSE
          WRITE(*,*) 'normalized to incident plane-wave amplitude'
       END IF
    ELSE IF(src%type==src_focus_x) THEN
       WRITE(*,*) 'source type: ', 'focus_x'
       WRITE(*,'(A,T16,E10.3)') ' focal length: ', src%focal
       WRITE(*,'(A,T14,E10.3)') ' beam waist: ', src%waist
       WRITE(*,'(A,T22,F6.2)') ' numerical aperture: ', src%napr
       IF(src%nfocus) THEN
          WRITE(*,*) 'normalized to maximum at focal plane'
       ELSE
          WRITE(*,*) 'normalized to incident plane-wave amplitude'
       END IF
    ELSE IF(src%type==src_focus_y) THEN
       WRITE(*,*) 'source type: ', 'focus_y'
       WRITE(*,'(A,T16,E10.3)') ' focal length: ', src%focal
       WRITE(*,'(A,T14,E10.3)') ' beam waist: ', src%waist
       WRITE(*,'(A,T22,F6.2)') ' numerical aperture: ', src%napr
       IF(src%nfocus) THEN
          WRITE(*,*) 'normalized to maximum at focal plane'
       ELSE
          WRITE(*,*) 'normalized to incident plane-wave amplitude'
       END IF
    ELSE IF(src%type==src_focus_azim) THEN
       WRITE(*,*) 'source type: ', 'focus_azim'
       WRITE(*,'(A,T16,E10.3)') ' focal length: ', src%focal
       WRITE(*,'(A,T14,E10.3)') ' beam waist: ', src%waist
       WRITE(*,'(A,T22,F6.2)') ' numerical aperture: ', src%napr
       IF(src%nfocus) THEN
          WRITE(*,*) 'normalized to maximum at focal plane'
       ELSE
          WRITE(*,*) 'normalized to incident plane-wave amplitude'
       END IF
    ELSE IF(src%type==src_focus_hg01) THEN
       WRITE(*,*) 'source type: ', 'focus_hg01'
       WRITE(*,'(A,T16,E10.3)') ' focal length: ', src%focal
       WRITE(*,'(A,T14,E10.3)') ' beam waist: ', src%waist
       WRITE(*,'(A,T22,F6.2)') ' numerical aperture: ', src%napr
       IF(src%nfocus) THEN
          WRITE(*,*) 'normalized to maximum at focal plane'
       ELSE
          WRITE(*,*) 'normalized to incident plane-wave amplitude'
       END IF
    ELSE IF(src%type==src_dipole) THEN
       WRITE(*,*) 'source type: ', 'dipole'
    END IF
  END SUBROUTINE print_source_info
END MODULE source
