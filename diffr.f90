! MODULE: diffr
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for computing diffracted power in periodic problems.
MODULE diffr
  USE source
  USE nfields
  USE common

  IMPLICIT NONE

  INTEGER, PARAMETER :: max_prdsrc = 4
  INTEGER, PARAMETER :: fresnel_nquad = 25

CONTAINS
  FUNCTION Gpff(r, rp, k, prd, i, j) RESULT(g)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r, rp
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), INTENT(IN) :: k

    ! Diffraction orders.
    INTEGER, INTENT(IN) :: i, j

    COMPLEX (KIND=dp) :: g

    REAL (KIND=dp), DIMENSION(2) :: kt
    COMPLEX (KIND=dp) :: kz, phasor
    REAL (KIND=dp) :: sgn, A

    IF(r(3)>rp(3)) THEN
       sgn = 1.0_dp
    ELSE
       sgn = -1.0_dp
    END IF

    g = 0.0_dp
          
    ! Lattice vector.
    kt = (/prd%coef(prd%cwl)%k0x + 2.0_dp*PI*(i/(prd%dx*prd%cp)&
         - j*prd%sp/(prd%dy*prd%cp)) ,&
         prd%coef(prd%cwl)%k0y + 2.0_dp*PI*j/prd%dy/)
    
    ! Skip evanescent waves.
    IF(REAL(k**2,KIND=dp)<dotr(kt,kt)) THEN
       RETURN
    END IF
    
    kz = SQRT(k**2 - dotr(kt,kt))
    
    phasor = EXP((0,1)*dotr(kt,r(1:2)))*EXP(-(0,1)*dotr(kt,rp(1:2)))*&
         EXP(sgn*(0,1)*kz*r(3))*EXP(-sgn*(0,1)*kz*rp(3))
    
    g = phasor/kz

    A = prd%dx*prd%dy*prd%cp

    g = g*(0,1)/(2*A)
  END FUNCTION Gpff

  FUNCTION gradGpff(r, rp, k, prd, i, j) RESULT(gg)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r, rp
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), INTENT(IN) :: k

    ! Diffraction orders.
    INTEGER, INTENT(IN) :: i, j

    COMPLEX (KIND=dp), DIMENSION(3) :: gg

    REAL (KIND=dp), DIMENSION(2) :: kt
    COMPLEX (KIND=dp) :: kz, phasor
    REAL (KIND=dp) :: sgn, A

    IF(r(3)>rp(3)) THEN
       sgn = 1.0_dp
    ELSE
       sgn = -1.0_dp
    END IF

    gg(:) = 0.0_dp
          
    ! Lattice vector.
    kt = (/prd%coef(prd%cwl)%k0x + 2.0_dp*PI*(i/(prd%dx*prd%cp)&
         - j*prd%sp/(prd%dy*prd%cp)) ,&
         prd%coef(prd%cwl)%k0y + 2.0_dp*PI*j/prd%dy/)
    
    ! Skip evanescent waves.
    IF(REAL(k**2,KIND=dp)<dotr(kt,kt)) THEN
       RETURN
    END IF
    
    kz = SQRT(k**2 - dotr(kt,kt))
    
    phasor = EXP((0,1)*dotr(kt,r(1:2)))*EXP(-(0,1)*dotr(kt,rp(1:2)))*&
         EXP(sgn*(0,1)*kz*r(3))*EXP(-sgn*(0,1)*kz*rp(3))
    
    gg(1:2) = kt*phasor/kz
    gg(3) = sgn*phasor

    A = prd%dx*prd%dy*prd%cp

    gg(:) = gg(:)/(2*A)
  END FUNCTION gradGpff

  SUBROUTINE diff_fields(mesh, ga, nf, x, nedgestot, omega, ri, prd, r, i, j, e, h)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nf, nedgestot, i, j
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: x
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3), INTENT(INOUT) :: e, h

    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: qpn
    INTEGER :: n, q, t, edgeind
    COMPLEX (KIND=dp) :: c1, c2, g, k
    COMPLEX (KIND=dp), DIMENSION(3) :: gg
    REAL (KIND=dp), DIMENSION(3) :: divfn
    REAL (KIND=dp), DIMENSION(3,SIZE(qw),3) :: fv
    REAL (KIND=dp) :: An

    k = ri*omega/c0    
    c1 = (0,1)*omega*mu0
    c2 = (0,1)*omega*(ri**2)*eps0

    e(:) = 0.0_dp
    h(:) = 0.0_dp

    DO n=1,mesh%nfaces
       An = mesh%faces(n)%area
       qpn = GLquad_points(n, mesh)

       DO q=1,3
          CALL vrwg(qpn(:,:),n,q,mesh,fv(:,:,q))
          divfn(q) = rwgDiv(n,q,mesh)
       END DO

       DO t=1,SIZE(qw)
          g = Gpff(r, qpn(:,t), k, prd, i, j)
          gg = gradGpff(r, qpn(:,t), k, prd, i, j)

          DO q=1,3
             edgeind = mesh%faces(n)%edge_indices(q)
             edgeind = mesh%edges(edgeind)%parent_index

             e = e + qw(t)*An*( c1*g*fv(:,t,q)*x(edgeind) + gg*divfn(q)*x(edgeind)/c2 +&
                  crossc(gg, CMPLX(fv(:,t,q),KIND=dp))*x(edgeind + nedgestot) )

             h = h + qw(t)*An*( c2*g*fv(:,t,q)*x(edgeind + nedgestot) +&
                  gg*divfn(q)*x(edgeind + nedgestot)/c1 -&
                  crossc(gg, CMPLX(fv(:,t,q),KIND=dp))*x(edgeind) )
          END DO
       END DO
    END DO
  END SUBROUTINE diff_fields

  FUNCTION diff_irradiance(mesh, ga, addsrc, src, x, nedgestot, omega, ri, ri_inc, prd, i, j)&
       RESULT(irr)
    TYPE(mesh_container), INTENT(IN) :: mesh
    LOGICAL, INTENT(IN) :: addsrc
    TYPE(srcdata), INTENT(IN) :: src
    COMPLEX (KIND=dp), INTENT(IN) :: ri, ri_inc
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot, i, j
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x

    REAL (KIND=dp) :: irr, pinc, eval_dist, k
    INTEGER :: nf
    REAL (KIND=dp), DIMENSION(3) :: dir
    COMPLEX (KIND=dp), DIMENSION(3) :: e, h, einc, hinc
    REAL (KIND=dp), DIMENSION(2) :: kt

    nf = 1

    ! Field evaluation distance. Arbitrary positive value.
    ! For good numerical accuracy, should be on the order of
    ! wavelength.
    eval_dist = 1e-6

    !dir = get_dir(pwtheta, pwphi)

    k = REAL(ri,KIND=dp)*omega/c0

    kt = (/prd%coef(prd%cwl)%k0x + 2.0_dp*PI*(i/(prd%dx*prd%cp)&
         - j*prd%sp/(prd%dy*prd%cp)) ,&
         prd%coef(prd%cwl)%k0y + 2.0_dp*PI*j/prd%dy/)

    ! Skip evanescent waves.
    IF(REAL(k**2,KIND=dp)<dotr(kt,kt)) THEN
       irr = 0.0_dp
       RETURN
    END IF

    dir = (/kt(1), kt(2), -SQRT(k**2 - dotr(kt,kt))/)

    dir = dir/normr(dir)

    CALL diff_fields(mesh, ga, nf, x(:,nf), nedgestot, omega, ri, prd, dir*eval_dist, i, j, e, h)

    IF(addsrc .AND. i==0 .AND. j==0) THEN
       CALL src_fields(src, omega, ri, dir*eval_dist, einc, hinc)
       
       e = e + einc
       h = h + hinc
    END IF

    pinc = REAL(ri_inc,KIND=dp)/(c0*mu0)
    
    ! The relative irradiance diffracted to 0th order in the given domain.
    irr = dotr(REAL(crossc(e, CONJG(h)), KIND=dp), dir)/pinc

  END FUNCTION diff_irradiance

  FUNCTION transmittance(mesh, ga, addsrc, src, x, nedgestot, omega, ri, ri_inc, prd,&
       z0, zsign) RESULT(power)
    TYPE(mesh_container), INTENT(IN) :: mesh
    LOGICAL, INTENT(IN) :: addsrc
    TYPE(srcdata), INTENT(IN) :: src
    COMPLEX (KIND=dp), INTENT(IN) :: ri, ri_inc
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: x
    REAL (KIND=dp), INTENT(IN) :: z0, zsign

    REAL (KIND=dp) :: power
    REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: qwx, ptx, qwy, pty
    COMPLEX (KIND=dp), DIMENSION(3,1) :: e, h
    COMPLEX (KIND=dp), DIMENSION(3) :: einc, hinc, poynting
    COMPLEX (KIND=dp) :: k
    REAL (KIND=dp), DIMENSION(3) :: pt
    REAL (KIND=dp) :: hdx, hdy, pinc
    INTEGER :: n, m, nx, ny
    REAL (KIND=dp), DIMENSION(3) :: xaxis, yaxis

    ! Wavenumber in diffraction medium.
    k = ri*omega/c0

    ! Select the number of integration points based on wavelength and period.
    !nx = NINT(prd%dx/b%sols(wlindex)%wl*20)
    !ny = NINT(prd%dy/b%sols(wlindex)%wl*20)
    nx = 51
    ny = 51

    ! Make sure that the numbers are odd.
    IF(MOD(nx,2)==0) THEN
       nx = nx + 1
    END IF

    IF(MOD(ny,2)==0) THEN
       ny = ny + 1
    END IF

    ALLOCATE(qwx(1:nx), ptx(1:nx), qwy(1:ny), pty(1:ny))

    hdx = prd%dx*0.5_dp
    hdy = prd%dy*0.5_dp

    xaxis = (/prd%cp, prd%sp, 0.0_dp/)
    yaxis = (/0.0_dp, 1.0_dp, 0.0_dp/)

    ! Compute weights and nodes from Simpson's rule.
    CALL get_simpsons_weights(-hdx, hdx, nx-1, qwx)
    CALL get_simpsons_points(-hdx, hdx, nx-1, ptx)
    CALL get_simpsons_weights(-hdy, hdy, ny-1, qwy)
    CALL get_simpsons_points(-hdy, hdy, ny-1, pty)

    power = 0.0_dp

             
    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(ny,nx,z0,xaxis,yaxis,ptx,pty,mesh,ga,x,nedgestot,omega,ri,prd,addsrc,src,qwx,qwy,zsign,power)&
    !$OMP PRIVATE(m,n,pt,einc,hinc,e,h,poynting)
    !$OMP DO REDUCTION(+:power) SCHEDULE(STATIC)
    DO m=1,ny
       DO n=1,nx

          pt = (/0.0_dp,0.0_dp,z0/) + xaxis*ptx(n) + yaxis*pty(m)
          
          CALL scat_fields(mesh, ga, x, nedgestot, omega, ri, prd, pt, e, h)

          IF(addsrc) THEN
             CALL src_fields(src, omega, ri, pt, einc, hinc)

             e(:,1) = e(:,1) + einc
             h(:,1) = h(:,1) + hinc
          END IF

          poynting = crossc(e(:,1), CONJG(h(:,1)))

          power = power + 0.5_dp*qwx(n)*qwy(m)*REAL(poynting(3)*zsign)
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
        
    ! cp is the Jacobian of the area integration.
    power = power*prd%cp

    pinc = 0.5_dp*prd%dx*prd%dy*prd%cp*REAL(ri_inc,KIND=dp)/(c0*mu0)
    
    ! Relative power.
    power = power/pinc

    DEALLOCATE(qwx, ptx, qwy, pty)

  END FUNCTION transmittance
END MODULE diffr
