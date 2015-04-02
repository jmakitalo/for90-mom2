! MODULE: ffields
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for computing scattered fields far away from the scatterer.
! Contains also higher level functions to compute radar cross-sections.
MODULE ffields
  USE rwgf
  USE symmetry
  USE dipole

  IMPLICIT NONE

CONTAINS
  ! Calculates the far-fields of E and H fields in spherical coordinates in the observation
  ! direction (theta,phi) at distance r.
  SUBROUTINE far_fields(mesh, nedgestot, omega, ri, ga, x, r, theta, phi, qd, et, ep, ht, hp)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    REAL (KIND=dp), INTENT(IN) :: r, theta, phi
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp), INTENT(INOUT) :: et, ep, ht, hp

    INTEGER :: n, q, index, nweights, t
    COMPLEX (KIND=dp) :: k, Nt, Np, Lt, Lp, quad1, quad2, quad3, quad4, phasor, prefix
    REAL (KIND=dp) :: st, ct, sp, cp, An, detj
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn
    REAL (KIND=dp), DIMENSION(3) :: rp, jg, uphi, utheta
    COMPLEX (KIND=dp), DIMENSION(3) :: jn, mn, ed, hd
    COMPLEX (KIND=dp), DIMENSION(nedgestot) :: alpha, beta
    INTEGER :: nf, ns
    COMPLEX (KIND=dp) :: gae

    st = SIN(theta)
    ct = COS(theta)
    sp = SIN(phi)
    cp = COS(phi)

    k = ri*omega/c0

    nweights = qd%num_nodes

    Nt = 0.0_dp
    Np = 0.0_dp
    Lt = 0.0_dp
    Lp = 0.0_dp

    DO nf=1,SIZE(ga)
       alpha = x(1:nedgestot,nf)
       beta = x((nedgestot+1):(2*nedgestot),nf)

       DO n=1,mesh%nfaces
          An = mesh%faces(n)%area
          qpn = quad_tri_points(qd, n, mesh)

          DO q=1,3
             index = mesh%faces(n)%edge_indices(q)
             index = mesh%edges(index)%parent_index

             DO ns=1,SIZE(ga)

                quad1 = 0.0_dp
                quad2 = 0.0_dp
                quad3 = 0.0_dp
                quad4 = 0.0_dp

                detj = ga(ns)%detj
                gae = ga(ns)%ef(nf)

                DO t=1,nweights
                   rp = MATMUL(ga(ns)%j, qpn(:,t))
                   jg = MATMUL(ga(ns)%j, rwg(qpn(:,t),n,q,mesh))

                   jn = jg*alpha(index)*gae
                   mn = jg*beta(index)*gae*detj

                   phasor = EXP(-(0,1)*k*(rp(1)*st*cp + rp(2)*st*sp + rp(3)*ct))
                   
                   quad1 = quad1 + qd%weights(t)*(ct*cp*jn(1) + ct*sp*jn(2) - st*jn(3))*phasor
                   quad2 = quad2 + qd%weights(t)*(-sp*jn(1) + cp*jn(2))*phasor
                   
                   quad3 = quad3 + qd%weights(t)*(ct*cp*mn(1) + ct*sp*mn(2) - st*mn(3))*phasor
                   quad4 = quad4 + qd%weights(t)*(-sp*mn(1) + cp*mn(2))*phasor
                END DO
                
                Nt = Nt + quad1*An
                Np = Np + quad2*An
                
                Lt = Lt + quad3*An
                Lp = Lp + quad4*An
             END DO
          END DO
       END DO
       
    END DO

    prefix = (0,1)*k*EXP((0,1)*k*r)/(4.0_dp*pi*r)

    et = prefix*(Lp + eta0/ri*Nt)
    ep = -prefix*(Lt - eta0/ri*Np)

    ht = -prefix*(Np - Lt*ri/eta0)
    hp = prefix*(Nt + Lp*ri/eta0)

    ! If we have dipole source, superpose its contribution.
    !IF(b%src_type==src_dipole) THEN
    !   CALL dipole_farfields(r*get_dir(theta,phi), b%dpos, b%dmom, omega, ri, ed, hd)

    !   uphi = (/-SIN(phi), COS(phi), 0.0_dp/)
    !   utheta = (/COS(phi)*COS(theta), SIN(phi)*COS(theta), -SIN(theta)/)

    !   ep = ep + dotc(ed, CMPLX(uphi,KIND=dp))
    !   et = et + dotc(ed, CMPLX(utheta,KIND=dp))

    !   hp = hp + dotc(hd, CMPLX(uphi,KIND=dp))
    !   ht = ht + dotc(hd, CMPLX(utheta,KIND=dp))
    !END IF
  END SUBROUTINE far_fields

  ! Radar cross-section integrated over solid angle theta_max
  ! around the z-axis in +z-direction.
  ! This routine is deprecated as it uses the fixed-order Simpson quadrature.
  SUBROUTINE rcs_solangle_old(mesh, nedgestot, omega, ri, ga, x, theta_max, qd, scatp)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    REAL (KIND=dp), INTENT(IN) :: theta_max
    TYPE(quad_data), INTENT(IN) :: qd

    REAL (KIND=dp), INTENT(OUT) :: scatp
    INTEGER, PARAMETER :: nquad = 11
    INTEGER :: n, m
    REAL (KIND=dp), DIMENSION(1:nquad) :: qwt, ptt, qwp, ptp
    REAL (KIND=dp) :: p
    COMPLEX (KIND=dp) :: et, ep, ht, hp

    ! Compute weights and nodes from Simpson's rule.
    CALL get_simpsons_weights(0.0_dp, theta_max, nquad-1, qwt)
    CALL get_simpsons_points(0.0_dp, theta_max, nquad-1, ptt)
    CALL get_simpsons_weights(0.0_dp, 2.0_dp*pi, nquad-1, qwp)
    CALL get_simpsons_points(0.0_dp, 2.0_dp*pi, nquad-1, ptp)

    scatp = 0.0_dp

    DO n=1,nquad
       DO m=1,nquad

          CALL far_fields(mesh, nedgestot, omega, ri, ga, x, 1.0_dp,&
               ptt(n), ptp(m), qd, et, ep, ht, hp)

          p = ABS(et)**2 + ABS(ep)**2

          scatp = scatp + qwt(n)*qwp(m)*p*SIN(ptt(n))
       END DO
    END DO

  END SUBROUTINE rcs_solangle_old

  ! Radar cross-section integrated over solid angle range [theta1,theta2].
  ! around the z-axis in +z-direction.
  SUBROUTINE rcs_solangle(mesh, nedgestot, omega, ri, ga, x, theta1, theta2, qd, scatp)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    REAL (KIND=dp), INTENT(IN) :: theta1, theta2
    TYPE(quad_data), INTENT(IN) :: qd

    REAL (KIND=dp), INTENT(OUT) :: scatp
    INTEGER, PARAMETER :: maxDepth = 10
    REAL (KIND=dp), PARAMETER :: maxerr = 1D-3
    COMPLEX (KIND=dp) :: scatpz

    CALL asqz2(integ, theta1, theta2, 0.0_dp, 2*PI, maxerr, maxDepth, scatpz)

    scatp = REAL(scatpz)

  CONTAINS
    FUNCTION integ(theta, phi) RESULT(res)
      REAL (KIND=dp), INTENT(IN) :: theta, phi
      COMPLEX (KIND=dp) :: res
      COMPLEX (KIND=dp) :: et, ep, ht, hp
      
      CALL far_fields(mesh, nedgestot, omega, ri, ga, x, 1.0_dp,&
           theta, phi, qd, et, ep, ht, hp)
      
      res = SIN(theta)*(ABS(et)**2 + ABS(ep)**2)
    END FUNCTION integ

  END SUBROUTINE rcs_solangle

  ! Bistatic radar-cross section a.k.a. radiated power per unit solid angle.
  ! fdir specifies a polarization filter direction vector if filter is true.
  SUBROUTINE rcs_ext(mesh, nedgestot, omega, ri, ga, x, ntheta, nphi, filter, fdir, qd, rcsdata)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ntheta, nphi
    LOGICAL, INTENT(IN) :: filter
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: fdir
    TYPE(quad_data), INTENT(IN) :: qd

    REAL (KIND=dp), DIMENSION(1:ntheta,1:nphi), INTENT(OUT) :: rcsdata

    INTEGER :: n, m
    COMPLEX (KIND=dp) :: et, ep, ht, hp
    REAL (KIND=dp) :: r=1, theta, phi, c
    REAL (KIND=dp), DIMENSION(3) :: vtheta, vphi
    COMPLEX (KIND=dp), DIMENSION(3) :: e

    !c = 4.0_dp*pi*(r**2)
    c = 1.0_dp

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(ntheta,nphi,r,rcsdata,c,mesh,nedgestot,omega,ri,ga,x,filter,fdir,qd)&
    !$OMP PRIVATE(n,m,theta,phi,et,ep,ht,hp,vtheta,vphi,e)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,ntheta
       DO m=1,nphi
          theta = pi*REAL((n-1),KIND=dp)/REAL((ntheta-1),KIND=dp)
          phi = 2.0_dp*pi*REAL((m-1),KIND=dp)/REAL((nphi-1),KIND=dp)

          CALL far_fields(mesh, nedgestot, omega, ri, ga, x, r, theta, phi, qd, et, ep, ht, hp)

          IF(filter==.FALSE.) THEN
             rcsdata(n,m) = c*(ABS(et)**2 + ABS(ep)**2)
          ELSE
             vtheta = (/COS(phi)*COS(theta), SIN(phi)*COS(theta), -SIN(theta)/)
             vphi = (/-SIN(phi), COS(phi), 0.0_dp/)
             e = et*vtheta + ep*vphi
             rcsdata(n,m) = c*ABS(dotc(CMPLX(fdir,KIND=dp), e))**2
          END IF
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  END SUBROUTINE rcs_ext

  ! Bistatic radar-cross section a.k.a. radiated power per unit solid angle.
  SUBROUTINE rcs(mesh, nedgestot, omega, ri, ga, x, ntheta, nphi, qd, rcsdata)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ntheta, nphi
    TYPE(quad_data), INTENT(IN) :: qd

    REAL (KIND=dp), DIMENSION(1:ntheta,1:nphi), INTENT(OUT) :: rcsdata

    CALL rcs_ext(mesh, nedgestot, omega, ri, ga, x, ntheta, nphi, .FALSE.,&
         (/0.0_dp,0.0_dp,0.0_dp/), qd, rcsdata)
  END SUBROUTINE rcs

END MODULE ffields
