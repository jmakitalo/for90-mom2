! MODULE: strat
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implementation of stratified medium Green's function. Work in progress.
MODULE strat
  USE mesh
  USE quad
  USE aux

  IMPLICIT NONE

CONTAINS
  SUBROUTINE zbesselj(sorder, norders, z, res)
    INTEGER, INTENT(IN) :: sorder, norders
    COMPLEX (KIND=dp), INTENT(IN) :: z

    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: res
    INTEGER :: nz, err
    DOUBLE PRECISION, DIMENSION(norders) :: yr, yi

    CALL ZBESJ(REAL(z), AIMAG(z), sorder, 1, norders, yr, yi, nz, err)

    !IF(err/=0) THEN
    !   WRITE(*,*) 'Error in Bessel function evaluation, code ', err
    !END IF

    res(:) = yr(:) + (0,1)*yi(:)
  END SUBROUTINE zbesselj

  SUBROUTINE zbesselh(sorder, norders, kind, z, res)
    INTEGER, INTENT(IN) :: sorder, norders, kind
    COMPLEX (KIND=dp), INTENT(IN) :: z

    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: res
    INTEGER :: nz, err
    DOUBLE PRECISION, DIMENSION(norders) :: yr, yi

    CALL ZBESH(REAL(z), AIMAG(z), sorder, 1, kind, norders, yr, yi, nz, err)

    !IF(err/=0) THEN
    !   WRITE(*,*) 'Error in Bessel function evaluation, code ', err
    !END IF

    res(:) = yr(:) + (0,1)*yi(:)
  END SUBROUTINE zbesselh

  ! Returns integral int_0^1 s*exp(x*s) ds.
  FUNCTION intChi(x) RESULT(res)
    COMPLEX (KIND=dp), INTENT(IN) :: x
    COMPLEX (KIND=dp) :: res, y

    IF(x==0) THEN
       res = 0.5_dp
    ELSE
       y = 1.0_dp/(x*x)
       res = EXP(x)*(1.0_dp/x - y) + y
    END IF
  END FUNCTION intChi

  ! Returns integral int_0^1 s^2*exp(x*s) ds.
  FUNCTION intXi(x) RESULT(res)
    COMPLEX (KIND=dp), INTENT(IN) :: x
    COMPLEX (KIND=dp) :: res

    IF(x==0) THEN
       res = 1.0_dp/3.0_dp
    ELSE
       res = (EXP(x) - 2.0_dp*intChi(x))/x
    END IF
  END FUNCTION intXi

  ! Return integral int_0^1 int_0^(1-s) s*exp(cs*s)*exp(ct*t) dtds.
  FUNCTION intExp(cs, ct) RESULT(res)
    COMPLEX (KIND=dp), INTENT(IN) :: cs, ct
    COMPLEX (KIND=dp) :: res

    IF(cs==0) THEN
       IF(ct==0) THEN
          res = 1.0_dp/6.0_dp
       ELSE
          res = (EXP(ct)*intChi(-ct) - 0.5_dp)/ct
       END IF
    ELSE IF(ct==0) THEN
       res = intChi(cs) - intXi(cs)
    ELSE
       res = (EXP(ct)*intChi(cs-ct) - intChi(cs))/ct
    END IF
  END FUNCTION intExp

  SUBROUTINE get_tri_vectors(mesh, m, me, s, t)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: m, me
    REAL (KIND=dp), DIMENSION(3), INTENT(INOUT) :: s, t

    IF(me==1) THEN
       s = mesh%nodes(mesh%faces(m)%node_indices(1))%p &
            - mesh%nodes(mesh%faces(m)%node_indices(3))%p

       t = mesh%nodes(mesh%faces(m)%node_indices(2))%p &
            - mesh%nodes(mesh%faces(m)%node_indices(3))%p
    ELSE IF(me==2) THEN
       s = mesh%nodes(mesh%faces(m)%node_indices(2))%p &
            - mesh%nodes(mesh%faces(m)%node_indices(1))%p

       t = mesh%nodes(mesh%faces(m)%node_indices(3))%p &
            - mesh%nodes(mesh%faces(m)%node_indices(1))%p
    ELSE IF(me==3) THEN
       s = mesh%nodes(mesh%faces(m)%node_indices(3))%p &
            - mesh%nodes(mesh%faces(m)%node_indices(2))%p

       t = mesh%nodes(mesh%faces(m)%node_indices(1))%p &
            - mesh%nodes(mesh%faces(m)%node_indices(2))%p
    ELSE
       WRITE(*,*) 'Invalid local edge index!'
       STOP
    END IF
  END SUBROUTINE get_tri_vectors

  ! Partial moment matrix element <fm,D(fn)> for RWG functions fn.
  ! Integrations are restricted over given triangles.
  ! Operator D is D(f) = int_S G*f dS' with Green's dyadic G for
  ! stratified medium.
  FUNCTION stratMoment(mesh, m, n, me, ne, k0, ri1, ri2, sigma_max) RESULT(mom)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: m, n, me, ne
    REAL (KIND=dp), INTENT(IN) :: sigma_max
    COMPLEX (KIND=dp), INTENT(IN) :: ri1, ri2, k0

    REAL (KIND=dp) :: cm, cn
    REAL (KIND=dp), DIMENSION(3) :: pm, pn, sm, tm, sn, tn
    COMPLEX (KIND=dp) :: intprop, intevan, mom, eps1, eps2, k1, k2, krho

    REAL (KIND=dp), PARAMETER :: eps = 1D-3
    INTEGER, PARAMETER :: maxDepth = 10

    INTEGER, PARAMETER :: npt = 200
    INTEGER :: i, j

    COMPLEX (KIND=dp), DIMENSION(npt,npt) :: plot
    REAL (KIND=dp), DIMENSION(npt,npt) :: x, y

    eps1 = ri1**2
    eps2 = ri2**2

    k1 = k0*ri1
    k2 = k0*ri2

    ! RWG function coefficients multiplied by Jacobian over triangle.
    cm = get_face_sign(m, me, mesh)*mesh%edges(mesh%faces(m)%edge_indices(me))%length
    cn = get_face_sign(n, ne, mesh)*mesh%edges(mesh%faces(n)%edge_indices(ne))%length

    ! Get RWG node for face m.
    IF(get_face_sign(m, me, mesh)>0) THEN
       pm = get_posit_bnode(m, me, mesh)
    ELSE
       pm = get_negat_bnode(m, me, mesh)
    END IF

    ! Get RWG node for face n.
    IF(get_face_sign(n, ne, mesh)>0) THEN
       pn = get_posit_bnode(n, ne, mesh)
    ELSE
       pn = get_negat_bnode(n, ne, mesh)
    END IF

    ! Get triangle edge vectors for RWG functions.
    CALL get_tri_vectors(mesh, m, me, sm, tm)
    CALL get_tri_vectors(mesh, n, ne, sn, tn)

    krho = 0

    DO i=1,npt
       DO j=1,npt
          !x(i,j) = REAL(i-1)/(npt-1)*sigma_max
          !y(i,j) = (REAL(j-1)/(npt-1)*sigma_max - sigma_max/2)*1e-3

          x(i,j) = REAL(i-1)/(npt-1)*sigma_max
          y(i,j) = REAL(j-1)/(npt-1)*2*pi
          plot(i,j) = fevan(x(i,j), y(i,j))
          !plot(i,j) = g( x(i,j) + (0,1)*y(i,j) )
          !plot(i,j) = int_fevan(x(i,j))
       END DO
    END DO

    CALL write_data('x.dat', x)
    CALL write_data('y.dat', y)
    CALL write_data('fevan_re.dat', REAL(plot))
    CALL write_data('fevan_im.dat', AIMAG(plot))

    RETURN

    ! Integrate.
    CALL asqz2(fprop, 0.0_dp, 0.5_dp*pi, 0.0_dp, 2.0_dp*pi, eps, maxDepth, intprop)
    CALL asqz2(fevan, 0.0_dp, sigma_max, 0.0_dp, 2.0_dp*pi, eps, maxDepth, intevan)

    mom = cm*cn*(intprop*(k1**2) + intevan)/(8*(pi**2))

  CONTAINS
    FUNCTION f(phi) RESULT(res)
      REAL (KIND=dp), INTENT(IN) :: phi
      COMPLEX (KIND=dp) :: res, kz1, kz2, rs, rp, intn, intm, ep0
      REAL (KIND=dp) :: sp, cp
      COMPLEX (KIND=dp), DIMENSION(3,3) :: Ms, Mp, M

      sp = SIN(phi)
      cp = COS(phi)

      ! Normal wave vector components.
      kz1 = SQRT(k1**2 - krho**2)

      IF(AIMAG(kz1)<0) THEN
         kz1 = -kz1
      END IF

      kz2 = SQRT(k2**2 - krho**2)

      IF(AIMAG(kz2)<0) THEN
         kz2 = -kz2
      END IF

      ! Reflection coefficients.
      rs = (kz1 - kz2)/(kz1 + kz2)
      rp = (eps2*kz1 - eps1*kz2)/(eps2*kz1 + eps1*kz2)

      ! s-polarized dyad. Column by column.
      Ms = (rs/kz1)*RESHAPE((/sp*sp, -cp*sp, 0.0_dp,  -cp*sp, cp*cp, 0.0_dp,&
           0.0_dp, 0.0_dp, 0.0_dp/), (/3,3/))

      ! p-polarized dyad.
      Mp = -(rp/(k1**2))*RESHAPE((/cp*cp*kz1, cp*sp*kz1, -cp*krho,&
           cp*sp*kz1, sp*sp*kz1, -sp*krho,&
           cp*krho, sp*krho, -(krho**2)/kz1/), (/3,3/))

      M = (Ms + Mp)

      intm = intExp((0,1)*krho*(cp*sm(1) + sp*sm(2)) + (0,1)*kz1*sm(3),&
           (0,1)*krho*(cp*tm(1) + sp*tm(2)) + (0,1)*kz1*tm(3))

      intn = intExp(-(0,1)*krho*(cp*sn(1) + sp*sn(2)) + (0,1)*kz1*sn(3),&
           -(0,1)*krho*(cp*tn(1) + sp*tn(2)) + (0,1)*kz1*tn(3))

      ep0 = EXP((0,1)*krho*(cp*(pm(1)-pn(1)) + sp*(pm(2)-pn(2))) + (0,1)*kz1*(pm(3)+pn(3)))

      res = ep0*intm*intn*(dotrc(sm, MATMUL(M, sn)) + dotrc(sm, MATMUL(M, tn)) +&
           dotrc(tm, MATMUL(M, sn)) + dotrc(tm, MATMUL(M, tn)))

    END FUNCTION f

    FUNCTION g(krho1) RESULT(res)
      COMPLEX (KIND=dp), INTENT(IN) :: krho1
      COMPLEX (KIND=dp) :: res

      krho = krho1

      res = asqz(f, 0.0_dp, 2.0_dp*pi, eps, maxDepth)*krho
    END FUNCTION g

    ! Propagating plane-wave integrand.
    FUNCTION fprop(theta, phi) RESULT(res)
      REAL (KIND=dp), INTENT(IN) :: theta, phi
      REAL (KIND=dp) :: st, ct, sp, cp
      COMPLEX (KIND=dp) :: res, intm, intn, rs, rp, kz1, kz2, ep0
      COMPLEX (KIND=dp), DIMENSION(3,3) :: Ms, Mp, M

      st = SIN(theta)
      ct = COS(theta)
      sp = SIN(phi)
      cp = COS(phi)

      ! Normal wave vector components.
      kz1 = k1*ct
      kz2 = SQRT(k2**2 - (k1*st)**2)

      ! Reflection coefficients.
      rs = (kz1 - kz2)/(kz1 + kz2)
      rp = (eps2*kz1 - eps1*kz2)/(eps2*kz1 + eps1*kz2)

      ! s-polarized dyad. Column by column.
      Ms = rs*RESHAPE((/sp*sp, -cp*sp, 0.0_dp,  -cp*sp, cp*cp, 0.0_dp,&
           0.0_dp, 0.0_dp, 0.0_dp/), (/3,3/))

      ! p-polarized dyad.
      Mp = -rp*RESHAPE((/cp*cp*ct*ct, cp*sp*ct*ct, -cp*st*ct,&
           cp*sp*ct*ct, sp*sp*ct*ct, -sp*st*ct,&
           cp*st*ct, sp*st*ct, -st*st/), (/3,3/))

      M = (Ms + Mp)*st/k1

      intm = intExp((0,1)*k1*st*(cp*sm(1) + sp*sm(2)) + (0,1)*k1*ct*sm(3),&
           (0,1)*k1*st*(cp*tm(1) + sp*tm(2)) + (0,1)*k1*ct*tm(3))

      intn = intExp(-(0,1)*k1*st*(cp*sn(1) + sp*sn(2)) + (0,1)*k1*ct*sn(3),&
           -(0,1)*k1*st*(cp*tn(1) + sp*tn(2)) + (0,1)*k1*ct*tn(3))

      ep0 = EXP((0,1)*k1*st*(cp*(pm(1)-pn(1)) + sp*(pm(2)-pn(2))) + (0,1)*k1*ct*(pm(3)+pn(3)))

      res = ep0*intm*intn*(dotrc(sm, MATMUL(M, sn)) + dotrc(sm, MATMUL(M, tn)) +&
           dotrc(tm, MATMUL(M, sn)) + dotrc(tm, MATMUL(M, tn)))
    END FUNCTION fprop

    ! Evanescent plane-wave integrand.
    !FUNCTION fevan(phi) RESULT(res)
    FUNCTION fevan(sigma, phi) RESULT(res)
      REAL (KIND=dp), INTENT(IN) :: sigma, phi
      !REAL (KIND=dp), INTENT(IN) :: phi
      REAL (KIND=dp) :: sp, cp
      COMPLEX (KIND=dp) :: res, x, y, c, intm, intn, kz1, kz2, rs, rp, ep0
      COMPLEX (KIND=dp), DIMENSION(3,3) :: Ms, Mp, M

      sp = SIN(phi)
      cp = COS(phi)
      x = (0,1)*sigma*sigma
      c = SQRT(sigma*sigma + k1*k1)
      y = c*sigma

      ! Normal wave vector components.
      kz1 = (0,1)*sigma
      kz2 = SQRT(k2**2 - k1**2 - sigma**2)

      ! Reflection coefficients.
      rs = (kz1 - kz2)/(kz1 + kz2)
      rp = (eps2*kz1 - eps1*kz2)/(eps2*kz1 + eps1*kz2)

      ! s-polarized dyad. Column by column.      
      Ms = -(0,1)*rs*RESHAPE((/sp*sp, -cp*sp, 0.0_dp,  -cp*sp, cp*cp, 0.0_dp,&
           0.0_dp, 0.0_dp, 0.0_dp/), (/3,3/))

      ! p-polarized dyad.
      Mp = -(rp/(k1**2))*RESHAPE((/cp*cp*x, cp*sp*x, -cp*y,&
           cp*sp*x, sp*sp*x, -sp*y,&
           cp*y, sp*y, (0,1)*(sigma*sigma + k1*k1)/), (/3,3/))

      M = Ms + Mp

      intm = intExp((0,1)*c*(cp*sm(1) + sp*sm(2)) - sigma*sm(3),&
           (0,1)*c*(cp*tm(1) + sp*tm(2)) - sigma*tm(3))

      intn = intExp(-(0,1)*c*(cp*sn(1) + sp*sn(2)) - sigma*sn(3),&
           -(0,1)*c*(cp*tn(1) + sp*tn(2)) - sigma*tn(3))

      ep0 = EXP((0,1)*c*(cp*(pm(1)-pn(1)) + sp*(pm(2)-pn(2))) - sigma*(pm(3)+pn(3)))

      res = ep0*intm*intn*(dotrc(sm, MATMUL(M, sn)) + dotrc(sm, MATMUL(M, tn)) +&
           dotrc(tm, MATMUL(M, sn)) + dotrc(tm, MATMUL(M, tn)))
    END FUNCTION fevan

!!$    FUNCTION int_fevan(sigma1) RESULT(res)
!!$      REAL (KIND=dp), INTENT(IN) :: sigma1
!!$      COMPLEX (KIND=dp) :: res
!!$
!!$      sigma = sigma1
!!$
!!$      res = asqz(fevan, 0.0_dp, 2.0_dp*pi, eps, maxDepth)
!!$    END FUNCTION int_fevan

  END FUNCTION stratMoment

  FUNCTION integGr(krho, rho, phi, zzp, k1, k2, eps1, eps2, mode) RESULT(Gdyad)
    COMPLEX (KIND=dp), INTENT(IN) :: krho, k1, k2, eps1, eps2
    REAL (KIND=dp), INTENT(IN) :: rho, phi, zzp
    INTEGER, INTENT(IN) :: mode

    COMPLEX (KIND=dp), DIMENSION(3,3) :: Gdyad
    COMPLEX (KIND=dp) :: kz1, kz2, rs, rp, J0, J1, J1b
    COMPLEX (KIND=dp), DIMENSION(2) :: bfs
    REAL (KIND=dp) :: sp, cp, s2p, c2p

    sp = SIN(phi)
    cp = COS(phi)

    s2p = SIN(2*phi)
    c2p = COS(2*phi)

    ! Normal wave vector components.
    kz1 = SQRT(k1**2 - krho**2)

    IF(AIMAG(kz1)<0) THEN
       kz1 = -kz1
    END IF

    kz2 = SQRT(k2**2 - krho**2)

    IF(AIMAG(kz2)<0) THEN
       kz2 = -kz2
    END IF
    
    ! Reflection coefficients.
    rs = (kz1 - kz2)/(kz1 + kz2)
    rp = (eps2*kz1 - eps1*kz2)/(eps2*kz1 + eps1*kz2)

    ! Bessel functions of orders 0, 1 and 2.
    ! bfs(1): order 0, bfs(2): order 1.

    IF(mode==1) THEN
       CALL zbesselj(0, 2, krho*rho, bfs)
    ELSE
       CALL zbesselh(0, 2, mode-1, krho*rho, bfs)
    END IF

    J0 = bfs(1)
    J1 = bfs(2)
    J1b = J1/(krho*rho)

    ! Dyad for p-polarized field.
    Gdyad = -RESHAPE((/(J0*cp*cp-J1b*c2p)*kz1, -(J1b-0.5_dp*J0)*s2p*kz1, -(0,1)*J1*cp*krho,&
         -(J1b-0.5_dp*J0)*s2p*kz1, (J0*sp*sp+J1b*c2p)*kz1, -(0,1)*J1*sp*krho,&
         (0,1)*J1*cp*krho, (0,1)*J1*sp*krho, -(krho**2)/kz1/), (/3,3/))*rp/(k1**2)

    ! Dyad for s-polarized field.
    Gdyad(1:2,1:2) = Gdyad(1:2,1:2) + RESHAPE((/J0*sp*sp+J1b*c2p, (J1b-0.5_dp*J0)*s2p,&
         (J1b-0.5_dp*J0)*s2p, J0*cp*cp-J1b*c2p/), (/2,2/))*rs/kz1

    Gdyad = Gdyad*EXP((0,1)*kz1*zzp)*krho

  END FUNCTION integGr

  SUBROUTINE plot_integGr()
    COMPLEX (KIND=dp) :: k0, ri1, ri2, krho, eps1, eps2, k1, k2
    REAL (KIND=dp) :: rho, phi, zzp, ka, t
    INTEGER, PARAMETER :: npt = 50
    REAL (KIND=dp), DIMENSION(npt*2,18) :: data
    COMPLEX (KIND=dp), DIMENSION(3,3) :: Gdyad
    INTEGER :: n

    k0 = 2*pi/1e-6
    ri1 = 1
    ri2 = 1.45
    eps1 = ri1**2
    eps2 = ri2**2

    rho = 1d-9
    phi = pi/4
    zzp = 0

    k1 = ri1*k0
    k2 = ri2*k0

    ka = (k1+k2)/2


    DO n=1,npt
       t = REAL(n-1)/(npt-1)*pi
       
       krho = ka*(1 - COS(t)) - (0,1)*ka*SIN(t)
       
       Gdyad = integGr(krho, rho, phi, zzp, k1, k2, eps1, eps2, 1)
       
       data(n,1:9) = REAL(RESHAPE(Gdyad,(/9/)))
       data(n,10:18) = AIMAG(RESHAPE(Gdyad,(/9/)))
    END DO


    IF(rho<zzp) THEN
       DO n=1,npt
          t = REAL(n-1)/(npt-1)
          
          krho = ka*2 + t*k1*1500
          
          Gdyad = integGr(krho, rho, phi, zzp, k1, k2, eps1, eps2, 1)
          
          data(n+npt,1:9) = REAL(RESHAPE(Gdyad,(/9/)))
          data(n+npt,10:18) = AIMAG(RESHAPE(Gdyad,(/9/)))
       END DO
    ELSE
       DO n=1,npt
          t = REAL(n-1)/(npt-1)
          
          krho = ka*2 + t*k1*(0,1)*1500
          
          Gdyad = integGr(krho, rho, phi, zzp, k1, k2, eps1, eps2, 2)/2
          Gdyad = Gdyad + integGr(CONJG(krho), rho, phi, zzp, k1, k2, eps1, eps2, 3)/2
          
          data(n+npt,1:9) = REAL(RESHAPE(Gdyad,(/9/)))
          data(n+npt,10:18) = AIMAG(RESHAPE(Gdyad,(/9/)))
       END DO
    END IF

    CALL write_data('strat.dat', data)
  END SUBROUTINE plot_integGr

  SUBROUTINE test_strat()
    TYPE(mesh_container) :: mesh
    INTEGER m, n, me, ne, i
    COMPLEX (KIND=dp) :: k0, ri1, ri2
    INTEGER, PARAMETER :: npt = 20
    REAL (KIND=dp), DIMENSION(npt,2) :: data
    REAL (KIND=dp) :: sigma_max

    k0 = 2*pi/1e-6
    ri1 = 1
    ri2 = 1.45

    m = 1442
    n = 1523
    me = 1
    ne = 2

    mesh = load_mesh('triangle.msh')
    CALL build_mesh(mesh, 1d-9)

    !DO i=1,npt
    !   sigma_max = k0/100*i
    !   data(i,1) = i
    !   data(i,2) = ABS(stratMoment(mesh, m, n, me, ne, k0, ri1, ri2, sigma_max))
    !END DO

    !CALL write_data('strat.dat', data)

    data(1,1) = stratMoment(mesh, m, n, me, ne, k0, ri1, ri2, 1200*REAL(k0))

    CALL delete_mesh(mesh)
  END SUBROUTINE test_strat
END MODULE strat
