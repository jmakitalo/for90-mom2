! MODULE: int
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Low-level routines for computing various singular integrals over
! line segments and triangle faces. The integrands consist of product of RWG
! function and kernel of the type R^m, where m may be negative integer.
! Based on article Hanninen 2006 Progress in Electromagnetic Research.
MODULE int
  USE mesh
  USE linalg

  IMPLICIT NONE

CONTAINS
  FUNCTION intLm1(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res
    REAL (KIND=dp) :: Rp, Rm, sp, sm, div1, div2
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, s

    IF(edgeind==1) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
       s = mesh%faces(faceind)%s(:,1)
    ELSEIF(edgeind==2) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p
       s = mesh%faces(faceind)%s(:,2)
    ELSEIF(edgeind==3) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
       s = mesh%faces(faceind)%s(:,3)
    ELSE
       WRITE(*,*) 'Error evaluating intLm1!'
       STOP
    END IF

    Rp = normr(r-p2)
    Rm = normr(r-p1)
    sp = dotr((p2-r), s)
    sm = dotr((p1-r), s)

    div1 = Rp-sp
    div2 = Rm+sm

    IF(ABS(div1)>ABS(div2)) THEN
       res = LOG((Rm-sm)/div1)
    ELSE
       res = LOG((Rp+sp)/div2)
    END IF
  END FUNCTION intLm1

  FUNCTION intL1(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res
    REAL (KIND=dp) :: Rp, Rm, sp, sm, h, t, R02
    REAL (KIND=dp), DIMENSION(3) :: p, p1, p2, s

    IF(edgeind==1) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
       p = 0.5_dp*(p1 + p2)
    ELSEIF(edgeind==2) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p
       p = 0.5_dp*(p1 + p2)
    ELSEIF(edgeind==3) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
       p = 0.5_dp*(p1 + p2)
    ELSE
       WRITE(*,*) 'Error evaluating intL1!'
       STOP
    END IF

    h = dotr(mesh%faces(faceind)%n, r-p)
    t = dotr(mesh%faces(faceind)%m(:,edgeind), r-p)
    R02 = (t**2) + (h**2)

    s = mesh%faces(faceind)%s(:,edgeind)

    Rp = normr(r-p2)
    Rm = normr(r-p1)
    sp = dotr((p2-r), s)
    sm = dotr((p1-r), s)

    IF(R02==0) THEN
       res = 0.5_dp*(sp*Rp - sm*Rm)
    ELSE
       res = 0.5_dp*R02*intLm1(r, faceind, edgeind, mesh) + 0.5_dp*(sp*Rp - sm*Rm)
    END IF
  END FUNCTION intL1

  FUNCTION intL3(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res
    REAL (KIND=dp) :: Rp, Rm, sp, sm, h, t, R02
    REAL (KIND=dp), DIMENSION(3) :: p, p1, p2, s

    IF(edgeind==1) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
       p = 0.5_dp*(p1 + p2)
    ELSEIF(edgeind==2) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p
       p = 0.5_dp*(p1 + p2)
    ELSEIF(edgeind==3) THEN
       p1 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p
       p2 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
       p = 0.5_dp*(p1 + p2)
    ELSE
       WRITE(*,*) 'Error evaluating intL3!'
       STOP
    END IF

    h = dotr(mesh%faces(faceind)%n, r-p)
    t = dotr(mesh%faces(faceind)%m(:,edgeind), r-p)
    R02 = (t**2) + (h**2)

    s = mesh%faces(faceind)%s(:,edgeind)

    Rp = normr(r-p2)
    Rm = normr(r-p1)
    sp = dotr((p2-r), s)
    sm = dotr((p1-r), s)

    res = 0.75_dp*R02*intL1(r, faceind, edgeind, mesh) + 0.25_dp*(sp*(Rp**3) - sm*(Rm**3))
  END FUNCTION intL3

  FUNCTION intSm3h(r, faceind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res, x, y, omega
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, p3, a1, a2, a3

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
    p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
    p3 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p

    a1 = (p1-r)/normr(p1-r)
    a2 = (p2-r)/normr(p2-r)
    a3 = (p3-r)/normr(p3-r)

    x = 1.0_dp + dotr(a1,a2) + dotr(a1,a3) + dotr(a2,a3)
    y = ABS(dotr(a1, crossr(a2,a3)))

    omega = 2.0_dp*ATAN2(y,x)

    IF(dotr(mesh%faces(faceind)%n, r-p1)>0.0_dp) THEN
       res = ABS(omega)
    ELSE
       res = -ABS(omega)
    END IF
  END FUNCTION intSm3h

  FUNCTION intSm3(r, faceind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res, h
    REAL (KIND=dp), DIMENSION(3) :: p1

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
    h = dotr(mesh%faces(faceind)%n, r-p1)

    IF(h==0) THEN
       res = 0
    ELSE
       res = intSm3h(r,faceind,mesh)/h
    END IF
  END FUNCTION intSm3

  FUNCTION intSm1(r, faceind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res, h, t1, t2, t3
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, p3
    REAL (KIND=dp) :: edgeint1, edgeint2, edgeint3

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
    p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
    p3 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p
    t1 = dotr(mesh%faces(faceind)%m(:,1), r-p1)
    t2 = dotr(mesh%faces(faceind)%m(:,2), r-p2)
    t3 = dotr(mesh%faces(faceind)%m(:,3), r-p3)

    h = dotr(mesh%faces(faceind)%n, r-p1)

    IF(h==0) THEN
       IF(t1==0) THEN
          edgeint1 = 0
       ELSE
          edgeint1 = intLm1(r, faceind, 1, mesh)
       END IF
    
       IF(t2==0) THEN
          edgeint2 = 0
       ELSE
          edgeint2 = intLm1(r, faceind, 2, mesh)
       END IF
       
       IF(t3==0) THEN
          edgeint3 = 0
       ELSE
          edgeint3 = intLm1(r, faceind, 3, mesh)
       END IF

       res = -(t1*edgeint1 + t2*edgeint2 + t3*edgeint3)

       ! Changed 15.8.2012 to permit observation points on edges.
       !res = -(t1*intLm1(r, faceind, 1, mesh) +&
       !     t2*intLm1(r, faceind, 2, mesh) + t3*intLm1(r, faceind, 3, mesh))
    ELSE
       res = -h*intSm3h(r, faceind, mesh) -(t1*intLm1(r, faceind, 1, mesh) +&
            t2*intLm1(r, faceind, 2, mesh) + t3*intLm1(r, faceind, 3, mesh))
    END IF
  END FUNCTION intSm1

  FUNCTION intS1(r, faceind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res, h, t1, t2, t3
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, p3

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
    p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
    p3 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p
    t1 = dotr(mesh%faces(faceind)%m(:,1), r-p1)
    t2 = dotr(mesh%faces(faceind)%m(:,2), r-p2)
    t3 = dotr(mesh%faces(faceind)%m(:,3), r-p3)

    h = dotr(mesh%faces(faceind)%n, r-p1)

    res = (h**2)/3.0_dp*intSm1(r, faceind, mesh) - 1.0_dp/3.0_dp*(t1*intL1(r, faceind, 1, mesh) +&
         t2*intL1(r, faceind, 2, mesh) + t3*intL1(r, faceind, 3, mesh))
  END FUNCTION intS1

  FUNCTION intK11(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res, L, A, sign

    sign = get_face_sign(faceind, edgeind, mesh)

    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    res = sign*L/A*intS1(r, faceind, mesh)
  END FUNCTION intK11

  FUNCTION intK1m1(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: res, L, A, sign

    sign = get_face_sign(faceind, edgeind, mesh)

    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    res = sign*L/A*intSm1(r, faceind, mesh)
  END FUNCTION intK1m1

  FUNCTION intK21(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: L, A, sign, h
    REAL (KIND=dp), DIMENSION(3) :: p, rho, res, p1
    REAL (KIND=dp), DIMENSION(3,3) :: m

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p

    sign = get_face_sign(faceind, edgeind, mesh)
    p = get_face_bnode(faceind, edgeind, mesh)

    h = dotr(mesh%faces(faceind)%n, r-p1)
    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    m = mesh%faces(faceind)%m

    rho = r - h*mesh%faces(faceind)%n

    res = sign*L/(2.0_dp*A)*(1.0_dp/3.0_dp*(m(:,1)*intL3(r, faceind, 1, mesh) +&
         m(:,2)*intL3(r, faceind, 2, mesh) +&
         m(:,3)*intL3(r, faceind, 3, mesh)) +&
         (rho-p)*intS1(r, faceind, mesh))
  END FUNCTION intK21

  FUNCTION intK2m1(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: L, A, sign, h
    REAL (KIND=dp), DIMENSION(3) :: p, rho, res, p1
    REAL (KIND=dp), DIMENSION(3,3) :: m

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p

    sign = get_face_sign(faceind, edgeind, mesh)
    p = get_face_bnode(faceind, edgeind, mesh)

    h = dotr(mesh%faces(faceind)%n, r-p1)
    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    m = mesh%faces(faceind)%m

    rho = r - h*mesh%faces(faceind)%n

    res = sign*L/(2.0_dp*A)*(m(:,1)*intL1(r, faceind, 1, mesh) +&
         m(:,2)*intL1(r, faceind, 2, mesh) +&
         m(:,3)*intL1(r, faceind, 3, mesh) +&
         (rho-p)*intSm1(r, faceind, mesh))
  END FUNCTION intK2m1

  FUNCTION intK31(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: h
    REAL (KIND=dp), DIMENSION(3) :: p, n, res
    REAL (KIND=dp), DIMENSION(3,3) :: m

    p = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p

    n = mesh%faces(faceind)%n
    h = dotr(n, r-p)
    m = mesh%faces(faceind)%m

    IF(h==0) THEN
       res = m(:,1)*intL1(r, faceind, 1, mesh) +&
            m(:,2)*intL1(r, faceind, 2, mesh) +&
            m(:,3)*intL1(r, faceind, 3, mesh)
    ELSE
       res = m(:,1)*intL1(r, faceind, 1, mesh) +&
            m(:,2)*intL1(r, faceind, 2, mesh) +&
            m(:,3)*intL1(r, faceind, 3, mesh) -&
            h*n*intSm1(r, faceind, mesh)
    END IF
  END FUNCTION intK31

  FUNCTION intK3m1(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: h
    REAL (KIND=dp), DIMENSION(3) :: p, n, res
    REAL (KIND=dp), DIMENSION(3,3) :: m

    p = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p

    n = mesh%faces(faceind)%n
    h = dotr(n, r-p)
    m = mesh%faces(faceind)%m

    res = m(:,1)*intLm1(r, faceind, 1, mesh) +&
         m(:,2)*intLm1(r, faceind, 2, mesh) +&
         m(:,3)*intLm1(r, faceind, 3, mesh) +&
         n*intSm3h(r, faceind, mesh)
  END FUNCTION intK3m1

  FUNCTION intK41(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3) :: p, res
    REAL (KIND=dp) :: L, A, sign

    sign = get_face_sign(faceind, edgeind, mesh)
    p = get_face_bnode(faceind, edgeind, mesh)

    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    res = -sign*L/(2.0_dp*A)*crossr((r-p), intK31(r, faceind, edgeind, mesh))
  END FUNCTION intK41

  FUNCTION intK4m1(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3) :: p, res
    REAL (KIND=dp) :: L, A, sign

    sign = get_face_sign(faceind, edgeind, mesh)
    p = get_face_bnode(faceind, edgeind, mesh)

    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    res = -sign*L/(2.0_dp*A)*crossr((r-p), intK3m1(r, faceind, edgeind, mesh))
  END FUNCTION intK4m1
END MODULE int
