! MODULE: rwgf
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for the evaluation of RWG basis functions.
MODULE rwgf
  USE mesh
  USE quad

  IMPLICIT NONE

CONTAINS
  FUNCTION rwg(r, faceind, edgeind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3) :: res, v
    REAL (KIND=dp) :: L, A

    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    IF(get_face_sign(faceind, edgeind, mesh)>0) THEN
       v = get_posit_bnode(faceind, edgeind, mesh)
       res = L/(2.0_dp*A)*(r-v)
    ELSE IF(get_face_sign(faceind, edgeind, mesh)<0) THEN
       v = get_negat_bnode(faceind, edgeind, mesh)
       res = L/(2.0_dp*A)*(v-r)
    ELSE
       WRITE(*,*) 'Invalid evaluation of RWG function!'
    END IF
  END FUNCTION rwg

  FUNCTION rwgDiv(faceind, edgeind, mesh) RESULT(res)
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: L, A, res, sign

    IF(get_face_sign(faceind, edgeind, mesh)>0) THEN
       sign = 1.0_dp
    ELSE IF(get_face_sign(faceind, edgeind, mesh)<0) THEN
       sign = -1.0_dp
    ELSE
       WRITE(*,*) 'Invalid evaluation of RWG divergence function!'
       STOP
    END IF

    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    res = sign*L/A
  END FUNCTION rwgDiv

  SUBROUTINE vrwg(r, faceind, edgeind, mesh, res)
    REAL (KIND=dp), DIMENSION(:,:), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3) :: v
    REAL (KIND=dp), DIMENSION(SIZE(r,1),SIZE(r,2)), INTENT(INOUT) :: res
    REAL (KIND=dp) :: L, A

    L = mesh%edges(mesh%faces(faceind)%edge_indices(edgeind))%length
    A = mesh%faces(faceind)%area

    IF(get_face_sign(faceind, edgeind, mesh)>0) THEN
       v = get_posit_bnode(faceind, edgeind, mesh)
       res(1,:) = L/(2.0_dp*A)*(r(1,:)-v(1))
       res(2,:) = L/(2.0_dp*A)*(r(2,:)-v(2))
       res(3,:) = L/(2.0_dp*A)*(r(3,:)-v(3))
    ELSE IF(get_face_sign(faceind, edgeind, mesh)<0) THEN
       v = get_negat_bnode(faceind, edgeind, mesh)
       res(1,:) = L/(2.0_dp*A)*(v(1)-r(1,:))
       res(2,:) = L/(2.0_dp*A)*(v(2)-r(2,:))
       res(3,:) = L/(2.0_dp*A)*(v(3)-r(3,:))
    ELSE
       WRITE(*,*) 'Invalid evaluation of RWG function!'
    END IF
  END SUBROUTINE vrwg

  FUNCTION solid_rwg(r, solidind, faceind, mesh) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: solidind, faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3) :: res, v
    REAL (KIND=dp) :: A, vol

    A = mesh%solid_faces(mesh%solids(solidind)%solid_face_indices(faceind))%area
    vol = mesh%solids(solidind)%volume

    IF(get_solid_face_sign(solidind, faceind, mesh)>0) THEN
       v = get_posit_solid_bnode(solidind, faceind, mesh)
       res = A/(3.0_dp*vol)*(r-v)
    ELSE IF(get_solid_face_sign(solidind, faceind, mesh)<0) THEN
       v = get_negat_solid_bnode(solidind, faceind, mesh)
       res = A/(3.0_dp*vol)*(v-r)
    ELSE
       WRITE(*,*) 'Invalid evaluation of solid RWG function!'
    END IF
  END FUNCTION solid_rwg

  FUNCTION solid_rwgDiv(solidind, faceind, mesh) RESULT(res)
    INTEGER, INTENT(IN) :: solidind, faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp) :: A, vol, res, sign

    IF(get_solid_face_sign(solidind, faceind, mesh)>0) THEN
       sign = 1.0_dp
    ELSE IF(get_solid_face_sign(solidind, faceind, mesh)<0) THEN
       sign = -1.0_dp
    ELSE
       WRITE(*,*) 'Invalid evaluation of solid RWG divergence function!'
       STOP
    END IF

    A = mesh%solid_faces(mesh%solids(solidind)%solid_face_indices(faceind))%area
    vol = mesh%solids(solidind)%volume

    res = sign*A/vol
  END FUNCTION solid_rwgDiv

  SUBROUTINE rwg_moments(mesh, F)
    TYPE(mesh_container), INTENT(IN) :: mesh

    COMPLEX (KIND=dp), DIMENSION(mesh%nedges,mesh%nedges), INTENT(INOUT) :: F
    INTEGER :: nbasis, nweights, n, m, q, p, r, t, s, faceind
    REAL (KIND=dp) :: A, int
    REAL (KIND=dp), DIMENSION(3) :: fm, fn
    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: qp

    WRITE(*,*) 'Building an F-matrix'

    F(:,:) = 0.0_dp
    nbasis = mesh%nedges
    nweights = SIZE(qw)

    DO m=1,mesh%nedges

       DO s=1,2
          faceind = mesh%edges(m)%face_indices(s)

          IF(faceind==-1) THEN
             CYCLE
          END IF

          qp = GLquad_points(faceind, mesh)
          A = mesh%faces(faceind)%area

          q = local_edge_index(mesh, faceind, m)
          
          int = 0.0_dp
          DO r=1,nweights
             fm = rwg(qp(:,r),faceind,q,mesh)
                
             int = int + qw(r)*dotr(fm, fm)
          END DO
          int = int*A
          
          F(m,m) = F(m,m) + int
       END DO

       DO n=m+1,mesh%nedges
          ! Find out the index of the common face.

          IF(mesh%edges(m)%face_indices(1)/=-1 .AND. (mesh%edges(m)%face_indices(1)==mesh%edges(n)%face_indices(1) .OR.&
               mesh%edges(m)%face_indices(1)==mesh%edges(n)%face_indices(2))) THEN
             faceind = mesh%edges(m)%face_indices(1)
          ELSE IF(mesh%edges(m)%face_indices(2)/=-1 .AND. (mesh%edges(m)%face_indices(2)==mesh%edges(n)%face_indices(1) .OR.&
               mesh%edges(m)%face_indices(2)==mesh%edges(n)%face_indices(2))) THEN
             faceind = mesh%edges(m)%face_indices(2)
          ELSE
             CYCLE
          END IF

          q = local_edge_index(mesh, faceind, m)
          p = local_edge_index(mesh, faceind, n)

          qp = GLquad_points(faceind, mesh)
          A = mesh%faces(faceind)%area

          int = 0.0_dp
          DO r=1,nweights
             fm = rwg(qp(:,r),faceind,q,mesh)
             fn = rwg(qp(:,r),faceind,p,mesh)
                
             int = int + qw(r)*dotr(fm, fn)
          END DO
          int = int*A

          F(m,n) = int
          F(n,m) = int
       END DO
    END DO
  END SUBROUTINE rwg_moments

  SUBROUTINE rwg_moments2(mesh, F)
    TYPE(mesh_container), INTENT(IN) :: mesh

    COMPLEX (KIND=dp), DIMENSION(mesh%nedges,mesh%nedges), INTENT(INOUT) :: F
    INTEGER :: nbasis, nweights, n, m, q, p, r, t, s, faceind
    REAL (KIND=dp) :: A, int
    REAL (KIND=dp), DIMENSION(3) :: fm, fn
    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: qp

    WRITE(*,*) 'Building an F-matrix'

    F(:,:) = 0.0_dp
    nbasis = mesh%nedges
    nweights = SIZE(qw)

    DO m=1,mesh%nedges

       DO n=m+1,mesh%nedges
          ! Find out the index of the common face.

          IF(mesh%edges(m)%face_indices(1)/=-1 .AND. (mesh%edges(m)%face_indices(1)==mesh%edges(n)%face_indices(1) .OR.&
               mesh%edges(m)%face_indices(1)==mesh%edges(n)%face_indices(2))) THEN
             faceind = mesh%edges(m)%face_indices(1)
          ELSE IF(mesh%edges(m)%face_indices(2)/=-1 .AND. (mesh%edges(m)%face_indices(2)==mesh%edges(n)%face_indices(1) .OR.&
               mesh%edges(m)%face_indices(2)==mesh%edges(n)%face_indices(2))) THEN
             faceind = mesh%edges(m)%face_indices(2)
          ELSE
             CYCLE
          END IF

          q = local_edge_index(mesh, faceind, m)
          p = local_edge_index(mesh, faceind, n)

          qp = GLquad_points(faceind, mesh)
          A = mesh%faces(faceind)%area

          int = 0.0_dp
          DO r=1,nweights
             fm = rwg(qp(:,r),faceind,q,mesh)
             fn = rwg(qp(:,r),faceind,p,mesh)
                
             int = int + qw(r)*dotr(crossr(mesh%faces(faceind)%n, fm), fn)
          END DO
          int = int*A

          F(m,n) = int
          F(n,m) = -int
       END DO
    END DO
  END SUBROUTINE rwg_moments2

  ! Evaluate a function expanded in RWG functions over
  ! the given surface mesh.
  FUNCTION rwg_exp(pt, mesh, faceind, coef) RESULT(res)
    INTEGER, INTENT(IN) :: faceind
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: coef
    TYPE(mesh_container), INTENT(IN) :: mesh

    COMPLEX (KIND=dp), DIMENSION(3) :: res
    REAL (KIND=dp), DIMENSION(3) :: fm
    INTEGER :: q, index

    res(:) = 0.0_dp

    DO q=1,3
       fm = rwg(pt, faceind, q, mesh)

       index = mesh%faces(faceind)%edge_indices(q)

       res = res + fm*coef(index)
    END DO

  END FUNCTION rwg_exp

  ! Evaluate the surface divergence of a function expanded in RWG functions over
  ! the given surface mesh.
  FUNCTION div_rwg_exp(mesh, faceind, coef) RESULT(res)
    INTEGER, INTENT(IN) :: faceind
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: coef
    TYPE(mesh_container), INTENT(IN) :: mesh

    COMPLEX (KIND=dp) :: res
    REAL (KIND=dp) :: fmDiv
    INTEGER :: q, index

    res = 0.0_dp

    DO q=1,3
       fmDiv = rwgDiv(faceind, q, mesh)

       index = mesh%faces(faceind)%edge_indices(q)

       res = res + fmDiv*coef(index)
    END DO

  END FUNCTION div_rwg_exp
END MODULE rwgf
