! MODULE: srcint
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for the evaluation of various boundary integral operators.
! These routines utilize the singularity subtraction routines to
! handle accurately cases, where observation point r is near or coincides
! with source point r'.
MODULE srcint
  USE quad
  USE int
  USE rwgf
  USE green
  USE greenprd
  USE symmetry
  USE OMP_LIB

  IMPLICIT NONE

CONTAINS
  ! Computes the singular part of intK1.
  FUNCTION intK1singular(r, faceind, edgeind, mesh, k, prd) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd

    COMPLEX (KIND=dp) :: res, aint, phase
    REAL (KIND=dp), DIMENSION(3) :: rho
    INTEGER :: n, m

    res = 0.0_dp

    IF(ASSOCIATED(prd)) THEN
       phase = 1.0_dp
       
       DO n=-1,1
          DO m=-1,1
             rho = (/n*prd%dx*prd%cp, m*prd%dy + n*prd%dx*prd%sp, 0.0_dp/)
             IF(prd%oblique) THEN
                phase = exp((0,1)*(prd%coef(prd%cwl)%k0x*rho(1) + prd%coef(prd%cwl)%k0y*rho(2)))
             END IF

             aint = 1.0_dp/(4.0_dp*pi)*intK1m1(r-rho, faceind, edgeind, mesh) &
                  - (k**2)/(8.0_dp*pi)*intK11(r-rho, faceind, edgeind, mesh)
             res = res + aint*phase
          END DO
       END DO
    ELSE
       res = 1.0_dp/(4.0_dp*pi)*intK1m1(r, faceind, edgeind, mesh) &
            - (k**2)/(8.0_dp*pi)*intK11(r, faceind, edgeind, mesh)
    END IF
  END FUNCTION intK1singular

  ! Returns int_T (O_g'G(r,r'))div'_t(f(r')) dS' for all three RWG functions f
  ! defined over triangle T. Singular case r in T is accurately computed.
  FUNCTION intK1(r, faceind, mesh, k, ga, prd, near, qd) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp), DIMENSION(3) :: res

    INTEGER :: t, nweights, edgeind
    REAL (KIND=dp) :: An
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn
    REAL (KIND=dp) :: divfn
    COMPLEX (KIND=dp), DIMENSION(qd%num_nodes) :: gv
    COMPLEX (KIND=dp) :: tmp, aint
    LOGICAL :: singular

    singular = .FALSE.

    !IF(near .AND. ga%id==gid_identity) THEN
    IF(near) THEN
       singular = .TRUE.
    END IF

    nweights = qd%num_nodes

    An = mesh%faces(faceind)%area
    qpn = quad_tri_points(qd, faceind, mesh)

    IF(ga%id/=gid_identity) THEN
       DO t=1,nweights
          qpn(:,t) = MATMUL(ga%j, qpn(:,t))
       END DO
    END IF

    IF(ASSOCIATED(prd)) THEN
       CALL vGp(r, qpn, prd, nweights, singular, gv)
    ELSE
       IF(singular) THEN
          CALL vGs(r, qpn, k, nweights, gv)
       ELSE
          CALL vGf(r, qpn, k, nweights, gv)
       END IF
    END IF

    tmp = An*SUM(qd%weights*gv)

    DO edgeind=1,3
       divfn = rwgDiv(faceind,edgeind,mesh)

       aint = 0.0_dp

       IF(singular) THEN
          aint = intK1singular(MATMUL(TRANSPOSE(ga%j),r), faceind, edgeind, mesh, k, prd)
       END IF

       res(edgeind) = tmp*divfn + aint
    END DO
  END FUNCTION intK1

  ! Computes the singular part of intK2.
  FUNCTION intK2singular(r, faceind, edgeind, mesh, k, prd) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd

    COMPLEX (KIND=dp), DIMENSION(3) :: res, aint
    COMPLEX (KIND=dp) :: phase
    REAL (KIND=dp), DIMENSION(3) :: rho
    INTEGER :: n, m

    res(:) = 0.0_dp

    IF(ASSOCIATED(prd)) THEN
       phase = 1.0_dp
       
       DO n=-1,1
          DO m=-1,1
             rho = (/n*prd%dx*prd%cp, m*prd%dy + n*prd%dx*prd%sp, 0.0_dp/)
             IF(prd%oblique) THEN
                phase = exp((0,1)*(prd%coef(prd%cwl)%k0x*rho(1) + prd%coef(prd%cwl)%k0y*rho(2)))
             END IF

             aint = 1.0_dp/(4.0_dp*pi)*intK2m1(r-rho, faceind, edgeind, mesh) &
                  - (k**2)/(8.0_dp*pi)*intK21(r-rho, faceind, edgeind, mesh)
             res = res + aint*phase
          END DO
       END DO
    ELSE
       res = 1.0_dp/(4.0_dp*pi)*intK2m1(r, faceind, edgeind, mesh) &
            - (k**2)/(8.0_dp*pi)*intK21(r, faceind, edgeind, mesh)
    END IF
  END FUNCTION intK2singular

  ! Returns int_T (O_g'G(r,r'))(J_g f(r')) dS' for all three RWG functions f
  ! defined over triangle T. Singular case r in T is accurately computed.
  FUNCTION intK2(r, faceind, mesh, k, ga, prd, near, qd) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp), DIMENSION(3,3) :: res
    COMPLEX (KIND=dp), DIMENSION(3) :: aint
    INTEGER :: edgeind, nweights, t
    REAL (KIND=dp) :: An
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn, qpn2
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: fv
    COMPLEX (KIND=dp), DIMENSION(qd%num_nodes) :: gv
    LOGICAL :: singular

    singular = .FALSE.

    !IF(near .AND. ga%id==gid_identity) THEN
    IF(near) THEN
       singular = .TRUE.
    END IF

    nweights = qd%num_nodes

    An = mesh%faces(faceind)%area
    qpn = quad_tri_points(qd, faceind, mesh)

    IF(ga%id/=gid_identity) THEN
       DO t=1,nweights
          qpn2(:,t) = MATMUL(ga%j, qpn(:,t))
       END DO
    ELSE
       qpn2 = qpn
    END IF

    IF(ASSOCIATED(prd)) THEN
       CALL vGp(r, qpn2, prd, nweights, singular, gv)
    ELSE
       IF(singular) THEN
          CALL vGs(r, qpn2, k, nweights, gv)
       ELSE
          CALL vGf(r, qpn2, k, nweights, gv)
       END IF
    END IF

    gv = gv*qd%weights*An

    DO edgeind=1,3
       CALL vrwg(qpn(:,:),faceind,edgeind,mesh,fv)

       aint = 0.0_dp

       IF(singular) THEN
          aint = intK2singular(MATMUL(TRANSPOSE(ga%j),r), faceind, edgeind, mesh, k, prd)
       END IF
       
       res(:,edgeind) = MATMUL(fv,gv) + aint

       IF(ga%id/=gid_identity) THEN
          res(:,edgeind) = MATMUL(ga%j, res(:,edgeind))
       END IF
    END DO
  END FUNCTION intK2

  ! Computes the singular part of intK3.
  FUNCTION intK3singular(r, faceind, edgeind, mesh, k, prd) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd

    COMPLEX (KIND=dp), DIMENSION(3) :: res, aint
    COMPLEX (KIND=dp) :: phase
    REAL (KIND=dp), DIMENSION(3) :: rho
    INTEGER :: n, m

    res(:) = 0.0_dp

    IF(ASSOCIATED(prd)) THEN
       phase = 1.0_dp
       
       DO n=-1,1
          DO m=-1,1
             rho = (/n*prd%dx*prd%cp, m*prd%dy + n*prd%dx*prd%sp, 0.0_dp/)
             IF(prd%oblique) THEN
                phase = exp((0,1)*(prd%coef(prd%cwl)%k0x*rho(1) + prd%coef(prd%cwl)%k0y*rho(2)))
             END IF

             aint = 1.0_dp/(4.0_dp*pi)*intK3m1(r-rho, faceind, edgeind, mesh) &
                  - (k**2)/(8.0_dp*pi)*intK31(r-rho, faceind, edgeind, mesh)
             res = res + aint*phase
          END DO
       END DO
    ELSE
       res = 1.0_dp/(4.0_dp*pi)*intK3m1(r, faceind, edgeind, mesh) &
            - (k**2)/(8.0_dp*pi)*intK31(r, faceind, edgeind, mesh)
    END IF
  END FUNCTION intK3singular

  ! Returns int_T (O_g'grad'G(r,r'))div'_tf(r') dS' for all three RWG functions f
  ! defined over triangle T. Singular case r in T is accurately computed.
  FUNCTION intK3(r, faceind, mesh, k, ga, prd, near, qd) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp), DIMENSION(3,3) :: res
    COMPLEX (KIND=dp), DIMENSION(3) :: aint
    INTEGER :: t, nweights, edgeind
    REAL (KIND=dp) :: An, divfn
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes) :: ggv
    COMPLEX (KIND=dp), DIMENSION(3) :: tmp
    LOGICAL :: singular

    singular = .FALSE.

    !IF(near .AND. ga%id==gid_identity) THEN
    IF(near) THEN
       singular = .TRUE.
    END IF

    nweights = qd%num_nodes

    An = mesh%faces(faceind)%area
    qpn = quad_tri_points(qd, faceind, mesh)

    IF(ga%id/=gid_identity) THEN
       DO t=1,nweights
          qpn(:,t) = MATMUL(ga%j, qpn(:,t))
       END DO
    END IF

    res(:,:) = 0.0_dp

    IF(ASSOCIATED(prd)) THEN
       CALL vgradGp(r, qpn, prd, nweights, singular, ggv)
    ELSE
       IF(singular) THEN
          CALL vgradGs(r, qpn, k, nweights, ggv)
       ELSE
          CALL vgradGf(r, qpn, k, nweights, ggv)
       END IF
    END IF

    tmp = MATMUL(ggv,qd%weights)*An

    DO edgeind=1,3
       divfn = rwgDiv(faceind,edgeind,mesh)

       aint = 0.0_dp

       IF(singular) THEN
          aint = intK3singular(MATMUL(TRANSPOSE(ga%j),r), faceind, edgeind, mesh, k, prd)*divfn
          aint = MATMUL(ga%j, aint)
       END IF

       res(:,edgeind) = tmp*divfn + aint
    END DO
  END FUNCTION intK3

  ! Computes the singular part of intK4.
  FUNCTION intK4singular(r, faceind, edgeind, mesh, k, prd, tfaceind) RESULT(res)
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, edgeind, tfaceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd

    COMPLEX (KIND=dp), DIMENSION(3) :: res, aint
    COMPLEX (KIND=dp) :: phase
    REAL (KIND=dp), DIMENSION(3) :: rho
    INTEGER :: n, m

    res(:) = 0.0_dp

    IF(ASSOCIATED(prd)) THEN
       phase = 1.0_dp
       
       DO n=-1,1
          DO m=-1,1
             IF(n/=0 .OR. m/=0 .OR. faceind/=tfaceind) THEN
                rho = (/n*prd%dx*prd%cp, m*prd%dy + n*prd%dx*prd%sp, 0.0_dp/)
                IF(prd%oblique) THEN
                   phase = exp((0,1)*(prd%coef(prd%cwl)%k0x*rho(1) + prd%coef(prd%cwl)%k0y*rho(2)))
                END IF

                aint = 1.0_dp/(4.0_dp*pi)*intK4m1(r-rho, faceind, edgeind, mesh) &
                     - (k**2)/(8.0_dp*pi)*intK41(r-rho, faceind, edgeind, mesh)
                res = res + aint*phase
             END IF
          END DO
       END DO
    ELSE
       res = 1.0_dp/(4.0_dp*pi)*intK4m1(r, faceind, edgeind, mesh) &
            - (k**2)/(8.0_dp*pi)*intK41(r, faceind, edgeind, mesh)
    END IF

  END FUNCTION intK4singular

  ! Returns int_T (O_g'grad'G(r,r'))x(J_g f(r')) dS' for all three RWG functions f
  ! defined over triangle T. Singular case r in T is accurately computed.
  FUNCTION intK4(r, faceind, mesh, k, ga, tfaceind, prd, near, qd) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: k
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: faceind, tfaceind
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    LOGICAL, INTENT(IN) :: near
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp), DIMENSION(3,3) :: res
    COMPLEX (KIND=dp), DIMENSION(3) :: aint
    INTEGER :: t, nweights, n, m, edgeind
    REAL (KIND=dp) :: An
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn, qpn2
    REAL (KIND=dp), DIMENSION(3) :: rho
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: fv
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes) :: ggv
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes) :: tmp
    COMPLEX (KIND=dp) :: phase
    LOGICAL :: singular

    IF(ga%id==gid_identity .AND. faceind==tfaceind .AND. ASSOCIATED(prd)==.FALSE.) THEN
       res(:,:) = 0.0_dp
       RETURN
    END IF

    singular = .FALSE.

    !IF(near .AND. ga%id==gid_identity) THEN
    IF(near) THEN
       singular = .TRUE.
    END IF

    nweights = qd%num_nodes

    An = mesh%faces(faceind)%area
    qpn = quad_tri_points(qd, faceind, mesh)

    IF(ga%id/=gid_identity) THEN
       DO t=1,nweights
          qpn2(:,t) = MATMUL(ga%j, qpn(:,t))
       END DO
    ELSE
       qpn2 = qpn
    END IF

    IF(ASSOCIATED(prd)) THEN
       CALL vgradGp(r, qpn2, prd, nweights, singular, ggv)
    ELSE
       IF(singular) THEN
          CALL vgradGs(r, qpn2, k, nweights, ggv)
       ELSE
          CALL vgradGf(r, qpn2, k, nweights, ggv)
       END IF
    END IF

    DO edgeind=1,3
       CALL vrwg(qpn,faceind,edgeind,mesh,fv)

       IF(ga%id/=gid_identity) THEN
          DO t=1,nweights
             fv(:,t) = MATMUL(ga%j, fv(:,t))
          END DO
       END IF

       CALL vcrosscr(ggv, fv, nweights, tmp)

       aint = 0.0_dp

       IF(singular) THEN
          aint = intK4singular(MATMUL(TRANSPOSE(ga%j),r), faceind, edgeind, mesh, k, prd, tfaceind)
          aint = MATMUL(ga%j, aint)/ga%detj
       END IF

       res(:,edgeind) = MATMUL(tmp,qd%weights)*An + aint
    END DO
  END FUNCTION intK4
END MODULE srcint
