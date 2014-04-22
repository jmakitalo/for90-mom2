! MODULE: nlsurf
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Surface second-order nonlinear source functions.
MODULE nlsurf
  USE rwgf
  USE symmetry
  USE bc
  USE srcint

  IMPLICIT NONE

  TYPE medium_nls
     COMPLEX (KIND=dp) :: chi2_nnn
     COMPLEX (KIND=dp) :: chi2_ntt
     COMPLEX (KIND=dp) :: chi2_ttn
  END TYPE medium_nls

CONTAINS
  ! Evaluate E-field on surface.
  FUNCTION efield(mesh, nedgestot, omega, ri, x, faceind, r) RESULT(e)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: nedgestot, faceind
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3) :: e, ms
    COMPLEX (KIND=dp) :: en, divj, eps
    INTEGER :: n, edgeind

    eps = (ri**2)*eps0

    divj = 0.0_dp
    ms(:) = 0.0_dp

    DO n=1,3
       edgeind = mesh%faces(faceind)%edge_indices(n)
       edgeind = mesh%edges(edgeind)%parent_index

       divj = divj + x(edgeind)*rwgDiv(faceind, n, mesh)

       ms = ms + x(nedgestot + edgeind)*rwg(r, faceind, n, mesh)
    END DO

    en = divj/((0,1)*omega*eps)

    e = crossc(CMPLX(mesh%faces(faceind)%n,KIND=dp), ms) + mesh%faces(faceind)%n*en
  END FUNCTION efield

  ! Computes the inverse action on field, i.e., E(J_g^-1 r) where
  ! J_g is a group action of index gai.
  FUNCTION efield_invmap(mesh, nedgestot, omega, ri, x, ga, faceind, gai, r) RESULT(e)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot, faceind, gai
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3) :: e, e2
    INTEGER :: n

    e(:) = 0.0_dp

    DO n=1,SIZE(ga)
       e2 = efield(mesh, nedgestot, omega, ri, x(:,n), faceind, r)

       e = e + CONJG(ga(gai)%ef(n))*MATMUL(TRANSPOSE(ga(gai)%j), e2)
    END DO
  END FUNCTION efield_invmap

  ! Computes the inverse action on nonlinear surface polarization Pnls.
  ! See comments for efield_invmap.
  FUNCTION Pnls_invmap(mesh, nedgestot, omega, ri, x, ga, faceind, gai, nls, r) RESULT(Pnls)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot, faceind, gai
    TYPE(medium_nls), INTENT(IN) :: nls
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3) :: Pnls, Pnls_tan, e, et
    COMPLEX (KIND=dp) :: Pnls_n, en
    REAL (KIND=dp), DIMENSION(3) :: nor, tan

    e = efield_invmap(mesh, nedgestot, omega, ri, x, ga, faceind, gai, r)

    nor = MATMUL(TRANSPOSE(ga(gai)%j), mesh%faces(faceind)%n)

    en = dotc(CMPLX(nor,KIND=dp), e)
    et = e - en*nor

    Pnls_n = eps0*(nls%chi2_nnn*(en**2) + nls%chi2_ntt*dotc(et,et))
    Pnls_tan = 2.0_dp*eps0*nls%chi2_ttn*et*en

    Pnls = Pnls_tan + Pnls_n*nor
  END FUNCTION Pnls_invmap

  ! Computes the nonlinear surface polarization Pnls in representation (fragment) nf.
  FUNCTION Pnls_frag(mesh, nedgestot, omega, ri, x, ga, faceind, nf, nls, r) RESULT(Pnls)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(medium_nls), INTENT(IN) :: nls
    INTEGER, INTENT(IN) :: nedgestot, faceind, nf
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: r

    COMPLEX (KIND=dp), DIMENSION(3) :: Pnls, Pnls2
    INTEGER :: n

    Pnls(:) = 0.0_dp

    DO n=1,SIZE(ga)
       Pnls2 = Pnls_invmap(mesh, nedgestot, omega, ri, x, ga, faceind, n, nls, r)

       Pnls = Pnls + ga(n)%ef(nf)*MATMUL(ga(n)%j, Pnls2)
    END DO

    Pnls = Pnls/REAL(SIZE(ga),KIND=dp)
  END FUNCTION Pnls_frag

  ! src_coef = (cM, cJ)
  ! sum_n cM_n <f_m,f_n> = <f_m,n x grad'Pn>/epsp
  ! sum_n cJ_n <f_m,f_n> = -i*Omega*<f_m,Pt>
  ! N.B. src_coef does not depend on direction of normal n.
  ! src_vec = (E0, H0)
  ! E0_n = -<div'f_m,Pn>/(2*epsp)
  ! H0_n = i*Omega*<f_m,Pt x n>/2
  SUBROUTINE nlsurf_coef(mesh, nedgestot, omegaff, riff, rish, epsp, xff, ga,&
       nf, nls, qd, phdx, phdy, src_coef, src_vec)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot, nf
    REAL (KIND=dp), INTENT(IN) :: omegaff
    COMPLEX (KIND=dp), INTENT(IN) :: riff, rish, phdx, phdy
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: epsp
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: xff
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(medium_nls), INTENT(IN) :: nls
    TYPE(quad_data), INTENT(IN) :: qd
   
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: src_coef, src_vec
    COMPLEX (KIND=dp), DIMENSION(mesh%nedges) :: coef
    COMPLEX (KIND=dp), DIMENSION(mesh%nedges,2) :: tmp
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: F
    INTEGER :: INFO, m, r, q, index, nweights, nbasis, m2, t, dim

    COMPLEX (KIND=dp) :: int1, int2, int3
    REAL (KIND=dp) :: A, fmDiv
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes) :: Pnls_tan
    COMPLEX (KIND=dp), DIMENSION(qd%num_nodes) :: Pnls_n
    COMPLEX (KIND=dp), DIMENSION(3) :: Pnls
    REAL (KIND=dp), DIMENSION(3) :: nor, fm, p1, p2, p
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qp
    REAL (KIND=dp) :: omegash
    REAL (KIND=dp), DIMENSION(2) :: pts
    LOGICAL :: contour_method
    INTEGER, DIMENSION(mesh%nedges*2) :: id
    COMPLEX (KIND=dp), DIMENSION(mesh%nedges*2) :: phase

    contour_method = .FALSE.

    omegash = 2.0_dp*omegaff

    WRITE(*,*) 'Computing second-order surface source basis coefficients'

    nweights = qd%num_nodes
    nbasis = mesh%nedges

    coef(:) = 0.0_dp
    src_vec(:) = 0.0_dp
    src_coef(:) = 0.0_dp

    DO m=1,mesh%nfaces
       qp = quad_tri_points(qd, m, mesh)
       A = mesh%faces(m)%area
       nor = mesh%faces(m)%n

       DO r=1,nweights
          Pnls = Pnls_frag(mesh, nedgestot, omegaff, riff, xff, ga, m, nf, nls, qp(:,r))

          Pnls_n(r) = dotc(CMPLX(nor,KIND=dp), Pnls)
          Pnls_tan(:,r) = Pnls - Pnls_n(r)*nor
       END DO

       DO q=1,3

          fmDiv = rwgDiv(m, q, mesh)
          
          int1 = 0.0_dp
          int2 = 0.0_dp
          int3 = 0.0_dp
          DO r=1,nweights
             fm = rwg(qp(:,r), m, q, mesh)

             int1 = int1 - qd%weights(r)*fmDiv*Pnls_n(r)
             int2 = int2 + qd%weights(r)*dotc(CMPLX(fm,KIND=dp), crossc(Pnls_tan(:,r), CMPLX(nor,KIND=dp)))
             int3 = int3 + qd%weights(r)*dotc(CMPLX(fm,KIND=dp), Pnls_tan(:,r))
          END DO

          index = mesh%faces(m)%edge_indices(q)

          coef(index) = coef(index) + int1*A/epsp(m)

          src_coef(index + nbasis) = src_coef(index + nbasis) + int3*A

          src_vec(index + nbasis) = src_vec(index + nbasis) + int2*A

       END DO
    END DO

    src_coef((nbasis+1):(2*nbasis)) = -src_coef((nbasis+1):(2*nbasis))*(0,1)*omegash

    src_vec(1:nbasis) = 0.5_dp*coef(1:nbasis)

    src_vec((nbasis+1):(2*nbasis)) = 0.5_dp*(0,1)*omegash*src_vec((nbasis+1):(2*nbasis))

    ALLOCATE(F(nbasis,nbasis))

    CALL rwg_moments(mesh, qd, F)

    ! Arrange two RHS of a linear system into one array.
    ! These vectors are symmetrical in the same way as the E-field.
    tmp(:,1) = coef
    tmp(:,2) = src_coef((nbasis+1):(2*nbasis))

    ! Impose boundary conditions on the system.
    CALL edge_bc(mesh, ga, phdx, phdy, nf, id, phase)
    CALL resolve_system_dependencies(F, tmp, id(1:nbasis), phase(1:nbasis))
    CALL reduce_system(F, tmp, dim, id(1:nbasis))

    ! Solve the system.
    CALL solve_multi_linsys(F(1:dim,1:dim), tmp(1:dim,:))

    ! Expand the solution vectors to size nbasis and place zeros
    ! where designated by symmetry.
    CALL expand_solution(dim, id(1:nbasis), phase(1:nbasis), tmp)

    ! Place the solved coefficients back.
    coef = tmp(:,1)
    src_coef((nbasis+1):(2*nbasis)) = tmp(:,2)

    ! This contour method was described by C. Forestiere, but is inconvenient
    ! from the points of view of symmetry.
    ! The other method is based on expanding nabla_t(P_n) in RWG-basis.
    IF(contour_method) THEN
       ! Points for n=2 Gauss-Legendre quadrature.
       pts(1) = -1.0_dp/SQRT(3.0_dp)
       pts(2) = -pts(1)

       DO m=1,mesh%nfaces
          nor = mesh%faces(m)%n
          
          DO q=1,3             
             ! Contour integration.
             
             int1 = 0.0_dp
             DO t=1,3
                p1 = mesh%nodes(mesh%faces(m)%node_indices(t))%p
                p2 = mesh%nodes(mesh%faces(m)%node_indices(indexrot3(t+1)))%p
                
                m2 = adjacent_face(mesh, m, t)
                
                DO r=1,2
                   p = p1 + (p2-p1)*(pts(r)*0.5_dp + 0.5_dp)
                   fm = rwg(p, m, q, mesh)
                   
                   Pnls = Pnls_frag(mesh, nedgestot, omegaff, riff, xff, ga, m, nf, nls, p)
                   Pnls_n(1) = dotc(CMPLX(nor,KIND=dp), Pnls)
                   
                   Pnls = Pnls_frag(mesh, nedgestot, omegaff, riff, xff, ga, m2, nf, nls, p)
                   Pnls_n(2) = dotc(CMPLX(mesh%faces(m2)%n,KIND=dp), Pnls)
                   
                   int1 = int1 + 0.5_dp*dotr(fm, p2-p1)*(Pnls_n(1) + Pnls_n(2))
                END DO
             END DO

             index = mesh%faces(m)%edge_indices(q)
             
             src_coef(index) = src_coef(index) + 0.5_dp*int1/epsp(m)
          END DO
       END DO
    ELSE
       DO m=1,mesh%nfaces
          qp = quad_tri_points(qd, m, mesh)
          A = mesh%faces(m)%area
          nor = mesh%faces(m)%n
          
          DO q=1,3
             
             int1 = 0.0_dp
             DO r=1,nweights
                fm = rwg(qp(:,r), m, q, mesh)
                
                int1 = int1 + qd%weights(r)*dotc(CMPLX(fm,KIND=dp),&
                     crossc(CMPLX(nor,KIND=dp), rwg_exp(qp(:,r), mesh, m, coef)))
             END DO
             
             index = mesh%faces(m)%edge_indices(q)
             
             src_coef(index) = src_coef(index) + int1*A
             
          END DO
       END DO       
    END IF

    ! Solve another system and use a tmp array to match dimensions for function calls.
    ! This RHS obeys H-field symmetry.
    tmp(:,1) = src_coef(1:nbasis)
    tmp(:,2) = 0.0_dp

    ! F was reduced by E-field symmetry previously, so need to compute it again.
    CALL rwg_moments(mesh, qd, F)

    ! Transform indexes from nbasis+n -> n.
    DO m=1,nbasis
       IF(id(nbasis+m)>nbasis) THEN
          id(nbasis+m) = id(nbasis+m) - nbasis
       END IF
    END DO

    CALL resolve_system_dependencies(F, tmp, id((nbasis+1):(2*nbasis)),&
         phase((nbasis+1):(2*nbasis)))
    CALL reduce_system(F, tmp, dim, id((nbasis+1):(2*nbasis)))
    CALL solve_linsys(F(1:dim,1:dim), tmp(1:dim,1))
    CALL expand_solution(dim, id((nbasis+1):(2*nbasis)), phase((nbasis+1):(2*nbasis)), tmp)

    src_coef(1:nbasis) = tmp(:,1)

    DEALLOCATE(F)
  END SUBROUTINE nlsurf_coef

  FUNCTION nls_surface1(mesh, pt, face_index, k, phi, qd) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt
    INTEGER, INTENT(IN) :: face_index
    COMPLEX (KIND=dp), INTENT(IN) :: k
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: phi
    TYPE(quad_data), INTENT(IN) :: qd
    
    COMPLEX (KIND=dp), DIMENSION(3) :: res
    COMPLEX (KIND=dp) :: g
    INTEGER :: n, t
    REAL (KIND=dp) :: An
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn
    REAL (KIND=dp), DIMENSION(3) :: nor
    
    res(:) = 0.0_dp
    
    DO n=1,mesh%nfaces
       IF(n==face_index) THEN
          CYCLE
       END IF
       
       qpn = quad_tri_points(qd, n, mesh)
       An = mesh%faces(n)%area
       nor = mesh%faces(n)%n
       
       DO t=1,qd%num_nodes
          g = Gf(pt, qpn(:,t), k)
          
          res = res + qd%weights(t)*g*phi(t,n)*nor*An
       END DO
    END DO
    
    res = res*(k**2)
  END FUNCTION nls_surface1

  FUNCTION nls_surface2(mesh, pt, face_index, k, phi, qd) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt
    INTEGER, INTENT(IN) :: face_index
    COMPLEX (KIND=dp), INTENT(IN) :: k
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: phi
    TYPE(quad_data), INTENT(IN) :: qd

    COMPLEX (KIND=dp) :: res
    COMPLEX (KIND=dp) :: g
    INTEGER :: n, t
    REAL (KIND=dp) :: An
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpn
    REAL (KIND=dp), DIMENSION(3) :: nor
    
    res = 0.0_dp
    
    DO n=1,mesh%nfaces
       IF(n==face_index) THEN
          CYCLE
       END IF
       
       qpn = quad_tri_points(qd, n, mesh)
       An = mesh%faces(n)%area
       nor = mesh%faces(n)%n
       
       DO t=1,qd%num_nodes
          g = dotc(CMPLX(nor,KIND=dp), gradGf(pt, qpn(:,t), k))
          
          res = res + qd%weights(t)*g*phi(t,n)*An
       END DO
    END DO
  END FUNCTION nls_surface2

  FUNCTION nls_contour1(mesh, pt, face_index, k, phi) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt
    INTEGER, INTENT(IN) :: face_index
    COMPLEX (KIND=dp), INTENT(IN) :: k
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: phi
    
    COMPLEX (KIND=dp), DIMENSION(3) :: res
    COMPLEX (KIND=dp) :: g
    INTEGER :: n, t, l
    REAL (KIND=dp), DIMENSION(3) :: nor, edgenor, p1, p2, p
    REAL (KIND=dp), DIMENSION(2) :: pts
    
    res(:) = 0.0_dp

    ! Points for n=2 Gauss-Legendre quadrature.
    pts(1) = -1.0_dp/SQRT(3.0_dp)
    pts(2) = -pts(1)
    
    DO n=1,mesh%nfaces
       IF(n==face_index) THEN
          CYCLE
       END IF
       
       nor = mesh%faces(n)%n
       
       DO l=1,3
          p1 = mesh%nodes(mesh%faces(n)%node_indices(l))%p
          p2 = mesh%nodes(mesh%faces(n)%node_indices(indexrot3(l+1)))%p
          
          edgenor = mesh%faces(n)%m(:,l)
          
          DO t=1,SIZE(pts)
             p = p1 + (p2-p1)*(pts(t)*0.5_dp + 0.5_dp)
             
             g = dotc(CMPLX(nor,KIND=dp), gradGf(pt, p, k))
             
             res = res + g*phi(t,l,n)*edgenor*normr(p2-p1)*0.5_dp
          END DO
       END DO
    END DO
  END FUNCTION nls_contour1

  FUNCTION nls_contour2(mesh, pt, face_index, k, phi) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt
    INTEGER, INTENT(IN) :: face_index
    COMPLEX (KIND=dp), INTENT(IN) :: k
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: phi
    
    COMPLEX (KIND=dp), DIMENSION(3) :: res
    COMPLEX (KIND=dp) :: g
    INTEGER :: n, t, l
    REAL (KIND=dp), DIMENSION(3) :: nor, edgenor, p1, p2, p
    REAL (KIND=dp), DIMENSION(2) :: pts
    
    res(:) = 0.0_dp

    ! Points for n=2 Gauss-Legendre quadrature.
    pts(1) = -1.0_dp/SQRT(3.0_dp)
    pts(2) = -pts(1)

    DO n=1,mesh%nfaces
       IF(n==face_index) THEN
          CYCLE
       END IF
       
       nor = mesh%faces(n)%n
       
       DO l=1,3
          p1 = mesh%nodes(mesh%faces(n)%node_indices(l))%p
          p2 = mesh%nodes(mesh%faces(n)%node_indices(indexrot3(l+1)))%p
          
          edgenor = mesh%faces(n)%m(:,l)
         
          DO t=1,SIZE(pts)
             p = p1 + (p2-p1)*(pts(t)*0.5_dp + 0.5_dp)
             
             g = dotc(CMPLX(edgenor,KIND=dp), gradGf(pt, p, k))
             
             res = res + g*phi(t,l,n)*nor*normr(p2-p1)*0.5_dp
          END DO
       END DO
    END DO
  END FUNCTION nls_contour2

  FUNCTION nls_contour3(mesh, pt, face_index, k, phi) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: pt
    INTEGER, INTENT(IN) :: face_index
    COMPLEX (KIND=dp), INTENT(IN) :: k
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(IN) :: phi
    
    COMPLEX (KIND=dp), DIMENSION(3) :: res
    COMPLEX (KIND=dp) :: g
    INTEGER :: n, t, l
    REAL (KIND=dp), DIMENSION(3) :: nor, p1, p2, p
    REAL (KIND=dp), DIMENSION(2) :: pts
    
    res(:) = 0.0_dp

    ! Points for n=2 Gauss-Legendre quadrature.
    pts(1) = -1.0_dp/SQRT(3.0_dp)
    pts(2) = -pts(1)

    DO n=1,mesh%nfaces
       DO l=1,3
          p1 = mesh%nodes(mesh%faces(n)%node_indices(l))%p
          p2 = mesh%nodes(mesh%faces(n)%node_indices(indexrot3(l+1)))%p
         
          DO t=1,SIZE(pts)
             p = p1 + (p2-p1)*(pts(t)*0.5_dp + 0.5_dp)
             
             g = Gf(pt, p, k)
             
             res = res + g*phi(t,l,n)*mesh%faces(n)%s(:,l)*normr(p2-p1)*0.5_dp
          END DO
       END DO
    END DO
  END FUNCTION nls_contour3

  SUBROUTINE nlsurf_srcvec(mesh, nedgestot, omegaff, riff, rish, epsp, xff, ga,&
       nf, nls, qd, src_vec)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot, nf
    REAL (KIND=dp), INTENT(IN) :: omegaff
    COMPLEX (KIND=dp), INTENT(IN) :: riff, rish
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: epsp
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: xff
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(medium_nls), INTENT(IN) :: nls
    TYPE(quad_data), INTENT(IN) :: qd
   
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: src_vec
    INTEGER :: m, m2, r, q, index, nweights, nbasis, t, i

    COMPLEX (KIND=dp) :: int1, k
    REAL (KIND=dp) :: A, fmDiv
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes) :: Pnls_tan
    COMPLEX (KIND=dp), DIMENSION(qd%num_nodes,mesh%nfaces) :: Pnls_n
    COMPLEX (KIND=dp), DIMENSION(2,3,mesh%nfaces) :: Pnls_n2
    COMPLEX (KIND=dp), DIMENSION(3) :: Pnls, int2, intaux2
    COMPLEX (KIND=dp), DIMENSION(3,3) :: intaux1
    REAL (KIND=dp), DIMENSION(3) :: nor, fm, p1, p2, p
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qp
    REAL (KIND=dp) :: omegash
    REAL (KIND=dp), DIMENSION(2) :: pts
    TYPE(prdnfo), POINTER :: prd

    prd => NULL()

    omegash = 2.0_dp*omegaff

    k = rish*omegash/c0

    ! Points for n=2 Gauss-Legendre quadrature.
    pts(1) = -1.0_dp/SQRT(3.0_dp)
    pts(2) = -pts(1)

    WRITE(*,*) 'Computing second-order surface source vector'

    nweights = qd%num_nodes
    nbasis = mesh%nedges

    ! Pre-compute polarization.
    ! Pre-divide by epsp.
    DO m=1,mesh%nfaces
       qp = quad_tri_points(qd, m, mesh)
       nor = mesh%faces(m)%n

       ! Surface integration points.
       DO r=1,nweights
          Pnls = Pnls_frag(mesh, nedgestot, omegaff, riff, xff, ga, m, nf, nls, qp(:,r))

          Pnls_n(r,m) = dotc(CMPLX(nor,KIND=dp), Pnls)/epsp(m)
          !Pnls_tan(:,r) = Pnls - Pnls_n(r)*nor
       END DO

       ! Contour integration points.
       DO i=1,3
          p1 = mesh%nodes(mesh%faces(m)%node_indices(i))%p
          p2 = mesh%nodes(mesh%faces(m)%node_indices(indexrot3(i+1)))%p

          DO r=1,SIZE(pts)
             p = p1 + (p2-p1)*(pts(r)*0.5_dp + 0.5_dp)

             Pnls = Pnls_frag(mesh, nedgestot, omegaff, riff, xff, ga, m, nf, nls, p)

             Pnls_n2(r,i,m) = dotc(CMPLX(nor,KIND=dp), Pnls)/epsp(m)
          END DO
       END DO
    END DO

    src_vec(:) = 0.0_dp

    ! Compute direct E-source.
    DO m=1,mesh%nfaces
       qp = quad_tri_points(qd, m, mesh)
       A = mesh%faces(m)%area

       DO q=1,3

          fmDiv = rwgDiv(m, q, mesh)
          
          int1 = 0.0_dp

          DO r=1,nweights
             fm = rwg(qp(:,r), m, q, mesh)

             int1 = int1 + qd%weights(r)*fmDiv*Pnls_n(r,m)
          END DO

          index = mesh%faces(m)%edge_indices(q)

          src_vec(index) = src_vec(index) - 0.5_dp*int1*A
       END DO
    END DO

    ! Compute integrated E-source.
    DO m=1,mesh%nfaces
       qp = quad_tri_points(qd, m, mesh)
       A = mesh%faces(m)%area

       DO q=1,3
          int1 = 0.0_dp

          DO r=1,nweights
             fm = rwg(qp(:,r), m, q, mesh)

             int2 = nls_surface1(mesh, qp(:,r), m, k, Pnls_n, qd) &
                  - nls_contour1(mesh, qp(:,r), m, k, Pnls_n2)&
                  + nls_contour2(mesh, qp(:,r), m, k, Pnls_n2)

             int1 = int1 + qd%weights(r)*dotc(CMPLX(fm,KIND=dp), int2)

             int1 = int1 + qd%weights(r)*rwgDiv(m, q, mesh)*nls_surface2(mesh, qp(:,r), m, k,&
                  Pnls_n, qd)
          END DO

          index = mesh%faces(m)%edge_indices(q)

          src_vec(index) = src_vec(index) + int1*A
       END DO
    END DO

    DO m=1,mesh%nfaces
       int2(:) = 0.0_dp

       DO i=1,3
          p1 = mesh%nodes(mesh%faces(m)%node_indices(i))%p
          p2 = mesh%nodes(mesh%faces(m)%node_indices(indexrot3(i+1)))%p
          
          DO r=1,SIZE(pts)
             p = p1 + (p2-p1)*(pts(r)*0.5_dp + 0.5_dp)

             int1 = nls_surface2(mesh, p, m, k, Pnls_n, qd)

             DO q=1,3
                fm = rwg(p, m, q, mesh)

                int2(q) = int2(q) + int1*dotr(fm,mesh%faces(m)%m(:,i))*normr(p2-p1)*0.5_dp

                int2(q) = int2(q) + 0.5_dp*Pnls_n2(r,i,m)*dotr(fm,mesh%faces(m)%m(:,i))*normr(p2-p1)*0.5_dp
             END DO
          END DO
       END DO

       DO q=1,3
          index = mesh%faces(m)%edge_indices(q)
          
          src_vec(index) = src_vec(index) - int2(q)
       END DO
    END DO

    ! Compute integrated H-source.
!!$    DO m2=1,mesh%nfaces
!!$       DO m=1,mesh%nfaces
!!$          int2(:) = 0.0_dp
!!$          
!!$          DO i=1,3
!!$             p1 = mesh%nodes(mesh%faces(m)%node_indices(i))%p
!!$             p2 = mesh%nodes(mesh%faces(m)%node_indices(indexrot3(i+1)))%p
!!$             
!!$             DO r=1,SIZE(pts)
!!$                p = p1 + (p2-p1)*(pts(r)*0.5_dp + 0.5_dp) + mesh%faces(m)%n*1d-10
!!$                
!!$                intaux1 = intK2(p, m2, mesh, k, ga(1), prd, .TRUE.)
!!$                
!!$                DO q=1,3
!!$                   int2(q) = int2(q) + dotc(CMPLX((p2-p1)*0.5_dp,KIND=dp),&
!!$                        intaux1(:,q))*Pnls_n2(r,i,m)
!!$                END DO
!!$             END DO
!!$          END DO
!!$          
!!$          DO q=1,3
!!$             index = mesh%faces(m2)%edge_indices(q) + nbasis
!!$             
!!$             src_vec(index) = src_vec(index) + int2(q)
!!$          END DO
!!$       END DO
!!$    END DO

    DO m2=1,mesh%nfaces
       qp = quad_tri_points(qd, m2, mesh)
       A = mesh%faces(m2)%area
       
       int2(:) = 0.0_dp

       DO r=1,qd%num_nodes
          intaux2 = nls_contour3(mesh, qp(:,r), m2, k, Pnls_n2)
             
          DO q=1,3
             fm = rwg(qp(:,r), m2, q, mesh)
             int2(q) = int2(q) + qd%weights(r)*dotc(CMPLX(fm,KIND=dp), intaux2)
          END DO
       END DO
       
       DO q=1,3
          index = mesh%faces(m2)%edge_indices(q) + nbasis
          
          src_vec(index) = src_vec(index) + A*int2(q)
       END DO
    END DO

    DO m2=1,mesh%nfaces
       qp = quad_tri_points(qd, m2, mesh)
       A = mesh%faces(m2)%area
       
       int2(:) = 0.0_dp

       DO r=1,qd%num_nodes
          intaux2(:) = 0.0_dp

          DO m=1,mesh%nfaces
             intaux1 = intK3(qp(:,r), m, mesh, k, ga(1), prd, .TRUE., qd)
             intaux2 = intaux2 + crossc(CMPLX(mesh%faces(m)%n,KIND=dp),&
                  intaux1(:,1)/rwgDiv(m, 1, mesh)*Pnls_n(1,m))
          END DO
             
          DO q=1,3
             fm = rwg(qp(:,r), m2, q, mesh)
             int2(q) = int2(q) + qd%weights(r)*dotc(CMPLX(fm,KIND=dp), intaux2)
          END DO
       END DO
       
       DO q=1,3
          index = mesh%faces(m2)%edge_indices(q) + nbasis
          
          src_vec(index) = src_vec(index) - A*int2(q)
       END DO
    END DO

    src_vec((nbasis+1):(2*nbasis)) = src_vec((nbasis+1):(2*nbasis))*(0,1)*omegash*(rish**2)*eps0
  END SUBROUTINE nlsurf_srcvec

END MODULE nlsurf
