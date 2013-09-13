! MODULE: nlsurf
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Surface second-order nonlinear source functions.
MODULE nlsurf
  USE rwgf
  USE symmetry

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

  ! Computes the inverse action on field, i.e., O_g^-1(E) where
  ! g is group element of index gai and O_g E(r) = E(p_g(r)) where
  ! p_g is group action on points.
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

       e = e + ga(gai)%ef(n)*MATMUL(TRANSPOSE(ga(gai)%j), e2)
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

       !Pnls = Pnls + CONJG(ga(n)%ef(nf)**2)*MATMUL(ga(n)%j, Pnls2)
       Pnls = Pnls + CONJG(ga(n)%ef(nf))*MATMUL(ga(n)%j, Pnls2)
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
       nf, nls, src_coef, src_vec)
    TYPE(mesh_container), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: nedgestot, nf
    REAL (KIND=dp), INTENT(IN) :: omegaff
    COMPLEX (KIND=dp), INTENT(IN) :: riff, rish
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: epsp
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: xff
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    TYPE(medium_nls), INTENT(IN) :: nls
   
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: src_coef, src_vec
    COMPLEX (KIND=dp), DIMENSION(mesh%nedges) :: coef
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: F
    INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
    INTEGER :: INFO, m, r, q, index, nweights, nbasis, m2, t

    COMPLEX (KIND=dp) :: int1, int2, int3
    REAL (KIND=dp) :: A, fmDiv
    COMPLEX (KIND=dp), DIMENSION(3,SIZE(qw)) :: Pnls_tan
    COMPLEX (KIND=dp), DIMENSION(SIZE(qw)) :: Pnls_n
    COMPLEX (KIND=dp), DIMENSION(3) :: Pnls
    REAL (KIND=dp), DIMENSION(3) :: nor, fm, p1, p2, p
    REAL (KIND=dp), DIMENSION(3,SIZE(qw)) :: qp
    REAL (KIND=dp) :: omegash
    REAL (KIND=dp), DIMENSION(2) :: pts
    LOGICAL :: contour_method

    contour_method = .FALSE.

    omegash = 2.0_dp*omegaff

    WRITE(*,*) 'Computing second-order surface source basis coefficients'

    nweights = SIZE(qw)
    nbasis = mesh%nedges

    coef(:) = 0.0_dp
    src_vec(:) = 0.0_dp
    src_coef(:) = 0.0_dp

    ! Points for n=2 Gauss-Legendre quadrature.
    pts(1) = -1.0_dp/SQRT(3.0_dp)
    pts(2) = -pts(1)

    DO m=1,mesh%nfaces
       qp = GLquad_points(m, mesh)
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

             int1 = int1 - qw(r)*fmDiv*Pnls_n(r)
             int2 = int2 + qw(r)*dotc(CMPLX(fm,KIND=dp), crossc(Pnls_tan(:,r), CMPLX(nor,KIND=dp)))
             int3 = int3 + qw(r)*dotc(CMPLX(fm,KIND=dp), Pnls_tan(:,r))
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

    ALLOCATE(F(nbasis,nbasis), IPIV(nbasis))

    CALL rwg_moments(mesh, F)

    CALL ZGETRF(nbasis, nbasis, F, nbasis, IPIV, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Matrix factorization failed!'
       STOP
    END IF

    CALL ZGETRS('N', nbasis, 1, F, nbasis, IPIV, coef, nbasis, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Solving of system failed!'
       STOP
    END IF

    CALL ZGETRS('N', nbasis, 1, F, nbasis, IPIV, src_coef((nbasis+1):(2*nbasis)), nbasis, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Solving of system failed!'
       STOP
    END IF

    IF(contour_method) THEN
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
          qp = GLquad_points(m, mesh)
          A = mesh%faces(m)%area
          nor = mesh%faces(m)%n
          
          DO q=1,3
             
             int1 = 0.0_dp
             DO r=1,nweights
                fm = rwg(qp(:,r), m, q, mesh)
                
                int1 = int1 + qw(r)*dotc(CMPLX(fm,KIND=dp),&
                     crossc(CMPLX(nor,KIND=dp), rwg_exp(qp(:,r), mesh, m, coef)))
             END DO
             
             index = mesh%faces(m)%edge_indices(q)
             
             src_coef(index) = src_coef(index) + int1*A
             
          END DO
       END DO       
    END IF

    CALL ZGETRS('N', nbasis, 1, F, nbasis, IPIV, src_coef(1:nbasis), nbasis, INFO)
    IF(INFO/=0) THEN
       WRITE(*,*) 'Solving of system failed!'
       STOP
    END IF

    DEALLOCATE(F, IPIV)
  END SUBROUTINE nlsurf_coef

END MODULE nlsurf
