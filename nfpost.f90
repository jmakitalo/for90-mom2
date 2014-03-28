! MODULE: nfpost
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Post-processing routines for near-field data (for e.g. visualization purposes).
MODULE nfpost
  USE nfields
  USE source

  IMPLICIT NONE

CONTAINS
  SUBROUTINE field_stream(name, mesh, scale, nedgestot, x, ga, omega, ri, prd, src, addsrc)
    CHARACTER (LEN=*), INTENT(IN) :: name
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: scale, omega
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(srcdata), INTENT(IN) :: src
    LOGICAL, INTENT(IN) :: addsrc

    INTEGER, PARAMETER :: max_points = 100, nsteps = 4
    COMPLEX (KIND=dp), DIMENSION(SIZE(x,1),SIZE(x,2),1) :: x2
    COMPLEX (KIND=dp), DIMENSION(3,1) :: ef, hf
    COMPLEX (KIND=dp), DIMENSION(3) :: einc, hinc
    REAL (KIND=dp), DIMENSION(3) :: efr
    REAL (KIND=dp), DIMENSION(max_points,mesh%nfaces,nsteps) :: efmag
    REAL (KIND=dp), DIMENSION(3,max_points,mesh%nfaces,nsteps) :: pt
    INTEGER, DIMENSION(mesh%nfaces,nsteps) :: npt
    REAL (KIND=dp) :: dx
    REAL (KIND=dp), DIMENSION(3) :: tan
    INTEGER :: n, m, s

    dx = mesh%avelen

    x2 = RESHAPE(x,(/SIZE(x,1),SIZE(x,2),1/))

    DO s=1,nsteps
       DO n=1,mesh%nfaces
          pt(:,1,n,s) = mesh%faces(n)%cp + mesh%faces(n)%n*dx/10
          
          DO m=1,(max_points-1)
             
             npt(n,s) = m
             
             CALL scat_fields(mesh, ga, x2, nedgestot,&
                  omega, ri, prd, pt(:,m,n,s), ef, hf)
             
             IF(addsrc) THEN
                CALL src_fields(src, omega, ri, pt(:,m,n,s), einc, hinc)
                
                ef(:,1) = ef(:,1) + einc
             END IF

             efr = REAL(ef(:,1)*EXP(-(0,1)*2*pi*REAL(s-1)/nsteps))
             
             efmag(m,n,s) = normr(efr)
             
             IF(efmag(m,n,s)<1d-2) THEN
                EXIT
             END IF

             IF(m==1 .AND. dotr(efr, mesh%faces(n)%n)<0) THEN
                EXIT
             END IF
             
             tan = efr/efmag(m,n,s)
             
             pt(:,m+1,n,s) = pt(:,m,n,s) + tan*dx
          END DO

          IF(npt(n,s)<2) THEN
             npt(n,s) = 0
          END IF
       END DO
    END DO

    CALL save_stream_fields_msh(name, mesh, npt, pt, efmag, scale)

  END SUBROUTINE field_stream

  SUBROUTINE field_domain(name, mesh, scale, nedgestot, x, ga, omega, ri, prd, src, addsrc,&
       origin, dsize, npoints)
    CHARACTER (LEN=*), INTENT(IN) :: name
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: scale, omega
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(srcdata), INTENT(IN) :: src
    LOGICAL, INTENT(IN) :: addsrc
    REAL (KIND=dp), DIMENSION(3), INTENT(IN) :: origin, dsize
    INTEGER, DIMENSION(3), INTENT(IN) :: npoints

    COMPLEX (KIND=dp), DIMENSION(SIZE(x,1),SIZE(x,2),1) :: x2
    CHARACTER (LEN=256) :: oname, numstr
    COMPLEX (KIND=dp), DIMENSION(3,1) :: ef, hf
    COMPLEX (KIND=dp), DIMENSION(3) :: einc, hinc
    COMPLEX (KIND=dp), DIMENSION(3,PRODUCT(npoints)) :: edata, hdata, edata2, hdata2
    REAL (KIND=dp), DIMENSION(3,PRODUCT(npoints)) :: pt, pt2
    REAL (KIND=dp) :: emax
    INTEGER :: n, m, l, index, npoints2

    x2 = RESHAPE(x,(/SIZE(x,1),SIZE(x,2),1/))

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(npoints,pt,origin,dsize,mesh,ga,x2,nedgestot,omega,ri,prd,addsrc,edata,hdata,src)&
    !$OMP PRIVATE(n,m,l,index,ef,hf,einc,hinc)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,npoints(1)
       DO m=1,npoints(2)
          DO l=1,npoints(3)
             index = l + (m-1)*npoints(3) + (n-1)*npoints(3)*npoints(2)

             pt(:,index) = origin - dsize/2 + dsize*(/REAL(n-1)/(npoints(1)-1),&
                  REAL(m-1)/(npoints(2)-1), REAL(l-1)/(npoints(3)-1)/)

             CALL scat_fields(mesh, ga, x2, nedgestot,&
                  omega, ri, prd, pt(:,index), ef, hf)

             edata(:,index) = ef(:,1)
             hdata(:,index) = hf(:,1)

             IF(addsrc) THEN
                CALL src_fields(src, omega, ri, pt(:,index), einc, hinc)

                edata(:,index) = edata(:,index) + einc
                hdata(:,index) = hdata(:,index) + hinc
             END IF
          END DO
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    ! Find maximum field magnitude.
    emax = 0.0_dp

    DO n=1,SIZE(pt,2)
       IF(normc(edata(:,n))>emax) THEN
          emax = normc(edata(:,n))
       END IF
    END DO

    ! Select only field values of sufficient magnitude.
    m = 0

    DO n=1,SIZE(pt,2)
       IF(normc(edata(:,n))>emax*0.01_dp) THEN
          m = m + 1
          pt2(:,m) = pt(:,n)
          edata2(:,m) = edata(:,n)
          hdata2(:,m) = hdata(:,n)
       END IF
    END DO

    CALL save_domain_vector_fields_msh(name, mesh, pt2(:,1:m), edata2(:,1:m), hdata2(:,1:m), scale)

  END SUBROUTINE field_domain

  SUBROUTINE field_mesh(name, mesh, scale, nedgestot, x, ga, omega, ri)
    CHARACTER (LEN=*), INTENT(IN) :: name
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: scale, omega
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    COMPLEX (KIND=dp), INTENT(IN) :: ri

    INTEGER :: n, n2, q, index, nf, nga, na
    COMPLEX (KIND=dp), DIMENSION(3,mesh%nfaces*SIZE(ga)) :: ef, hf
    COMPLEX (KIND=dp), DIMENSION(3) :: et, ht
    COMPLEX (KIND=dp) :: en, hn
    COMPLEX (KIND=dp) :: eps, gae
    REAL (KIND=dp), DIMENSION(3) :: fn, nor
    REAL (KIND=dp) :: detj
    CHARACTER (LEN=256) :: oname, numstr
    TYPE(mesh_container) :: mesh2

    WRITE(*,*) 'Computing near fields on particle mesh.'

    eps = (ri**2)*eps0

    nga = SIZE(ga)

    mesh2%nnodes = mesh%nnodes*nga
    mesh2%nfaces = mesh%nfaces*nga
    ALLOCATE(mesh2%nodes(mesh2%nnodes))
    ALLOCATE(mesh2%faces(mesh2%nfaces))

    DO na=1,nga
       DO n=1,mesh%nnodes
          mesh2%nodes(n + mesh%nnodes*(na-1))%p = MATMUL(ga(na)%j, mesh%nodes(n)%p)
       END DO

       DO n=1,mesh%nfaces
          mesh2%faces(n + mesh%nfaces*(na-1))%node_indices(:) = &
               mesh%faces(n)%node_indices(:) + mesh%nnodes*(na-1)

          mesh2%faces(n + mesh%nfaces*(na-1))%n = MATMUL(ga(na)%j, mesh%faces(n)%n)
       END DO
    END DO

    DO na=1,nga
       detj = ga(na)%detj

       DO n=1,mesh%nfaces

          n2 = n + mesh%nfaces*(na-1)

          et(:) = 0.0_dp
          en = 0.0_dp

          ht(:) = 0.0_dp
          hn = 0.0_dp

          nor = MATMUL(ga(na)%j, mesh%faces(n)%n)
          
          DO nf=1,nga
             gae = ga(na)%ef(nf)

             DO q=1,3
                index = mesh%faces(n)%edge_indices(q)
                index = mesh%edges(index)%parent_index

                fn = MATMUL(ga(na)%j, crossr(mesh%faces(n)%n, rwg(mesh%faces(n)%cp, n, q, mesh)))

                ht = ht - fn*x(index, nf)*gae*detj
                hn = hn + rwgDiv(n, q, mesh)*x(nedgestot + index, nf)*gae*detj
                
                et = et + fn*x(nedgestot + index, nf)*gae
                en = en + rwgDiv(n, q, mesh)*x(index, nf)*gae
             END DO
          END DO

          en = en/((0,1)*omega*eps)
          hn = hn/((0,1)*omega*mu0)

          ef(:,n2) = et + nor*en
          hf(:,n2) = ht + nor*hn
          
       END DO

    END DO

    CALL save_vector_fields_msh(name, mesh2, ef, hf, scale)

    WRITE(*,*) 'Maximum of |E| is ',&
         MAXVAL(SQRT(ABS(ef(1,:))**2 + ABS(ef(2,:))**2 + ABS(ef(3,:))**2))


    !CALL save_field_msh(name, mesh2, en, eta, scale)

    DEALLOCATE(mesh2%nodes, mesh2%faces)

  END SUBROUTINE field_mesh

  SUBROUTINE field_external_mesh(name, mesh, scale, nedgestot, x, ga, omega, ri, prd,&
       addsrc, src, extmesh)
    CHARACTER (LEN=*), INTENT(IN) :: name
    TYPE(mesh_container), INTENT(IN) :: mesh, extmesh
    REAL (KIND=dp), INTENT(IN) :: scale, omega
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    LOGICAL, INTENT(IN) :: addsrc
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(srcdata), INTENT(IN) :: src

    INTEGER :: n, n2, q, index, nf, nga, na
    COMPLEX (KIND=dp), DIMENSION(3,extmesh%nfaces) :: ef, hf
    COMPLEX (KIND=dp), DIMENSION(3) :: einc, hinc
    CHARACTER (LEN=256) :: oname, numstr
    COMPLEX (KIND=dp), DIMENSION(SIZE(x,1),SIZE(x,2),1) :: x2
    REAL (KIND=dp), DIMENSION(3) :: pt

    WRITE(*,*) 'Computing near fields on external mesh.'

    nga = SIZE(ga)

    x2 = RESHAPE(x,(/SIZE(x,1),SIZE(x,2),1/))

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(extmesh,mesh,ga,x2,nedgestot,omega,ri,prd,ef,hf,src,addsrc)&
    !$OMP PRIVATE(n,pt,einc,hinc)
    !$OMP DO SCHEDULE(STATIC)
    DO n=1,extmesh%nfaces
       pt = extmesh%faces(n)%cp

       CALL scat_fields(mesh, ga, x2, nedgestot,&
            omega, ri, prd, pt, ef(:,n), hf(:,n))

       IF(addsrc) THEN
          CALL src_fields(src, omega, ri, pt, einc, hinc)
          
          ef(:,n) = ef(:,n) + einc
          hf(:,n) = hf(:,n) + hinc
       END IF

    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    CALL save_vector_fields_msh(name, extmesh, ef, hf, scale)

  END SUBROUTINE field_external_mesh

  SUBROUTINE gradPnls_mesh(name, mesh, scale, nedgestot, x, ga, omega, ri)
    CHARACTER (LEN=*), INTENT(IN) :: name
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), INTENT(IN) :: scale, omega
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: x
    TYPE(group_action), DIMENSION(:), INTENT(IN) :: ga
    INTEGER, INTENT(IN) :: nedgestot
    COMPLEX (KIND=dp), INTENT(IN) :: ri

    INTEGER :: n, m, p, q, index, nf, nga, na, nbasis
    COMPLEX (KIND=dp), DIMENSION(3,mesh%nfaces*SIZE(ga)) :: gradPn
    COMPLEX (KIND=dp), DIMENSION(mesh%nedges) :: b
    COMPLEX (KIND=dp), DIMENSION(3) :: et, ht
    COMPLEX (KIND=dp) :: en
    COMPLEX (KIND=dp) :: eps
    REAL (KIND=dp), DIMENSION(3) :: fn
    REAL (KIND=dp) :: A, fmDiv
    CHARACTER (LEN=256) :: oname, numstr
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: F

    WRITE(*,*) 'Computing near fields on particle mesh.'

    nbasis = mesh%nedges

    eps = (ri**2)*eps0

    b(:) = 0.0_dp

    DO m=1,mesh%nfaces
       A = mesh%faces(m)%area

       en = 0.0_dp
       
       DO p=1,3
          index = mesh%faces(m)%edge_indices(p)
          index = mesh%edges(index)%parent_index
          
          en = en + rwgDiv(m, p, mesh)*x(index, 1)/((0,1)*omega*eps)
       END DO

       DO q=1,3

          fmDiv = rwgDiv(m, q, mesh)

          index = mesh%faces(m)%edge_indices(q)

          b(index) = b(index) + A*(en**2)*fmDiv
       END DO
    END DO

    ALLOCATE(F(nbasis,nbasis))

    CALL rwg_moments(mesh, F)

    CALL solve_linsys(F, b)

    DEALLOCATE(F)

    DO n=1,mesh%nfaces

       gradPn(:,n) = 0.0_dp
              
       DO q=1,3
          index = mesh%faces(n)%edge_indices(q)
          
          fn = rwg(mesh%faces(n)%cp, n, q, mesh)
          
          gradPn(:,n) = gradPn(:,n) + fn*b(index)
       END DO
    END DO
        
    CALL save_vector_fields_msh(name, mesh, gradPn, gradPn, scale)

  END SUBROUTINE gradPnls_mesh

END MODULE nfpost
