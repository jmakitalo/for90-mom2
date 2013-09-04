! MODULE: constants
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Various constants for general use.
MODULE constants
  IMPLICIT NONE
  INTRINSIC SQRT
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
  REAL (KIND=dp), PARAMETER :: &
       pi = 3.141592653589, &
       eps0 = 8.8541878176D-12, &
       mu0 = 4*pi*1D-7, &
       c0 = 2.997924580003D8, &
       eta0 = 3.767303134622D2
  REAL (KIND=dp), PARAMETER :: radtodeg = 180.0_dp/pi,&
       degtorad = pi/180.0_dp
  INTEGER, PARAMETER :: prd_none = 1,&
       prd_2d = 2
END MODULE constants
