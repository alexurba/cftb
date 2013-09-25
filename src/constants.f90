module constants

  !--------------------------------------------------------------------!
  !         just the library of constants used in the TB tools         !
  !--------------------------------------------------------------------!

  implicit none

  !---------------------- OS dependent constants ----------------------!
  ! CURDIR    : OS dependent symbol for the current directory          !
  ! DIRSEP    : OS dependent symbol for the directory deparator        !
  !--------------------------------------------------------------------!

  character, parameter :: CURDIR = '.'
  character, parameter :: DIRSEP = '/'


  !---------------------- fortran 90 type kinds -----------------------!
  ! DP    : double precision kind                                      !
  !--------------------------------------------------------------------!

  integer,   parameter :: DP = kind(1.0d0)
  

  !----------------- mathematical/numerical constants -----------------!
  ! see below for explanations                                         !
  !--------------------------------------------------------------------!

  double precision, parameter :: PI          = 3.141592653589793d0
  double precision, parameter :: SQRT_PI     = 1.772453850905516d0
  double precision, parameter :: PI_INV      = 1.0d0/PI
  double precision, parameter :: SQRT_PI_INV = 1.0d0/SQRT_PI
  double precision, parameter :: PI2         = 2.0d0*PI
  double precision, parameter :: PI2_INV     = 1.0d0/PI2
  double precision, parameter :: SQRT2       = 1.4142135623731d0
  double precision, parameter :: SQRT2_INV   = 1.0d0/SQRT2

  !------------------------- unit conversion --------------------------!

  double precision, parameter :: eV2Ha       =  0.0367493089d0
  double precision, parameter :: Ha2eV       = 27.2113961d0
  double precision, parameter :: Ry2Ha       =  0.5d0
  double precision, parameter :: Ha2Ry       =  2.0d0
  double precision, parameter :: Ang2Bohr    =  1.88972598501d0

end module constants
