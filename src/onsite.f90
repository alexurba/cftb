program onsite

  !--------------------------------------------------------------------!
  ! onsite -- calculate diagonal elements of the real space TB         !
  !           Hamilton matrix                                          !
  !                                                                    !
  ! Usage: ./onsite.x [options]                                        !
  !                                                                    !
  ! Options:                                                           !
  ! -h  --help     show usage information                              !
  ! -d  --dir      specify directory with tight-binding parameters     !
  !                (default: ./tbparam)                                !
  ! -i  --input    specify input file name (default: INP)              !
  ! -o  --output   specify output file name (default: onsite.out)      !
  ! -f  --format   specify format of the input file (default: mbpp)    !
  ! -p  --param    specify TB parameter set (default: proj01)          !
  !--------------------------------------------------------------------!
  ! 2010-10-29 Alexander Urban (AU)                                    !
  ! 2011-03-28 AU --- updated tbio interface                           !
  !--------------------------------------------------------------------!

  use geometry,  only: geo_init,               &
                       geo_final,              &
                       atomType,               &
                       atomTypeName,           &
                       cooLattice,             &
                       latticeVec,             &
                       nAtoms,                 &
                       nAtomsOfType,           &
                       nTypes

  use pbc,       only: pbc_d_max_origin,       &
                       pbc_number_of_tvecs,    &
                       pbc_compute_tvecs
                 
  use tbio,      only: tbio_init,              &
                       tbio_final,             &
                       tbio_args
                 
  use tbmatrix,  only: mat_init,               &
                       mat_final,              &
                       mat_r_max,              &
                       mat_get_nOrbitals,      &
                       mat_setup_0,            &
                       mat_setup_onsite,       &
                       mat_print_onsite,       &
                       mat_aoNum,              &
                       mat_nOrbsMax,           &
                       mat_H_0, mat_S_0,       &
                       mat_print_matrix

  use tbparam,   only: tbp_init,               &
                       tbp_final,              &
                       tbp_H, tbp_S,           &
                       tbp_O, tbp_pot,         &
                       nl, l_of_il

  implicit none

  !--------------------------------------------------------------------!
  !                             constants                              !
  !                                                                    !
  ! NLMAX     : max. number of different angular momentum channels per !
  !             atom type                                              !
  ! CURDIR    : OS dependent symbol for the current directory          !
  ! DIRSEP    : OS dependent symbol for the directory deparator        !
  !--------------------------------------------------------------------!

  integer,   parameter :: NLMAX  = 3
  character, parameter :: CURDIR = '.'
  character, parameter :: DIRSEP = '/'

  !--------------------------------------------------------------------!
  ! r_max          : max. radial cut-off for all TB matrix elements    !
  ! d_max          : max. cart. distance of any atom from the origin   !
  ! nTvecs         : number of real space translation vectors T within !
  !                  radial cut-off r_max                              !
  ! Tvec(i,ivec)   : i-th component of the ivec-th T vector            !
  ! TvecLen(ivec)  : norm of the ivec-th T vector                      !
  !--------------------------------------------------------------------!
  ! infile         : name of the input file                            !
  ! outfile        : name of the output file                           !
  ! filetype       : type/format of the input file                     !
  ! paramset       : name of the TB parameters set                     !
  ! paramdir       : directory with the TB parameters files            !
  !--------------------------------------------------------------------!
  
  double precision                              :: r_max, d_max
  integer                                       :: nTvecs
  integer,          dimension(:,:), allocatable :: Tvec
  double precision, dimension(:),   allocatable :: TvecLen
  integer                                       :: nev

  character(len=100) :: infile, outfile, filetype, paramset, paramdir
  
  !----------------------------- options ------------------------------!

  call tbio_init(name  = 'onsite.x',                                  &
                 descr = 'only calculate the on-site elements of'     &
                       //' the Hamiltonian',                          & 
                 args  = (/ 'in    ', 'out   ', 'format', 'param ',   &
                            'dir   ' /),                              &
                 extra = '')

  call tbio_args('in',     infile)
  call tbio_args('out',    outfile,  default='onsite.out')
  call tbio_args('format', filetype)
  call tbio_args('param',  paramset)
  call tbio_args('dir',    paramdir)

  call tbio_final()

  !-------------------------- initialization --------------------------!

  call geo_init(infile, filetype)
  call tbp_init(NLMAX, nTypes, atomTypeName, paramset, paramdir)

!  call mat_init(nTypes, nAtoms, atomType, nl, NLMAX, l_of_il, .false., &
!                'real')
  call mat_init(nTypes, nAtoms, atomType, nl, NLMAX, l_of_il, .true., &
                'real')
  nev = mat_get_nOrbitals()

  r_max = mat_r_max
  d_max = pbc_d_max_origin(nAtoms, cooLattice, latticeVec)

  ! translation vectors within the interaction radius:
  nTvecs = pbc_number_of_tvecs(r_max, d_max, latticeVec)
  allocate(Tvec(3,nTvecs), TvecLen(nTvecs))
  call pbc_compute_tvecs(r_max, d_max, latticeVec, nTvecs, Tvec, Tveclen)


  !--------------------------- calculation ----------------------------!

  mat_H_0(:,:) = 0.0d0

  ! add distance dependent on-site contributions:
  call mat_setup_onsite(nTypes, nAtoms, atomType, mat_aoNum,   &
                        cooLattice, latticeVec, nTvecs, Tvec,  &
                        nl, NLMAX, l_of_il, mat_nOrbsMax, nev, &
                        tbp_O, mat_H_0, mat_S_0)
  

  !------------------------------ output ------------------------------!

  call mat_print_onsite(nTypes, nAtoms, atomType, atomTypeName, &
                        nl, NLMAX, l_of_il, trim(outfile))


  !--------------------------- finalization ---------------------------!

  call mat_final()
  call tbp_final()
  call geo_final()

  deallocate(Tvec, TvecLen)

end program onsite
