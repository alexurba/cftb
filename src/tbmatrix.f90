module tbmatrix

  !--------------------------------------------------------------------!
  ! Handling of the tight-binding Hamilton and overlap matrices in     !
  ! various flavours: full, complex matrices, full real matrices,      !
  ! diagonal H matrix, triangular (packed) matrices.                   !
  !   For the set-up of the matrices the real space bond integrals are !
  ! looked up from the bond integral tables and for the expansion of   !
  ! cubic harmonics for the bond geometries.                           ! 
  !   The module additionally provides interfaces to linear algebra    !
  ! libraries (LAPACK) for the solution of the TB eigenvalue problem.  !
  !                                                                    !
  ! Note: the bond integral curves must be initialized before this     !
  !       module !!!                                                   !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! 2010-10-29 Alexander Urban (AU)                                    !
  ! 2011-03-31 AU --- new matrix set-up routines                       !
  !--------------------------------------------------------------------!

  use bondint,   only: bi_bondint,    &
                       bi_onsite,     &
                       bi_get_isInit, &
                       bi_get_l_max,  &
                       bi_get_r_max

  use harmonics, only: rot_init,      &
                       rot_final,     &
                       rot_integral,  &
                       PI, PI2

  use io,        only: io_lower

  implicit none
  save

  public  :: mat_init,            &
             mat_final,           &
             mat_get_nOrbitals,   &
             mat_setup_0,         &
             mat_setup_d,         &
             mat_setup_z,         &
             mat_setup_onsite,    &
             mat_eigen,           &
             mat_eigen_d,         &
             mat_eigen_z,         &
             mat_elem_0,          &
             mat_elem_d,          &
             mat_elem_z,          &
             mat_elem_ons,        &
             mat_print_H,         &
             mat_print_S,         &
             mat_print_onsite,    &
             mat_print_matrix

  private :: mat_bondgeom,                  &
             mat_bondint,                   &
             lapack_solve_gevp_real,        &
             lapack_solve_gevp_real_packed


  !--------------------------- constants ------------------------------!
  ! EPS   : something like the machine epsilon, numerical zero         !
  ! DP    : double precision fortran kind                              !
  ! M_*   : matrix types (see below)                                   !
  !--------------------------------------------------------------------!

  double precision, parameter, private :: EPS = 1.0d-10
  integer,          parameter, private :: DP  = kind(1.0d0)
  integer,          parameter, private :: M_NONE         = 0
  integer,          parameter, private :: M_COMPLEX_FULL = 1
  integer,          parameter, private :: M_REAL_FULL    = 2
  integer,          parameter, private :: M_REAL_DIA     = 3
  
  !----------------------------- public -------------------------------!
  ! nOrbitals  : total number of orbitals in the system, e.g. the size !
  !              of the full H and S matrices is nOrbitals x nOrbitals !
  !----------------------------- private ------------------------------!
  ! mat_isInit : .true., if the module has been initialized            !
  ! mat_has_S  : .true., if the overlap matrix has been initialized    !
  !                                                                    !
  ! mat_mat    : type of the matrix set in memory                      !
  !              = M_NONE         --> no matrices in memory            !
  !                M_COMPLEX_FULL --> complete, complex matrices       !
  !                M_REAL_FULL    --> complete, real mat. (gamma-pnt.) !
  !                M_REAL_DIA     --> diagonal only, real (on-site)    !
  ! mat_M_?    : for M in {H,S}  and ? in {z, d, d_dia}                !
  !              arrays for matrices of different kinds (see above)    !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! mat_aoNum(iat) : AO index for the first orbital of atom iat        !
  ! mat_nOrbsMax   : max number of orbitals per atom                   !
  !--------------------------------------------------------------------!
  ! mat_l_max  : max. angular momentum in the system                   !
  ! mat_r_max  : upper radial cut-off for ALL matrix elements          !
  !--------------------------------------------------------------------!

  integer,                                       public  :: nOrbitals

  logical,                                       private :: mat_isInit
  logical,                                       private :: mat_has_S

  integer,                                       public  :: mat_mat
  double precision, dimension(:,:), allocatable, public  :: mat_H_0
  double precision, dimension(:,:), allocatable, public  :: mat_S_0
  complex(kind=DP), dimension(:,:), allocatable, public  :: mat_H_z
  complex(kind=DP), dimension(:,:), allocatable, public  :: mat_S_z
  double precision, dimension(:,:), allocatable, public  :: mat_H_d
  double precision, dimension(:,:), allocatable, public  :: mat_S_d

  integer,          dimension(:),   allocatable, public  :: mat_aoNum
  integer,          dimension(:),   allocatable, public  :: mat_nOrbs
  integer,                                       public  :: mat_nOrbsMax

  integer,                                       public  :: mat_l_max
  double precision,                              public  :: mat_r_max

  !--------------------------------------------------------------------!
  !        solution of the eigenvalue problem (complex or real)        !
  !--------------------------------------------------------------------!
  
  interface mat_eigen
     module procedure mat_eigen_d, mat_eigen_z
  end interface

  !--------------------------------------------------------------------!
  !                simplistic matrix output procedures                 !
  !--------------------------------------------------------------------!

  interface mat_print_matrix
     module procedure mat_print_matrix_d, mat_print_matrix_z
  end interface

contains

  subroutine mat_init(nTypes, nAtoms, atomType, nl, nlmax, l_of_il, &
                      overlap, matType)

    implicit none

    !------------------------------------------------------------------!
    ! nTypes           : total number of atomic species                !
    ! nAtoms           : total number of atoms                         !
    ! atomType(i)      : atomic species of atom i                      !
    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
    ! nlmax            : max number of angular momentum channels       !
    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
    ! overlap          : .true., if overlap matrix shall be set-up     !
    ! matType          : matrix type;  'complex', 'real'               !
    !------------------------------------------------------------------!

    integer,                                    intent(in) :: nTypes, nAtoms
    integer,          dimension(nAtoms),        intent(in) :: atomType
    integer,          dimension(nTypes),        intent(in) :: nl
    integer,                                    intent(in) :: nlmax
    integer,          dimension(nlmax, nTypes), intent(in) :: l_of_il
    logical,                                    intent(in) :: overlap
    character(len=*),                           intent(in) :: matType

    integer :: iatom, itype

    if (mat_isInit) then
       write(0,*) "Warning: module `tbmatrix' has already been initialized."
       return
    end if

    ! initialize spherical harmonics rotation routines:
    if (bi_get_isInit()) then
       mat_l_max = bi_get_l_max()
       mat_r_max = bi_get_r_max()
    else
       write(0,*) 'Error: bond integral tables not available in ' &
            // "`mat_init()'."
       stop
    end if
    call rot_init(mat_l_max)

    ! matrix type:
    select case(trim(io_lower(matType)))
    case('complex', 'k-point')
       mat_mat = M_COMPLEX_FULL
    case('real', 'gamma')
       mat_mat = M_REAL_FULL
    case default
       mat_mat = M_NONE
       write(0,*) "Error: unknown matrix type in `mat_init()': ", trim(matType)
       stop
    end select

    ! some indices and the dimension of the (full) matrices:
    allocate(mat_aoNum(nAtoms), mat_nOrbs(nTypes))
    mat_nOrbsMax = 0
    do itype = 1, nTypes
       mat_nOrbs(itype) = sum(2*l_of_il(1:nl(itype), itype) + 1)
       mat_nOrbsMax     = max(mat_nOrbsMax, mat_nOrbs(itype))
    end do
    nOrbitals = 0
    do iatom = 1, nAtoms
       mat_aoNum(iatom) = nOrbitals + 1
       itype            = atomType(iatom)
       nOrbitals        = nOrbitals +  mat_nOrbs(itype)
    end do

    ! allocate memory:
    mat_has_S   = overlap
    select case(mat_mat)
    case(M_COMPLEX_FULL)
       if (mat_has_S) then
          allocate( mat_H_0(nOrbitals, nOrbitals), &
                    mat_S_0(nOrbitals, nOrbitals), &
                    mat_H_z(nOrbitals, nOrbitals),  &
                    mat_S_z(nOrbitals, nOrbitals) )
       else
          allocate( mat_H_0(nOrbitals, nOrbitals), &
                    mat_H_z(nOrbitals, nOrbitals)   )
       end if
    case(M_REAL_FULL)
       if (mat_has_S) then
          allocate( mat_H_0(nOrbitals, nOrbitals), &
                    mat_S_0(nOrbitals, nOrbitals), &
                    mat_H_d(nOrbitals, nOrbitals),  &
                    mat_S_d(nOrbitals, nOrbitals) )
       else
          allocate( mat_H_0(nOrbitals, nOrbitals), &
                    mat_H_d(nOrbitals, nOrbitals)   )
       end if
    end select

    mat_isInit = .true.

  end subroutine mat_init

  !--------------------------------------------------------------------!

  subroutine mat_final()

    implicit none

    if (.not. mat_isInit) return

    deallocate(mat_aoNum, mat_nOrbs)
    
    if (allocated(mat_H_z))     deallocate(mat_H_z)
    if (allocated(mat_S_z))     deallocate(mat_S_z)
    if (allocated(mat_H_d))     deallocate(mat_H_d)
    if (allocated(mat_S_d))     deallocate(mat_S_d)
    if (allocated(mat_H_0))     deallocate(mat_H_0)
    if (allocated(mat_S_0))     deallocate(mat_S_0)

    call rot_final()

    mat_mat    = M_NONE
    mat_isInit = .false.

  end subroutine mat_final

  !--------------------------------------------------------------------!
  !                         property requests                          !
  ! See comments at variable declarations for explanation.             !
  !--------------------------------------------------------------------!

  function mat_get_nOrbitals() result(dval)
    implicit none
    double precision :: dval
    dval = nOrbitals
  end function mat_get_nOrbitals


  !====================== matrix set-up routines ======================!


  !--------------------------------------------------------------------!
  !           Hamilton and overlap matrices for T = (0,0,0)            !
  !--------------------------------------------------------------------!

  subroutine mat_setup_0(nTypes, nAtoms, atomType, aoNum, cooLattice,  &
                  latticeVec, nl, nlmax, l_of_il, nloc, norbs, H_0, S_0)

    implicit none

    !------------------------------------------------------------------!
    ! nTypes           : total number of atomic species                !
    ! nAtoms           : total number of atoms                         !
    ! atomType(i)      : atomic species of atom i                      !
    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
    !                    atom j                                        !
    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
    ! nlmax            : max number of angular momentum channels       !
    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
    ! norbs            : number of basis functions (matrix dimension)  !
    ! H_0, S_0         : Hamilton and overlap matrices for T = (0,0,0) !
    !------------------------------------------------------------------!

    integer,                                    intent(in)  :: nTypes, nAtoms
    integer,          dimension(nAtoms),        intent(in)  :: atomType
    integer,          dimension(nAtoms),        intent(in)  :: aoNum
    double precision, dimension(3,nAtoms),      intent(in)  :: cooLattice
    double precision, dimension(3,3),           intent(in)  :: latticeVec
    integer,          dimension(nTypes),        intent(in)  :: nl
    integer,                                    intent(in)  :: nlmax
    integer,          dimension(nlmax, nTypes), intent(in)  :: l_of_il
    integer,                                    intent(in)  :: nloc
    integer,                                    intent(in)  :: norbs
    double precision, dimension(norbs,norbs),   intent(out) :: H_0, S_0

    double precision, dimension(nloc,nloc) :: H_loc, S_loc
    integer,          dimension(nlmax)     :: lvec1, lvec2
    double precision, dimension(3)         :: coo1, coo2

    integer          :: iat1, itype1, nl1, norbs1, iao1, iao1f
    integer          :: iat2, itype2, nl2, norbs2, iao2, iao2f

    H_0(:,:) = 0.0d0
    S_0(:,:) = 0.0d0

    atom1 : do iat1 = 1, nAtoms
       itype1       = atomType(iat1)
       coo1         = cooLattice(1:3,iat1)
       nl1          = nl(itype1)
       lvec1(1:nl1) = l_of_il(1:nl1, itype1)
       norbs1       = sum(2*lvec1(1:nl1) + 1)
       iao1         = aoNum(iat1)
       iao1f        = iao1 + norbs1 - 1

       atom2 : do iat2 = iat1 + 1, nAtoms
          itype2       = atomType(iat2)
          coo2         = cooLattice(1:3,iat2)
          nl2          = nl(itype2)
          lvec2(1:nl2) = l_of_il(1:nl2, itype2)
          norbs2       = sum(2*lvec2(1:nl2) + 1)
          iao2         = aoNum(iat2)
          iao2f        = iao2 + norbs2 - 1
       
          call mat_elem_0(itype1, coo1, nl1, lvec1, &
                          itype2, coo2, nl2, lvec2, &
                          latticeVec, nloc, H_loc, S_loc)

          H_0(iao2:iao2f,iao1:iao1f) = H_loc(1:norbs2,1:norbs1)
          S_0(iao2:iao2f,iao1:iao1f) = S_loc(1:norbs2,1:norbs1)
                
       end do atom2
    end do atom1

    ! symmetrize matrix:
    H_0 = H_0 + transpose(H_0)
    S_0 = S_0 + transpose(S_0)

  end subroutine mat_setup_0

  !--------------------------------------------------------------------!
  !            add distance dependent on-site contributions            !
  !--------------------------------------------------------------------!

  subroutine mat_setup_onsite(nTypes, nAtoms, atomType, aoNum,   &
       cooLattice, latticeVec, nTvecs, Tvec, nl, nlmax, l_of_il, &
       nloc, norbs, has_O, H, S)

    implicit none

    !------------------------------------------------------------------!
    ! nTypes           : total number of atomic species                !
    ! nAtoms           : total number of atoms                         !
    ! atomType(i)      : atomic species of atom i                      !
    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
    !                    atom j                                        !
    ! nTvecs           : number of real space translation vectors T    !
    ! Tvec(i,j)        : i-th component of the j-th T vector           !
    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
    ! nlmax            : max number of angular momentum channels       !
    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
    ! norbs            : number of basis functions (matrix dimension)  !
    ! has_O            : .true. --> distance dependent on-site params  !
    ! H,S              : real valued matrices                          !
    !------------------------------------------------------------------!

    integer,                                    intent(in)    :: nTypes, nAtoms
    integer,          dimension(nAtoms),        intent(in)    :: atomType
    integer,          dimension(nAtoms),        intent(in)    :: aoNum
    double precision, dimension(3,nAtoms),      intent(in)    :: cooLattice
    double precision, dimension(3,3),           intent(in)    :: latticeVec
    integer,                                    intent(in)    :: nTvecs
    integer,          dimension(nTvecs),        intent(in)    :: Tvec
    integer,          dimension(nTypes),        intent(in)    :: nl
    integer,                                    intent(in)    :: nlmax
    integer,          dimension(nlmax, nTypes), intent(in)    :: l_of_il
    integer,                                    intent(in)    :: nloc
    integer,                                    intent(in)    :: norbs
    logical,                                    intent(in)    :: has_O
    double precision, dimension(norbs,norbs),   intent(inout) :: H, S

    double precision, dimension(nloc,nloc) :: H_loc
    integer,          dimension(nlmax)     :: lvec1
    double precision, dimension(3)         :: coo1, coo2

    integer          :: iat1, itype1, nl1, norbs1, iao1, iao1f
    integer          :: iat2, itype2
    integer          :: iorb, il1, m


    atom1 : do iat1 = 1, nAtoms
       itype1        = atomType(iat1)
       coo1          = cooLattice(1:3,iat1)
       nl1           = nl(itype1)
       lvec1(1:nl1)  = l_of_il(1:nl1, itype1)
       norbs1        = sum(2*lvec1(1:nl1) + 1)
       iao1          = aoNum(iat1)
       iao1f         = iao1 + norbs1 - 1

       iorb = iao1
       do il1 = 1, nl1
          do m = -lvec1(il1), lvec1(il1)
             H(iorb,iorb) = H(iorb,iorb) + bi_onsite(lvec1(il1),itype1)
             S(iorb,iorb) = 1.0d0
             iorb = iorb + 1
          end do
       end do

       if (has_O) then
       atom2 : do iat2 = 1, nAtoms
          itype2 = atomType(iat2)
          coo2   = cooLattice(1:3,iat2)
       
          call mat_elem_ons(itype1, coo1, nl1, lvec1, itype2, coo2, &
                            latticeVec, nTvecs, Tvec, nloc, H_loc)

          H(iao1:iao1f,iao1:iao1f) = H(iao1:iao1f,iao1:iao1f) &
                                   + H_loc(1:norbs1,1:norbs1)
                
       end do atom2
       end if

    end do atom1

  end subroutine mat_setup_onsite

  !--------------------------------------------------------------------!
  !        set-up complex Bloch matrix at a particular k-point         !
  !--------------------------------------------------------------------!

  subroutine mat_setup_z(nTypes, nAtoms, atomType, aoNum, cooLattice,  &
       latticeVec, nTvecs, Tvec, nl, nlmax, l_of_il, rkpt, nloc, norbs, &
       H_0, S_0, H, S)

    implicit none

    !------------------------------------------------------------------!
    ! nTypes           : total number of atomic species                !
    ! nAtoms           : total number of atoms                         !
    ! atomType(i)      : atomic species of atom i                      !
    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
    !                    atom j                                        !
    ! nTvecs           : number of real space translation vectors T    !
    ! Tvec(i,j)        : i-th component of the j-th T vector           !
    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
    ! nlmax            : max number of angular momentum channels       !
    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
    ! rkpt(1:3)        : k-point in reciprocal coordinates             !
    ! norbs            : number of basis functions (matrix dimension)  !
    ! H_0, S_0         : real valued matrices for T = (0,0,0)          !
    ! H, S             : complex Bloch matrices (output)               !
    !------------------------------------------------------------------!

    integer,                                    intent(in)  :: nTypes, nAtoms
    integer,          dimension(nAtoms),        intent(in)  :: atomType
    integer,          dimension(nAtoms),        intent(in)  :: aoNum
    double precision, dimension(3,nAtoms),      intent(in)  :: cooLattice
    double precision, dimension(3,3),           intent(in)  :: latticeVec
    integer,                                    intent(in)  :: nTvecs
    integer,          dimension(nTvecs),        intent(in)  :: Tvec
    integer,          dimension(nTypes),        intent(in)  :: nl
    integer,                                    intent(in)  :: nlmax
    integer,          dimension(nlmax, nTypes), intent(in)  :: l_of_il
    double precision, dimension(3),             intent(in)  :: rkpt
    integer,                                    intent(in)  :: nloc
    integer,                                    intent(in)  :: norbs
    double precision, dimension(norbs,norbs),   intent(in)  :: H_0, S_0
    complex(kind=DP), dimension(norbs,norbs),   intent(out) :: H, S

    complex(kind=DP), dimension(nloc,nloc) :: H_loc, S_loc
    integer,          dimension(nlmax)     :: lvec1, lvec2
    double precision, dimension(3)         :: coo1, coo2

    integer          :: iat1, itype1, nl1, norbs1, iao1, iao1f
    integer          :: iat2, itype2, nl2, norbs2, iao2, iao2f

    H(:,:) = (0.0d0, 0.0d0)
    S(:,:) = (0.0d0, 0.0d0)

    atom1 : do iat1 = 1, nAtoms
       itype1        = atomType(iat1)
       coo1          = cooLattice(1:3,iat1)
       nl1           = nl(itype1)
       lvec1(1:nl1)  = l_of_il(1:nl1, itype1)
       norbs1        = sum(2*lvec1(1:nl1) + 1)
       iao1          = aoNum(iat1)
       iao1f         = iao1 + norbs1 - 1

       call mat_elem_z(itype1, coo1, nl1, lvec1, &
                       itype1, coo1, nl1, lvec1, &
                       latticeVec, nTvecs, Tvec, &
                       rkpt, nloc, H_loc, S_loc)
         
       ! the factor 1/2 is for easier symmetrization:
       H(iao1:iao1f,iao1:iao1f) = 0.5d0*H_loc(1:norbs1,1:norbs1)
       S(iao1:iao1f,iao1:iao1f) = 0.5d0*S_loc(1:norbs1,1:norbs1)

       atom2 : do iat2 = iat1 + 1, nAtoms
          itype2       = atomType(iat2)
          coo2         = cooLattice(1:3,iat2)
          nl2          = nl(itype2)
          lvec2(1:nl2) = l_of_il(1:nl2, itype2)
          norbs2       = sum(2*lvec2(1:nl2) + 1)
          iao2         = aoNum(iat2)
          iao2f        = iao2 + norbs2 - 1
       
          call mat_elem_z(itype1, coo1, nl1, lvec1, &
                          itype2, coo2, nl2, lvec2, &
                          latticeVec, nTvecs, Tvec, &
                          rkpt, nloc, H_loc, S_loc)

          H(iao2:iao2f,iao1:iao1f) = H_loc(1:norbs2,1:norbs1)
          S(iao2:iao2f,iao1:iao1f) = S_loc(1:norbs2,1:norbs1)
                
       end do atom2
    end do atom1

    ! Hermitization:
    H = H + conjg(transpose(H))
    S = S + conjg(transpose(S))

    ! add T = (0,0,0) matrices:
    H = H + dcmplx(H_0, 0.0d0)
    S = S + dcmplx(S_0, 0.0d0)

  end subroutine mat_setup_z

  !--------------------------------------------------------------------!
  !        set-up real valued Bloch matrices at the gamma point        !
  !--------------------------------------------------------------------!

  subroutine mat_setup_d(nTypes, nAtoms, atomType, aoNum, cooLattice,  &
            latticeVec, nTvecs, Tvec, nl, nlmax, l_of_il, nloc, norbs, &
            H_0, S_0, H, S)

    implicit none

    !------------------------------------------------------------------!
    ! nTypes           : total number of atomic species                !
    ! nAtoms           : total number of atoms                         !
    ! atomType(i)      : atomic species of atom i                      !
    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
    !                    atom j                                        !
    ! nTvecs           : number of real space translation vectors T    !
    ! Tvec(i,j)        : i-th component of the j-th T vector           !
    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
    ! nlmax            : max number of angular momentum channels       !
    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
    ! rkpt(1:3)        : k-point in reciprocal coordinates             !
    ! norbs            : number of basis functions (matrix dimension)  !
    ! H_0, S_0         : real valued matrices for T = (0,0,0)          !
    ! H, S             : real valued Bloch matrices (output)           !
    !------------------------------------------------------------------!

    integer,                                    intent(in)  :: nTypes, nAtoms
    integer,          dimension(nAtoms),        intent(in)  :: atomType
    integer,          dimension(nAtoms),        intent(in)  :: aoNum
    double precision, dimension(3,nAtoms),      intent(in)  :: cooLattice
    double precision, dimension(3,3),           intent(in)  :: latticeVec
    integer,                                    intent(in)  :: nTvecs
    integer,          dimension(nTvecs),        intent(in)  :: Tvec
    integer,          dimension(nTypes),        intent(in)  :: nl
    integer,                                    intent(in)  :: nlmax
    integer,          dimension(nlmax, nTypes), intent(in)  :: l_of_il
    integer,                                    intent(in)  :: nloc
    integer,                                    intent(in)  :: norbs
    double precision, dimension(norbs,norbs),   intent(in)  :: H_0, S_0
    double precision, dimension(norbs,norbs),   intent(out) :: H, S

    double precision, dimension(nloc,nloc) :: H_loc, S_loc
    integer,          dimension(nlmax)     :: lvec1, lvec2
    double precision, dimension(3)         :: coo1, coo2

    integer          :: iat1, itype1, nl1, norbs1, iao1, iao1f
    integer          :: iat2, itype2, nl2, norbs2, iao2, iao2f

    H(:,:) = 0.0d0
    S(:,:) = 0.0d0

    atom1 : do iat1 = 1, nAtoms
       itype1        = atomType(iat1)
       coo1          = cooLattice(1:3,iat1)
       nl1           = nl(itype1)
       lvec1(1:nl1)  = l_of_il(1:nl1, itype1)
       norbs1        = sum(2*lvec1(1:nl1) + 1)
       iao1          = aoNum(iat1)
       iao1f         = iao1 + norbs1 - 1

       call mat_elem_d(itype1, coo1, nl1, lvec1, &
                       itype1, coo1, nl1, lvec1, &
                       latticeVec, nTvecs, Tvec, &
                       nloc, H_loc, S_loc)
        
       ! the factor 1/2 is for easier symmetrization:
       H(iao1:iao1f,iao1:iao1f) = 0.5d0*H_loc(1:norbs1,1:norbs1)
       S(iao1:iao1f,iao1:iao1f) = 0.5d0*S_loc(1:norbs1,1:norbs1)

       atom2 : do iat2 = iat1 + 1, nAtoms
          itype2       = atomType(iat2)
          coo2         = cooLattice(1:3,iat2)
          nl2          = nl(itype2)
          lvec2(1:nl2) = l_of_il(1:nl2, itype2)
          norbs2       = sum(2*lvec2(1:nl2) + 1)
          iao2         = aoNum(iat2)
          iao2f        = iao2 + norbs2 - 1
       
          call mat_elem_d(itype1, coo1, nl1, lvec1, &
                          itype2, coo2, nl2, lvec2, &
                          latticeVec, nTvecs, Tvec, &
                          nloc, H_loc, S_loc)

          H(iao2:iao2f,iao1:iao1f) = H_loc(1:norbs2,1:norbs1)
          S(iao2:iao2f,iao1:iao1f) = S_loc(1:norbs2,1:norbs1)
                
       end do atom2
    end do atom1

    ! symmetrization:
    H = H + transpose(H)
    S = S + transpose(S)

    ! add T = (0,0,0) matrices:
    H = H + H_0
    S = S + S_0

  end subroutine mat_setup_d



  !============================== I / O ===============================!

  !--------------------------------------------------------------------!
  !     simplistic matrix output procedures (mainly for debugging)     !
  !     --> can be accessed via the `mat_print_matrix()' interface     !
  !--------------------------------------------------------------------!

  subroutine mat_print_matrix_d(A, unit)

    implicit none
    
    double precision, dimension(:,:), intent(in) :: A
    integer,                optional, intent(in) :: unit

    integer :: m, n, u
    integer :: i, j

    if (present(unit)) then
       u = unit
    else
       u = 6
    end if

    m = size(A(:,1))
    n = size(A(1,:))

    do i = 1, m
    do j = 1, n
       write(u,'(1x,F6.3)', advance='no') A(i,j)
    end do
    write(u,*)
    end do
    write(u,*)

  end subroutine mat_print_matrix_d

  !--------------------------------------------------------------------!

  subroutine mat_print_matrix_z(A, unit)

    implicit none
    
    complex(kind=DP), dimension(:,:), intent(in) :: A
    integer,                optional, intent(in) :: unit

    integer :: m, n, u
    integer :: i, j

    if (present(unit)) then
       u = unit
    else
       u = 6
    end if

    m = size(A(:,1))
    n = size(A(1,:))

    do i = 1, m
    do j = 1, n
       write(u,'(1x,F6.3)', advance='no') real(A(i,j))
    end do
    write(u,*)
    end do
    write(u,*)

    do i = 1, m
    do j = 1, n
       write(u,'(1x,F6.3)', advance='no') aimag(A(i,j))
    end do
    write(u,*)
    end do
    write(u,*)

  end subroutine mat_print_matrix_z



  !--------------------------------------------------------------------!
  !                          output routines                           !
  !                                                                    !
  ! mat_print_H()      : write Hamilton matrix to an output file       !
  ! mat_print_S()      : write overlap matrix to an output file        !
  ! mat_print_onsite() : write out onsite elements                     !
  !--------------------------------------------------------------------!

  subroutine mat_print_H(fname)

    implicit none

    character(len=*), intent(in) :: fname

    integer,           parameter :: u_H = 87
    integer                      :: iorb1, iorb2

    if (.not. mat_isInit) then
       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in `mat_print_H()'."
       stop
    else if (mat_mat == M_NONE) then
       write(0,*) "Error: no matrix set up for printing in `mat_print_H()'."
       stop
    end if

    open(u_H, file=trim(fname), status='replace', action='write')

    select case(mat_mat)
    case(M_COMPLEX_FULL)
       write(u_H, '("# real part of the H matrix:")')
       do iorb1 = 1, nOrbitals
          do iorb2 = 1, nOrbitals
             write(u_H, '(1x,F12.8)', advance='no') dble(mat_H_z(iorb2, iorb1))
          end do
          write(u_H, *)
       end do
       write(u_H, *)
       write(u_H, '("# imaginary part of the H matrix:")')
       do iorb1 = 1, nOrbitals
          do iorb2 = 1, nOrbitals
             write(u_H, '(1x,F12.8)', advance='no') aimag(mat_H_z(iorb2, iorb1))
          end do
          write(u_H, *)
       end do
    case(M_REAL_FULL)
       do iorb1 = 1, nOrbitals
          do iorb2 = 1, nOrbitals
             write(u_H, '(1x,F12.8)', advance='no') mat_H_d(iorb2, iorb1)
          end do
          write(u_H, *)
       end do
    end select

    close(u_H)

  end subroutine mat_print_H

  !--------------------------------------------------------------------!

  subroutine mat_print_S(fname)

    implicit none

    character(len=*), intent(in) :: fname

    integer,           parameter :: u_S = 88
    integer                      :: iorb1, iorb2

    if (.not. mat_isInit) then
       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in `mat_print_S()'."
       stop
    else if ((mat_mat == M_NONE) .or. (mat_mat == M_REAL_DIA)) then
       write(0,*) "Error: no matrix set up for printing in `mat_print_S()'."
       stop
    else if (.not. mat_has_S) then
       write(0,*) "Error: overlap matrix not available in `mat_print_S()'."
       stop
    end if

    open(u_S, file=trim(fname), status='replace', action='write')

    select case(mat_mat)
    case(M_COMPLEX_FULL)
       write(u_S, '("# real part of the S matrix:")')
       do iorb1 = 1, nOrbitals
          do iorb2 = 1, nOrbitals
             write(u_S, '(1x,F12.8)', advance='no') dble(mat_S_z(iorb2, iorb1))
          end do
          write(u_S, *)
       end do
       write(u_S, *)
       write(u_S, '("# imaginary part of the S matrix:")')
       do iorb1 = 1, nOrbitals
          do iorb2 = 1, nOrbitals
             write(u_S, '(1x,F12.8)', advance='no') aimag(mat_S_z(iorb2, iorb1))
          end do
          write(u_S, *)
       end do
    case(M_REAL_FULL)
       do iorb1 = 1, nOrbitals
          do iorb2 = 1, nOrbitals
             write(u_S, '(1x,F12.8)', advance='no') mat_S_d(iorb2, iorb1)
          end do
          write(u_S, *)
       end do
    end select

    close(u_S)

  end subroutine mat_print_S

  !--------------------------------------------------------------------!

  subroutine mat_print_onsite(nTypes, nAtoms, atomType, atomTypeName, &
                              nl, nlmax, l_of_il, fname)

    implicit none

    integer,                                    intent(in) :: nTypes, nAtoms
    integer,          dimension(nAtoms),        intent(in) :: atomType
    character(len=2), dimension(nTypes),        intent(in) :: atomTypeName
    integer,          dimension(nTypes),        intent(in) :: nl
    integer,                                    intent(in) :: nlmax
    integer,          dimension(nlmax, nTypes), intent(in) :: l_of_il
    character(len=*),                           intent(in) :: fname

    integer, parameter :: u_ons = 89

    integer          :: iatom1, itype1, il1, iorb1
    integer          :: l1, m1

    if (.not. mat_isInit) then
       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in" &
                  // " `mat_print_onsite()'."
       stop
    end if

    open(u_ons, file=trim(fname), status='replace', action='write')

    iorb1 = 0
    atom1 : do iatom1 = 1, nAtoms
       itype1 = atomType(iatom1)
       write(u_ons, '(1x,I7,1x,A2)', advance='no') iatom1, atomTypeName(itype1)
       angmom1 : do il1 = 1, nl(itype1)
          l1 = l_of_il(il1, itype1)
          mquantnum1 : do m1 = -l1, l1
             iorb1 = iorb1 + 1
             write(u_ons, '(1x,F12.8)', advance='no') mat_H_0(iorb1, iorb1)
          end do mquantnum1
       end do angmom1
       write(u_ons, *)
    end do atom1

    close(u_ons)

  end subroutine mat_print_onsite



  !===================== local atom pair matrices =====================!


  !--------------------------------------------------------------------!
  !                atom pair block of the Bloch matrix                 !
  !--------------------------------------------------------------------!

  subroutine mat_elem_0(itype1, coo1, nl1, l1, itype2, coo2, nl2, l2, &
                        latticeVec, nloc, H_loc, S_loc)

    implicit none

    integer,                                intent(in)  :: itype1, itype2
    double precision, dimension(3),         intent(in)  :: coo1, coo2
    integer,                                intent(in)  :: nl1, nl2
    integer,          dimension(nl1),       intent(in)  :: l1
    integer,          dimension(nl2),       intent(in)  :: l2
    double precision, dimension(3,3),       intent(in)  :: latticeVec
    integer,                                intent(in)  :: nloc
    double precision, dimension(nloc,nloc), intent(out) :: H_loc, S_loc

    double precision :: Hval, Sval
    double precision :: alpha, cosbeta, r
    integer          :: il1, m1, iao1
    integer          :: il2, m2, iao2

    ! T = (0,0,0)

    call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)

    iao1 = 0
    do il1 = 1, nl1
    do m1  = -l1(il1), l1(il1)
       iao1 = iao1 + 1
       iao2 = 0
       do il2 = 1, nl2
       do m2  = -l2(il2), l2(il2)
          iao2 = iao2 + 1

          call mat_bondint(itype1, itype2, l1(il1), l2(il2), m1, m2, &
                           alpha, cosbeta, r, H=Hval, S=Sval)
          H_loc(iao2,iao1) = Hval
          S_loc(iao2,iao1) = Sval

       end do
       end do
    end do
    end do

  end subroutine mat_elem_0

  !--------------------------------------------------------------------!

  subroutine mat_elem_ons(itype1, coo1, nl1, l1, itype2, coo2, &
                          latticeVec, nTvecs, Tvec, nloc, H_loc)

    implicit none

    integer,                                intent(in)  :: itype1, itype2
    double precision, dimension(3),         intent(in)  :: coo1, coo2
    integer,                                intent(in)  :: nl1
    integer,          dimension(nl1),       intent(in)  :: l1
    double precision, dimension(3,3),       intent(in)  :: latticeVec
    integer,                                intent(in)  :: nTvecs
    integer,          dimension(3,nTvecs),  intent(in)  :: Tvec
    integer,                                intent(in)  :: nloc
    double precision, dimension(nloc,nloc), intent(out) :: H_loc

    double precision :: Oval
    double precision :: alpha, cosbeta, r
    integer          :: il1, m1, iao1
    integer          :: il2, m2, iao2
    integer          :: iT, sgn

    ! distance dependent on-site contributions

    H_loc(:,:) = 0.0d0

    call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)

    ononsite : if (r > EPS) then
    iao1 = 0
    do il1 = 1, nl1
    do m1  = -l1(il1), l1(il1)
       iao1 = iao1 + 1
       iao2 = 0
       do il2 = 1, nl1
       do m2  = -l1(il2), l1(il2)
          iao2 = iao2 + 1
          call mat_bondint(itype1, itype2, l1(il1), l1(il2), m1, m2, &
                           alpha, cosbeta, r, O=Oval)
          H_loc(iao2,iao1) = Oval
       end do
       end do
    end do
    end do
    end if ononsite

    ! translations:
    transl : do iT  = 1, nTvecs
    vz     : do sgn = -1, 1, 2
       call mat_bondgeom(coo1, coo2+dble(sgn*Tvec(1:3,iT)), latticeVec, &
                         alpha, cosbeta, r)      
       iao1 = 0
       do il1 = 1, nl1
       do m1  = -l1(il1), l1(il1)
          iao1 = iao1 + 1
          iao2 = 0
          do il2 = 1, nl1
          do m2  = -l1(il2), l1(il2)
             iao2 = iao2 + 1
             call mat_bondint(itype1, itype2, l1(il1), l1(il2), m1, m2, &
                              alpha, cosbeta, r, O=Oval)
             H_loc(iao2,iao1) = H_loc(iao2,iao1) + Oval
          end do
          end do
       end do
       end do
    end do vz
    end do transl

  end subroutine mat_elem_ons

  !--------------------------------------------------------------------!

  subroutine mat_elem_z(itype1, coo1, nl1, l1, itype2, coo2, nl2, l2, &
                        latticeVec, nTvecs, Tvec, rkpt, nloc, H_loc, S_loc)

    implicit none

    integer,                                intent(in)  :: itype1, itype2
    double precision, dimension(3),         intent(in)  :: coo1, coo2
    integer,                                intent(in)  :: nl1, nl2
    integer,          dimension(nl1),       intent(in)  :: l1
    integer,          dimension(nl2),       intent(in)  :: l2
    double precision, dimension(3,3),       intent(in)  :: latticeVec
    integer,                                intent(in)  :: nTvecs
    integer,          dimension(3,nTvecs),  intent(in)  :: Tvec
    double precision, dimension(3),         intent(in)  :: rkpt
    integer,                                intent(in)  :: nloc
    complex(kind=DP), dimension(nloc,nloc), intent(out) :: H_loc, S_loc

    double precision :: alpha, cosbeta, r
    double precision :: KxT
    complex(kind=DP) :: sfac
    double precision :: Hval, Sval
    integer          :: il1, m1, iao1
    integer          :: il2, m2, iao2
    integer          :: iT, sgn

    H_loc(:,:) = (0.0d0, 0.0d0)
    S_loc(:,:) = (0.0d0, 0.0d0)

    ! other T (but only half of the T vectors)

    transl : do iT = 1, nTvecs
    vz     : do sgn = -1, 1, 2
       KxT  = sum(rkpt(1:3)*dble(sgn*Tvec(1:3,iT)))
       sfac = exp(dcmplx(0.0d0, PI2*kxT))
!DEBUG
!       write(55,*) 'kpt           = ', rkpt
!       write(55,*) 'Tvec          = ', iT, sgn*Tvec(1:3,iT)
!       write(55,*) 'sfac/(i*2*pi) = ', KxT
!END DEBUG
       call mat_bondgeom(coo1, coo2+dble(sgn*Tvec(1:3,iT)), latticeVec, &
                         alpha, cosbeta, r)
       iao1 = 0
       do il1 = 1, nl1
       do m1  = -l1(il1), l1(il1)
          iao1 = iao1 + 1
          iao2 = 0
          do il2 = 1, nl2
          do m2  = -l2(il2), l2(il2)
             iao2 = iao2 + 1
             call mat_bondint(itype1, itype2, l1(il1), l2(il2), m1, m2, &
                              alpha, cosbeta, r, H=Hval, S=Sval)
             H_loc(iao2,iao1) = H_loc(iao2,iao1) + sfac*dcmplx(Hval,0.0d0)
             S_loc(iao2,iao1) = S_loc(iao2,iao1) + sfac*dcmplx(Sval,0.0d0)
!DEBUG
!             if (all(rkpt(:) == (/ 0.0d0, 0.25d0, 0.5d0 /))) then
!                write(55,'(1x,6(I2,2x),5(F11.8,2x))')      &
!                     iao1, iao2, l1(il1), m1, l2(il2), m2, &
!                     r, alpha, cosbeta, Hval, Sval
!             end if
!END DEBUG
          end do
          end do
       end do
       end do
    end do vz
    end do transl

!!$    ! Hermitization
!!$    ! also includes the other half of the T vectors: H_ij(+T) = H_ji(-T)
!!$    H_loc = H_loc + conjg(transpose(H_loc))
!!$    S_loc = S_loc + conjg(transpose(S_loc))

  end subroutine mat_elem_z

  !--------------------------------------------------------------------!

  subroutine mat_elem_d(itype1, coo1, nl1, l1, itype2, coo2, nl2, l2, &
                        latticeVec, nTvecs, Tvec, nloc, H_loc, S_loc)

    implicit none

    integer,                                intent(in)  :: itype1, itype2
    double precision, dimension(3),         intent(in)  :: coo1, coo2
    integer,                                intent(in)  :: nl1, nl2
    integer,          dimension(nl1),       intent(in)  :: l1
    integer,          dimension(nl2),       intent(in)  :: l2
    double precision, dimension(3,3),       intent(in)  :: latticeVec
    integer,                                intent(in)  :: nTvecs
    integer,          dimension(3,nTvecs),  intent(in)  :: Tvec
    integer,                                intent(in)  :: nloc
    double precision, dimension(nloc,nloc), intent(out) :: H_loc, S_loc

    double precision :: alpha, cosbeta, r
    double precision :: Hval, Sval
    integer          :: il1, m1, iao1
    integer          :: il2, m2, iao2
    integer          :: iT, sgn

    H_loc(:,:) = (0.0d0, 0.0d0)
    S_loc(:,:) = (0.0d0, 0.0d0)

    ! other T (but only half of the T vectors)

    transl : do iT = 1, nTvecs
    vz     : do sgn = -1, 1, 2
       call mat_bondgeom(coo1, coo2+dble(sgn*Tvec(1:3,iT)), latticeVec, &
                         alpha, cosbeta, r)
       iao1 = 0
       do il1 = 1, nl1
       do m1  = -l1(il1), l1(il1)
          iao1 = iao1 + 1
          iao2 = 0
          do il2 = 1, nl2
          do m2  = -l2(il2), l2(il2)
             iao2 = iao2 + 1
             call mat_bondint(itype1, itype2, l1(il1), l2(il2), m1, m2, &
                              alpha, cosbeta, r, H=Hval, S=Sval)
             H_loc(iao2,iao1) = H_loc(iao2,iao1) + Hval
             S_loc(iao2,iao1) = S_loc(iao2,iao1) + Sval
          end do
          end do
       end do
       end do
    end do vz
    end do transl

!!$    ! Symmetrization
!!$    ! also includes the other half of the T vectors: H_ij(+T) = H_ji(-T)
!!$    H_loc = H_loc + transpose(H_loc)
!!$    S_loc = S_loc + transpose(S_loc)

  end subroutine mat_elem_d

  !--------------------------------------------------------------------!

  subroutine mat_elem_full_z(itype1, coo1, nl1, l1, itype2, coo2, nl2,   &
       l2, latticeVec, nTvecs, Tvec, rkpt, nloc, H_loc, S_loc)

    implicit none

    integer,                                intent(in)  :: itype1, itype2
    double precision, dimension(3),         intent(in)  :: coo1, coo2
    integer,                                intent(in)  :: nl1, nl2
    integer,          dimension(nl1),       intent(in)  :: l1
    integer,          dimension(nl2),       intent(in)  :: l2
    double precision, dimension(3,3),       intent(in)  :: latticeVec
    integer,                                intent(in)  :: nTvecs
    integer,          dimension(3,nTvecs),  intent(in)  :: Tvec
    double precision, dimension(3),         intent(in)  :: rkpt
    integer,                                intent(in)  :: nloc
    complex(kind=DP), dimension(nloc,nloc), intent(out) :: H_loc, S_loc

    double precision                       :: alpha, cosbeta, r
    integer                                :: il1, m1, iao1
    complex(kind=DP), dimension(nloc,nloc) :: H_z, S_z
    double precision, dimension(nloc,nloc) :: H_d, S_d

    H_loc(:,:) = (0.0d0, 0.0d0)
    S_loc(:,:) = (0.0d0, 0.0d0)

    call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)

    if (r < 1.0d-5) then  ! on-site element block

       ! constant diagonal elements:
       iao1 = 0
       do il1 = 1, nl1
       do m1  = -l1(il1), l1(il1)
          iao1 = iao1 + 1
          H_loc(iao1,iao1) = dcmplx(bi_onsite(l1(il1),itype1), 0.0d0)
          S_loc(iao1,iao1) = (1.0d0, 0.0d0)
       end do
       end do

    else ! no on-site block

       ! T = (0,0,0)
       call mat_elem_0(itype1, coo1, nl1, l1,    &
                       itype2, coo2, nl2, l2,    &
                       latticeVec, nloc, H_d, S_d)
       H_loc(:,:) = H_loc(:,:) + dcmplx(H_d(:,:), 0.0d0)
       S_loc(:,:) = S_loc(:,:) + dcmplx(S_d(:,:), 0.0d0)

       ! other T
       call mat_elem_z(itype1, coo1, nl1, l1,    &
                       itype2, coo2, nl2, l2,    &
                       latticeVec, nTvecs, Tvec, &
                       rkpt, nloc, H_z, S_z)
       H_loc(:,:) = H_loc(:,:) + H_z(:,:)
       S_loc(:,:) = S_loc(:,:) + S_z(:,:)

    end if

  end subroutine mat_elem_full_z

  !--------------------------------------------------------------------!

  subroutine mat_elem_full_d(itype1, coo1, nl1, l1, itype2, coo2, nl2,   &
       l2, latticeVec, nTvecs, Tvec, nloc, H_loc, S_loc)

    implicit none

    integer,                                intent(in)  :: itype1, itype2
    double precision, dimension(3),         intent(in)  :: coo1, coo2
    integer,                                intent(in)  :: nl1, nl2
    integer,          dimension(nl1),       intent(in)  :: l1
    integer,          dimension(nl2),       intent(in)  :: l2
    double precision, dimension(3,3),       intent(in)  :: latticeVec
    integer,                                intent(in)  :: nTvecs
    integer,          dimension(3,nTvecs),  intent(in)  :: Tvec
    integer,                                intent(in)  :: nloc
    double precision, dimension(nloc,nloc), intent(out) :: H_loc, S_loc

    double precision                       :: alpha, cosbeta, r
    integer                                :: il1, m1, iao1
    double precision, dimension(nloc,nloc) :: H_d, S_d

    H_loc(:,:) = 0.0d0
    S_loc(:,:) = 0.0d0

    call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)

    if (r < 1.0d-5) then  ! on-site block

       ! constant diagonal elements:
       iao1 = 0
       do il1 = 1, nl1
       do m1  = -l1(il1), l1(il1)
          iao1 = iao1 + 1
          H_loc(iao1,iao1) = bi_onsite(l1(il1),itype1)
          S_loc(iao1,iao1) = 1.0d0
       end do
       end do

    else ! off-diagonal block:

       ! T = (0,0,0)
       call mat_elem_0(itype1, coo1, nl1, l1,    &
                       itype2, coo2, nl2, l2,    &
                       latticeVec, nloc, H_d, S_d)
       H_loc(:,:) = H_loc(:,:) + H_d(:,:)
       S_loc(:,:) = S_loc(:,:) + S_d(:,:)

       ! other T
       call mat_elem_d(itype1, coo1, nl1, l1,    &
                       itype2, coo2, nl2, l2,    &
                       latticeVec, nTvecs, Tvec, &
                       nloc, H_d, S_d)
       H_loc(:,:) = H_loc(:,:) + H_d(:,:)
       S_loc(:,:) = S_loc(:,:) + S_d(:,:)

    end if

  end subroutine mat_elem_full_d




  !================== derivatives of matrix elements ==================!




  subroutine mat_delem_full_z(itype1, coo1, nl1, l1, itype2, coo2, nl2, &
       l2, latticeVec, recLattVec, nTvecs, Tvec, rkpt, nloc, tbp_O,     &
       H_deriv, S_deriv, O1_deriv, O2_deriv)

    implicit none

    integer,                                  intent(in)  :: itype1, itype2
    double precision, dimension(3),           intent(in)  :: coo1, coo2
    integer,                                  intent(in)  :: nl1, nl2
    integer,          dimension(nl1),         intent(in)  :: l1
    integer,          dimension(nl2),         intent(in)  :: l2
    double precision, dimension(3,3),         intent(in)  :: latticeVec
    double precision, dimension(3,3),         intent(in)  :: recLattVec
    integer,                                  intent(in)  :: nTvecs
    integer,          dimension(3,nTvecs),    intent(in)  :: Tvec
    double precision, dimension(3),           intent(in)  :: rkpt
    integer,                                  intent(in)  :: nloc
    logical,                                  intent(in)  :: tbp_O
    complex(kind=DP), dimension(nloc,nloc,3), intent(out) :: H_deriv, S_deriv
    double precision, dimension(nloc,nloc,3), intent(out) :: O1_deriv, O2_deriv

    double precision,            parameter :: displ = 1.0d-6

    double precision, dimension(3)         :: shft
    complex(kind=DP), dimension(nloc,nloc) :: H1, H2, S1, S2
    double precision, dimension(nloc,nloc) :: O1, O2
    integer                                :: i
    complex(kind=DP)                       :: norm

    norm = dcmplx(0.5d0/displ, 0.0d0)

    dim : do i = 1, 3

       shft(1:3) = 0.0d0
       shft(i)   = displ/PI2
       shft(1:3) = shft(1)*recLattVec(1:3,1) &
                 + shft(2)*recLattVec(1:3,2) &
                 + shft(3)*recLattVec(1:3,3)

       call mat_elem_full_z(itype1, coo1-shft, nl1, l1,           &
                            itype2, coo2,      nl2, l2,           &
                            latticeVec, nTvecs, Tvec, rkpt, nloc, &
                            H1, S1)
       call mat_elem_full_z(itype1, coo1+shft, nl1, l1,           &
                            itype2, coo2,      nl2, l2,           &
                            latticeVec, nTvecs, Tvec, rkpt, nloc, &
                            H2, S2)
       H_deriv(:,:,i) = (H2(:,:) - H1(:,:))*norm
       S_deriv(:,:,i) = (S2(:,:) - S1(:,:))*norm

       if (tbp_O) then

          call mat_elem_ons(itype1, coo1-shft, nl1, l1, &
                            itype2, coo2,               &
                            latticeVec, nTvecs, Tvec, nloc, O1)
          call mat_elem_ons(itype1, coo1+shft, nl1, l1, &
                            itype2, coo2,               &
                            latticeVec, nTvecs, Tvec, nloc, O2)
          O1_deriv(:,:,i) = (O2(:,:) - O1(:,:))*norm

          call mat_elem_ons(itype2, coo2,      nl2, l2, &
                            itype1, coo1-shft,          &
                            latticeVec, nTvecs, Tvec, nloc, O1)
          call mat_elem_ons(itype2, coo2,      nl2, l2, &
                            itype1, coo1+shft,          &
                            latticeVec, nTvecs, Tvec, nloc, O2)
          O2_deriv(:,:,i) = (O2(:,:) - O1(:,:))*norm

       end if

    end do dim

  end subroutine mat_delem_full_z

  !--------------------------------------------------------------------!

  subroutine mat_delem_full_d(itype1, coo1, nl1, l1, itype2, coo2, nl2, &
       l2, latticeVec, recLattVec, nTvecs, Tvec, nloc, tbp_O, &
       H_deriv, S_deriv, O1_deriv, O2_deriv)

    implicit none

    integer,                                  intent(in)  :: itype1, itype2
    double precision, dimension(3),           intent(in)  :: coo1, coo2
    integer,                                  intent(in)  :: nl1, nl2
    integer,          dimension(nl1),         intent(in)  :: l1
    integer,          dimension(nl2),         intent(in)  :: l2
    double precision, dimension(3,3),         intent(in)  :: latticeVec
    double precision, dimension(3,3),         intent(in)  :: recLattVec
    integer,                                  intent(in)  :: nTvecs
    integer,          dimension(3,nTvecs),    intent(in)  :: Tvec
    integer,                                  intent(in)  :: nloc
    logical,                                  intent(in)  :: tbp_O
    double precision, dimension(nloc,nloc,3), intent(out) :: H_deriv, S_deriv
    double precision, dimension(nloc,nloc,3), intent(out) :: O1_deriv, O2_deriv

    double precision,            parameter :: displ = 1.0d-6

    double precision, dimension(3)         :: shft
    double precision, dimension(nloc,nloc) :: H1, H2, S1, S2, O1, O2
    integer                                :: i
    double precision                       :: norm

    norm = 0.5d0/displ

    dim : do i = 1, 3
       shft(1:3) = 0.0d0
       shft(i)   = displ/PI2
       shft(1:3) = shft(1)*recLattVec(1:3,1) &
                 + shft(2)*recLattVec(1:3,2) &
                 + shft(3)*recLattVec(1:3,3)

       call mat_elem_full_d(itype1, coo1-shft, nl1, l1,     &
                            itype2, coo2,      nl2, l2,     &
                            latticeVec, nTvecs, Tvec, nloc, &
                            H1, S1)
       call mat_elem_full_d(itype1, coo1+shft, nl1, l1,     &
                            itype2, coo2,      nl2, l2,     &
                            latticeVec, nTvecs, Tvec, nloc, &
                            H2, S2)
       H_deriv(:,:,i) = (H2(:,:) - H1(:,:))*norm
       S_deriv(:,:,i) = (S2(:,:) - S1(:,:))*norm

       if (tbp_O) then

          call mat_elem_ons(itype1, coo1-shft, nl1, l1, &
                            itype2, coo2,               &
                            latticeVec, nTvecs, Tvec, nloc, O1)
          call mat_elem_ons(itype1, coo1+shft, nl1, l1, &
                            itype2, coo2,               &
                            latticeVec, nTvecs, Tvec, nloc, O2)
          O1_deriv(:,:,i) = (O2(:,:) - O1(:,:))*norm

          call mat_elem_ons(itype2, coo2,      nl2, l2, &
                            itype1, coo1-shft,          &
                            latticeVec, nTvecs, Tvec, nloc, O1)
          call mat_elem_ons(itype2, coo2,      nl2, l2, &
                            itype1, coo1+shft,          &
                            latticeVec, nTvecs, Tvec, nloc, O2)
          O2_deriv(:,:,i) = (O2(:,:) - O1(:,:))*norm

       end if

    end do dim

  end subroutine mat_delem_full_d




  !====================== Slater/Koster rotation ======================!



  !--------------------------------------------------------------------!
  !               bond geometry in spherical coordinates               !
  !--------------------------------------------------------------------!

    subroutine mat_bondgeom(coo1, coo2, avec, alpha, cosbeta, r)
    
    implicit none

    !------------------------------------------------------------------!
    ! coo1, coo2  : lattice coordinates of the two atoms               !
    ! avec(i,j)   : i-th ccomponent of the j-th lattice vector         !
    ! r, alpha, cosbeta : spherical coordinates (on exit)              !
    !               note: cosbeta is cos(beta), of course              !
    !------------------------------------------------------------------!

    double precision, dimension(3),   intent(in)  :: coo1, coo2
    double precision, dimension(3,3), intent(in)  :: avec
    double precision,                 intent(out) :: alpha, cosbeta, r

    double precision, dimension(3)                :: vec

    double precision, parameter :: EPS_here = 1.0d-12

    vec(1:3) = coo2(1:3) - coo1(1:3)
    vec(1:3) = vec(1)*avec(1:3,1) &
             + vec(2)*avec(1:3,2) &
             + vec(3)*avec(1:3,3)

    r        = sqrt( sum( vec(1:3)*vec(1:3) ) )
    vec      = vec/r
    cosbeta  = vec(3)
    alpha    = atan2(vec(2),vec(1))

    if (alpha < 0.0d0) then
       alpha = alpha + PI2
    else if (alpha >= PI2) then
       alpha = alpha - PI2
    end if

    if ((abs(alpha) <= EPS_here) .or. (abs(alpha-PI2) <= EPS_here)) then
       alpha = 0.0d0
    end if

  end subroutine mat_bondgeom

  !--------------------------------------------------------------------!
  !                 look up and rotate matrix elements                 !
  !                                                                    !
  ! Matrix elements are integrals of the kind                          !
  !                                                                    !
  !   M(l1,m1,l2,m2) = int { K(r1,l1,m1) Operator K(r2,l2,m2) d_tau }  !
  !                                                                    !
  ! where K(l,m) are cubic harmonics centered at r.                    !
  ! This routine looks up distance (|r2-r1|) dependently tabulated     !
  ! values of matrix elements for the trivial geometries               !
  ! M(r,l1,l2,|m|).  The values for other geometries are calculated    !
  ! via an expansion in rotated cubic harmonics..                      !
  !--------------------------------------------------------------------!

  subroutine mat_bondint(itype1, itype2, l1, l2, m1, m2, alpha, &
                         cosbeta, r, H, S, O)
    
    implicit none

    !------------------------------------------------------------------!
    ! type1, type2  : atomic species numbers                           !
    ! l1, l2        : angular momenta                                  !
    ! m1, m2        : magnetic quantum numbers                         !
    ! r, alpha, cosbeta : spherical coordinates describing the bond    !
    !                 geometry (see `mat_bondgeom')                    !
    ! H, S, O       : values for Hamilton, overlap and on-site curves  !
    !                 optional, on exit only                           !
    !------------------------------------------------------------------!

    integer,                     intent(in)  :: itype1, itype2
    integer,                     intent(in)  :: l1, l2, m1, m2
    double precision,            intent(in)  :: alpha, cosbeta, r
    double precision,  optional, intent(out) :: H, S, O

    integer                                  :: im, n_m
    double precision, dimension(mat_l_max+1) :: c_m
    double precision                         :: Hval, Sval, Oval

    if (r > mat_r_max) then
       if (present(H)) H = 0.0d0
       if (present(S)) S = 0.0d0
       if (present(O)) O = 0.0d0
       return
    end if

    n_m = min(l1, l2) + 1
    call rot_integral(alpha, cosbeta, l1, m1, l2, m2, n_m, c_m)

    if (present(H)) H = 0.0d0
    if (present(S)) S = 0.0d0
    if (present(O)) O = 0.0d0

    do im = 1, n_m
    if (abs(c_m(im)) > EPS) then
       if (present(H) .and. present(S) .and. present(O)) then
          call bi_bondint(r, itype1, itype2, l1, l2, im-1, H=Hval, S=Sval, O=Oval)
          H = H + c_m(im)*Hval
          S = S + c_m(im)*Sval
          O = O + c_m(im)*Oval
       else if (present(H) .and. present(S)) then
          call bi_bondint(r, itype1, itype2, l1, l2, im-1, H=Hval, S=Sval)
          H = H + c_m(im)*Hval
          S = S + c_m(im)*Sval
       else if (present(H) .and. present(O)) then
          call bi_bondint(r, itype1, itype2, l1, l2, im-1, H=Hval, O=Oval)
          H = H + c_m(im)*Hval
          O = O + c_m(im)*Oval
       else if (present(H)) then
          call bi_bondint(r, itype1, itype2, l1, l2, im-1, H=Hval)
          H = H + c_m(im)*Hval
       else if (present(S) .and. present(O)) then
          call bi_bondint(r, itype1, itype2, l1, l2, im-1, S=Sval, O=Oval)
          S = S + c_m(im)*Sval
          O = O + c_m(im)*Oval
       else if (present(S)) then
          call bi_bondint(r, itype1, itype2, l1, l2, im-1, S=Sval)
          S = S + c_m(im)*Sval
       else if (present(O)) then
          call bi_bondint(r, itype1, itype2, l1, l2, im-1, O=Oval)
          O = O + c_m(im)*Oval
       end if
    end if
    end do

  end subroutine mat_bondint





  !====================================================================!
  !                         eigenvalue problem                         !
  !====================================================================!





  subroutine mat_eigen_d(nev, eval, evec)

    implicit none

    !------------------------------------------------------------------!
    ! nev        : the min(nev,nOrbitals) lowest eigenvalues and       !
    !              eigenvectors shall be determined                    !
    ! eval(i)    : i-th eigenvalue (on exit)                           !
    ! evec(i,j)  : i-th component of the j-th eigenvector              !
    !------------------------------------------------------------------!

    integer,                              intent(inout) :: nev
    double precision, dimension(nev),     intent(out)   :: eval
    double precision, dimension(nev,nev), intent(out)   :: evec

    if (.not. mat_isInit) then
       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in `mat_eigen()'."
       stop
    end if

    if (mat_mat == M_NONE) then
       write(0,*) "Error: no suitable Hamilton matrix available in `mat_eigen()'."
       stop
    end if

    nev = min(nev, nOrbitals)

    select case(mat_mat)
    case(M_REAL_FULL)
       if (mat_has_S .and. (nev == nOrbitals)) then
          call lapack_solve_gevp_real(mat_H_d, mat_S_d, nOrbitals, eval, evec)
       end if
    end select

  end subroutine mat_eigen_d

  !--------------------------------------------------------------------!

  subroutine mat_eigen_z(nev, eval, evec)

    implicit none

    !------------------------------------------------------------------!
    ! nev        : the min(nev,nOrbitals) lowest eigenvalues and       !
    !              eigenvectors shall be determined                    !
    ! eval(i)    : i-th eigenvalue (on exit)                           !
    ! evec(i,j)  : i-th component of the j-th eigenvector              !
    !------------------------------------------------------------------!

    integer,                              intent(inout) :: nev
    double precision, dimension(nev),     intent(out)   :: eval
    complex(kind=DP), dimension(nev,nev), intent(out)   :: evec

    if (.not. mat_isInit) then
       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in `mat_eigen()'."
       stop
    end if

    if (mat_mat == M_NONE) then
       write(0,*) "Error: no suitable Hamilton matrix available in `mat_eigen()'."
       stop
    end if

    nev = min(nev, nOrbitals)

    select case(mat_mat)
    case(M_COMPLEX_FULL)
       if (mat_has_S .and. (nev == nOrbitals)) then
          call lapack_solve_gevp_cmplx(mat_H_z, mat_S_z, nOrbitals, eval, evec)
       end if
    end select

  end subroutine mat_eigen_z

  !--------------------------------------------------------------------!
  !                        interface to LAPACK                         !
  !--------------------------------------------------------------------!

  subroutine lapack_solve_gevp_cmplx(H, S, n, eval, evec)

    implicit none

    integer,                          intent(in)     :: n
    complex(kind=DP), dimension(n,n), intent(in)     :: H, S
    double precision, dimension(n),   intent(out)    :: eval
    complex(kind=DP), dimension(n,n), intent(out)    :: evec

    complex(kind=DP), dimension(n,n)              :: B
    complex(kind=DP), dimension(:), allocatable   :: work
    double precision, dimension(:), allocatable   :: rwork
    integer                                       :: info

    complex(kind=DP) :: memsize
    integer          :: lwork 

    external :: ZHEGV

    evec(1:n,1:n) = H(1:n,1:n)
    B(1:n,1:n)    = S(1:n,1:n)

    ! memory size query:
    lwork = -1
    call ZHEGV(1, 'V', 'U', n, evec, n, B, n, eval, memsize, lwork, rwork, info)
    lwork = ceiling(dble(memsize))

    ! solge GEVP:
    allocate(work(lwork), rwork(3*n - 2))
    call ZHEGV(1, 'V', 'U', n, evec, n, B, n, eval, work, lwork, rwork, info)
    deallocate(work, rwork)

    if (info /= 0) then
       write(0,*) 'Error: LAPACK error in ZHEGV: ', info
       stop
    end if

  end subroutine lapack_solve_gevp_cmplx

  !--------------------------------------------------------------------!

  subroutine lapack_solve_gevp_real(H, S, n, eval, evec)

    implicit none

    integer,                          intent(in)     :: n
    double precision, dimension(n,n), intent(in)     :: H, S
    double precision, dimension(n),   intent(out)    :: eval
    double precision, dimension(n,n), intent(out)    :: evec

    double precision, dimension(n,n)              :: B
    double precision, dimension(:), allocatable   :: work
    integer                                       :: info

    double precision :: memsize
    integer :: lwork 

    external :: DSYGV

    evec(1:n,1:n) = H(1:n,1:n)
    B(1:n,1:n)    = S(1:n,1:n)

    ! memory size query:
    lwork = -1
    call DSYGV(1, 'V', 'U', n, evec, n, B, n, eval, memsize, lwork, info)
    lwork = ceiling(memsize)

    ! solge GEVP:
    allocate(work(lwork))
    call DSYGV(1, 'V', 'U', n, evec, n, B, n, eval, work, lwork, info)
    deallocate(work)

    if (info /= 0) then
       write(0,*) 'Error: LAPACK error in DSYGV: ', info
       stop
    end if

  end subroutine lapack_solve_gevp_real

  !--------------------------------------------------------------------!

  subroutine lapack_solve_gevp_real_packed(AP, BP, n, m, val, vec)

    implicit none

    integer,                          intent(in)     :: n, m
    double precision, dimension(n),   intent(inout)  :: AP, BP
    double precision, dimension(m),   intent(out)    :: val
    double precision, dimension(m,m), intent(out)    :: vec

    double precision, dimension(4*(n+m)+1)        :: work
    integer,          dimension(5*m+3)            :: iwork
    integer                                       :: info

    external :: DSPGVD

    call DSPGVD(1, 'V', 'L', m, AP, BP, val, vec, m, work, &
                        4*(n+m)+1, iwork, 5*m+3, info )

    if (info /= 0) then
       write(0,*) 'Error: LAPACK error in DSPGVD: ', info
       stop
    end if

  end subroutine lapack_solve_gevp_real_packed









  !====================================================================!
  !                                                                    !
  !         alternate matrix set-up (taken from the MBPP code)         !
  !                                                                    !
  !   the subroutine mbpp_tblocmat() was implemented by Bernd Meyer    !
  !                                                                    !
  !====================================================================!






  subroutine mat_setup_full_d_alt(nTypes, nAtoms, atomType, aoNum,     &
       cooLattice, latticeVec, nTvecs, Tvec, nl, nlmax, l_of_il, nloc, &
       norbs, H, S)

    implicit none

    !------------------------------------------------------------------!
    ! nTypes           : total number of atomic species                !
    ! nAtoms           : total number of atoms                         !
    ! atomType(i)      : atomic species of atom i                      !
    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
    !                    atom j                                        !
    ! nTvecs           : number of real space translation vectors T    !
    ! Tvec(i,j)        : i-th component of the j-th T vector           !
    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
    ! nlmax            : max number of angular momentum channels       !
    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
    ! rkpt(1:3)        : k-point in reciprocal coordinates             !
    ! norbs            : number of basis functions (matrix dimension)  !
    ! H, S             : real valued Bloch matrices (output)           !
    !------------------------------------------------------------------!

    integer,                                    intent(in)  :: nTypes, nAtoms
    integer,          dimension(nAtoms),        intent(in)  :: atomType
    integer,          dimension(nAtoms),        intent(in)  :: aoNum
    double precision, dimension(3,nAtoms),      intent(in)  :: cooLattice
    double precision, dimension(3,3),           intent(in)  :: latticeVec
    integer,                                    intent(in)  :: nTvecs
    integer,          dimension(nTvecs),        intent(in)  :: Tvec
    integer,          dimension(nTypes),        intent(in)  :: nl
    integer,                                    intent(in)  :: nlmax
    integer,          dimension(nlmax, nTypes), intent(in)  :: l_of_il
    integer,                                    intent(in)  :: nloc
    integer,                                    intent(in)  :: norbs
    double precision, dimension(norbs,norbs),   intent(out) :: H, S

    integer,          dimension(NLMAX) :: l1, l2
    double precision, dimension(9,9)   :: H_loc, S_loc
    double precision, dimension(3)     :: coo1, coo2
    integer :: iat1, nl1, itype1, iao1, iao1f, norbs1, iorb1, iorb2
    integer :: iat2, nl2, itype2, iao2, iao2f, norbs2

    if (nloc > 9) then
       write(0,*) "Error: interface only available for l <= 2"
       stop
    end if

    atom1 : do iat1 = 1, nAtoms
       itype1    = atomType(iat1)
       nl1       = nl(itype1)
       l1(1:nl1) = l_of_il(1:nl1,itype1)
       coo1(:)   = cooLattice(:,iat1)
       iao1      = aoNum(iat1)
       norbs1    = sum(2*l1(1:nl1) + 1)
       iao1f     = iao1 + norbs1 - 1

       atom2 : do iat2 = 1, nAtoms
          itype2    = atomType(iat2)
          nl2       = nl(itype2)
          l2(1:nl2) = l_of_il(1:nl2,itype2)
          coo2(:)   = cooLattice(:,iat2)
          iao2      = aoNum(iat2)
          norbs2    = sum(2*l2(1:nl2) + 1)
          iao2f     = iao2 + norbs2 - 1
       
          call mat_elem_full_d_alt(itype1, coo1, nl1, l1,    &
                                   itype2, coo2, nl2, l2,    &
                                   latticeVec, nTvecs, Tvec, &
                                   9, H_loc, S_loc)
          do iorb1 = 1, norbs1
          do iorb2 = 1, norbs2
             H(iao2+iorb2-1,iao1+iorb1-1) = H_loc(iorb2,iorb1)
             S(iao2+iorb2-1,iao1+iorb1-1) = S_loc(iorb2,iorb1)
          end do
          end do 

       end do atom2
    end do atom1

  end subroutine mat_setup_full_d_alt

  !--------------------------------------------------------------------!

  subroutine mat_elem_full_d_alt(itype1, coo1, nl1, l1, itype2, coo2, nl2, &
       l2, latticeVec, nTvecs, Tvec, nloc, H_loc, S_loc)

    implicit none

    integer,                                intent(in)  :: itype1, itype2
    double precision, dimension(3),         intent(in)  :: coo1, coo2
    integer,                                intent(in)  :: nl1, nl2
    integer,          dimension(nl1),       intent(in)  :: l1
    integer,          dimension(nl2),       intent(in)  :: l2
    double precision, dimension(3,3),       intent(in)  :: latticeVec
    integer,                                intent(in)  :: nTvecs
    integer,          dimension(3,nTvecs),  intent(in)  :: Tvec
    integer,                                intent(in)  :: nloc
    double precision, dimension(nloc,nloc), intent(out) :: H_loc, S_loc

    double precision, dimension(9,9) :: H, S
    double precision, dimension(3)   :: tn1, tn2
    integer                          :: iT, sgn, iao1, iao2
    integer, dimension(9), parameter :: reorder = (/1,4,2,3,9,7,5,6,8/)

    if (nloc > 9) then
       write(0,*) "Error: interface only available for l <= 2"
       stop
    end if

    tn1(1:3) = coo1(1)*latticeVec(1:3,1) &
             + coo1(2)*latticeVec(1:3,2) &
             + coo1(3)*latticeVec(1:3,3)
    tn2(1:3) = coo2(1)*latticeVec(1:3,1) &
             + coo2(2)*latticeVec(1:3,2) &
             + coo2(3)*latticeVec(1:3,3)

    ! T = (0,0,0)
    call mbpp_tblocmat(itype1, tn1, nl1, l1, itype2, tn2, nl2, l2, H, S)
    do iao1 = 1, nloc
    do iao2 = 1, nloc
       H_loc(iao2,iao1) = H(reorder(iao1),reorder(iao2))
       S_loc(iao2,iao1) = S(reorder(iao1),reorder(iao2))
    end do
    end do
    
    do iT = 1, nTvecs
    do sgn = -1, 1, 2
       tn2(:) = coo2(:)+dble(sgn*Tvec(1:3,iT))
       tn2(1:3) = tn2(1)*latticeVec(1:3,1) &
                + tn2(2)*latticeVec(1:3,2) &
                + tn2(3)*latticeVec(1:3,3)
       call mbpp_tblocmat(itype1, tn1, nl1, l1, &
                          itype2, tn2, nl2, l2, H, S)
       do iao1 = 1, nloc
       do iao2 = 1, nloc
          H_loc(iao2,iao1) = H_loc(iao2,iao1) &
                           + H(reorder(iao2),reorder(iao1))
          S_loc(iao2,iao1) = S_loc(iao2,iao1) &
                           + S(reorder(iao2),reorder(iao1))
       end do
       end do
    end do
    end do

  end subroutine mat_elem_full_d_alt

  !--------------------------------------------------------------------!

  subroutine mbpp_tblocmat(ityp1, coo1, nl1, lvec1, &
                           ityp2, coo2, nl2, lvec2, &
                           hlocmat, slocmat)

    implicit none

    integer,                          intent(in)  :: ityp1, ityp2
    double precision, dimension(3),   intent(in)  :: coo1, coo2
    integer,                          intent(in)  :: nl1, nl2
    integer,          dimension(nl1), intent(in)  :: lvec1
    integer,          dimension(nl2), intent(in)  :: lvec2
    double precision, dimension(9,9), intent(out) :: slocmat, hlocmat

    !------------------------------------------------------------------!
    ! Setup the (l1m1,l2m2) block of the real-space Hamilton and       !
    ! overlap matrix for the atom pair (ityp1,iat1)--(ityp2,iat2,tvec).!
    !------------------------------------------------------------------!

    double precision, parameter :: sqrt3o4 = 0.86602540378443864676d0

    integer          :: il1, ipos1, l1, mmax1, m1, norb1
    integer          :: il2, ipos2, l2, mmax2, norb2

    double precision :: costheta, sintheta, cosphi, sinphi, cos2theta, &
                        sin2theta, ctct, stst, cos2phi, sin2phi

    double precision                 :: tnorm, ssk, hsk              
    double precision, dimension(3,3) :: rotp
    double precision, dimension(5,5) :: rotd
    double precision, dimension(9,9) :: sskmat, swork
    double precision, dimension(9,9) :: hskmat, hwork

    double precision, dimension(3)   :: tn
    double precision                 :: dist

    tn   = coo2(:) - coo1(:)
    dist = sqrt(sum(tn*tn))

    ! treat onsite matrix elements separately:
    if (dist < 1.d-2) then
       slocmat(1:9,1:9) = 0.d0
       hlocmat(1:9,1:9) = 0.d0
       ipos1 = 1
       do il1 = 1, nl1
          l1 = lvec1(il1)
          do m1 = 1, 2*l1+1
             slocmat(ipos1,ipos1) = 1.d0
             hlocmat(ipos1,ipos1) = bi_onsite(l1,ityp1)
             ipos1 = ipos1 + 1
          enddo
       enddo
       return
    end if      

    ! two-center integrals:
    tnorm = 1.d0/dist
    tn(1:3) = tnorm*tn(1:3)

    ! set-up rotation matrices so that 'tn' becomes the new z-axis:
    if (tn(3)-1.d0 > -1.d-6) then
       costheta = 1.d0
       sintheta = 0.d0
       cosphi   = 1.d0
       sinphi   = 0.d0
    else if (tn(3)+1.d0 < 1.d-6) then
       costheta = -1.d0
       sintheta = 0.d0
       cosphi   = 1.d0
       sinphi   = 0.d0
    else
       costheta = tn(3)
       sintheta = sqrt(1.d0 - tn(3)*tn(3))
       cosphi   = tn(1)/sintheta
       sinphi   = tn(2)/sintheta
    endif

    ! p rotation matrix.
    rotp(1,1) = costheta
    rotp(1,2) = sintheta*cosphi
    rotp(1,3) = sintheta*sinphi
    rotp(2,1) = -sintheta
    rotp(2,2) = costheta*cosphi
    rotp(2,3) = costheta*sinphi
    rotp(3,1) = 0.d0
    rotp(3,2) = -sinphi
    rotp(3,3) = cosphi

    ! set-up d rotation matrix only when needed:
    if ( (mat_nOrbs(ityp1) >= 5) .or. (mat_nOrbs(ityp2) >= 5) ) then

       cos2theta = costheta*costheta - sintheta*sintheta
       sin2theta = 2.d0*sintheta*costheta
       ctct      = 0.5d0*(costheta*costheta + 1.d0)
       stst      = sintheta*sintheta
       cos2phi   = cosphi*cosphi - sinphi*sinphi
       sin2phi   = 2.d0*sinphi*cosphi

       rotd(1,1) = 3.d0*ctct-2.d0
       rotd(1,2) = sqrt3o4*sin2theta*cosphi
       rotd(1,3) = sqrt3o4*sin2theta*sinphi
       rotd(1,4) = sqrt3o4*stst*cos2phi
       rotd(1,5) = sqrt3o4*stst*sin2phi
       rotd(2,1) = -sqrt3o4*sin2theta
       rotd(2,2) = cos2theta*cosphi
       rotd(2,3) = cos2theta*sinphi
       rotd(2,4) = 0.5d0*sin2theta*cos2phi
       rotd(2,5) = 0.5d0*sin2theta*sin2phi
       rotd(3,1) = 0.d0
       rotd(3,2) = -costheta*sinphi
       rotd(3,3) = costheta*cosphi
       rotd(3,4) = -sintheta*sin2phi
       rotd(3,5) = sintheta*cos2phi
       rotd(4,1) = sqrt3o4*stst
       rotd(4,2) = -0.5d0*sin2theta*cosphi
       rotd(4,3) = -0.5d0*sin2theta*sinphi
       rotd(4,4) = ctct*cos2phi
       rotd(4,5) = ctct*sin2phi
       rotd(5,1) = 0.d0
       rotd(5,2) = sintheta*sinphi
       rotd(5,3) = -sintheta*cosphi
       rotd(5,4) = -costheta*sin2phi
       rotd(5,5) = costheta*cos2phi

    end if

    ! initialize SK representation of local Hamilton and overlap matrix:
    hskmat(1:9,1:9) = 0.D0
    sskmat(1:9,1:9) = 0.D0

    ! take SK parameters from universial (tabulated) curves:
    ipos2 = 1
    il2_loop : do il2 = 1, nl2
       l2    = lvec2(il2)
       mmax2 = 2*l2 + 1

       ipos1 = 1
       il1_loop : do il1 = 1, nl1
          l1    = lvec1(il1)
          mmax1 = 2*l1 + 1

          ! \sigma matrix elements:
          call bi_bondint(dist, ityp1, ityp2, l1, l2, 0, H=hsk, S=ssk)
          sskmat(ipos1,ipos2) = ssk
          hskmat(ipos1,ipos2) = hsk

          ! \pi matrix elements:
          if ( (l1 > 0) .and. (l2 > 0) ) then
             call bi_bondint(dist, ityp1, ityp2, l1, l2, 1, H=hsk, S=ssk)
             sskmat(ipos1+1,ipos2+1) = ssk
             sskmat(ipos1+2,ipos2+2) = ssk
             hskmat(ipos1+1,ipos2+1) = hsk
             hskmat(ipos1+2,ipos2+2) = hsk
          end if

          ! \delta matrix elements:
          if ( (l1 > 1) .and. (l2 > 1) ) then
             call bi_bondint(dist, ityp1, ityp2, l1, l2, 2, H=hsk, S=ssk)
             sskmat(ipos1+3,ipos2+3) = ssk
             sskmat(ipos1+4,ipos2+4) = ssk
             hskmat(ipos1+3,ipos2+3) = hsk
             hskmat(ipos1+4,ipos2+4) = hsk
          end if

          ipos1 = ipos1 + mmax1
       end do il1_loop
         
       ipos2 = ipos2 + mmax2
    end do il2_loop

    ! rotate from Slater-Koster representation to global Carthesian coordinates:

    norb2 = mat_nOrbs(ityp2)
    ipos1 = 1
    do il1 = 1, nl1
       l1 = lvec1(il1)
       mmax1 = 2*l1 + 1
       if (l1 == 0) then
          swork(1,1:norb2) = sskmat(1,1:norb2)
          hwork(1,1:norb2) = hskmat(1,1:norb2)
       else if (l1 == 1) then
          call DGEMM('T','N',mmax1,norb2,mmax1,1.d0,rotp,mmax1, &
                     sskmat(ipos1,1),9,0.d0,swork(ipos1,1),9)
          call DGEMM('T','N',mmax1,norb2,mmax1,1.d0,rotp,mmax1, &
                     hskmat(ipos1,1),9,0.d0,hwork(ipos1,1),9)
       else if (l1 .eq. 2) then
          call DGEMM('T','N',mmax1,norb2,mmax1,1.d0,rotd,mmax1, &
                     sskmat(ipos1,1),9,0.d0,swork(ipos1,1),9)
          call DGEMM('T','N',mmax1,norb2,mmax1,1.d0,rotd,mmax1, &
                     hskmat(ipos1,1),9,0.d0,hwork(ipos1,1),9)
       else
          write(0,*) "Error: l > 2 not available in `tbmatrix'."
          stop
       end if
       ipos1 = ipos1 + mmax1
    end do

    norb1 = mat_nOrbs(ityp1)
    ipos2 = 1
    do il2 = 1, nl2
       l2 = lvec2(il2)
       mmax2 = 2*l2 + 1
       if (l2 == 0) then
          slocmat(1:norb1,1) = swork(1:norb1,1)
          hlocmat(1:norb1,1) = hwork(1:norb1,1)
       else if (l2 == 1) then
          call DGEMM('N','N',norb1,mmax2,mmax2,1.d0,swork(1,ipos2), &
                     9,rotp,mmax2,0.d0,slocmat(1,ipos2),9)
          call DGEMM('N','N',norb1,mmax2,mmax2,1.d0,hwork(1,ipos2), &
                     9,rotp,mmax2,0.d0,hlocmat(1,ipos2),9)
       else if (l2 == 2) then
          call DGEMM('N','N',norb1,mmax2,mmax2,1.d0,swork(1,ipos2), &
                     9,rotd,mmax2,0.d0,slocmat(1,ipos2),9)
          call DGEMM('N','N',norb1,mmax2,mmax2,1.d0,hwork(1,ipos2), &
                     9,rotd,mmax2,0.d0,hlocmat(1,ipos2),9)
       else
          write(0,*) "Error: l > 2 not available in `tbmatrix'."
          stop
       end if
       ipos2 = ipos2 + mmax2
    end do

  end subroutine mbpp_tblocmat







end module tbmatrix






  !======================= DEPRICATED PROCEDURES ======================!





!!$  !--------------------------------------------------------------------!
!!$  !           Hamilton and overlap matrices for T = (0,0,0)            !
!!$  !--------------------------------------------------------------------!
!!$
!!$  subroutine mat_setup_0(nTypes, nAtoms, atomType, aoNum, cooLattice,  &
!!$                  latticeVec, nl, nlmax, l_of_il, nloc, norbs, H_0, S_0)
!!$
!!$    implicit none
!!$
!!$    !------------------------------------------------------------------!
!!$    ! nTypes           : total number of atomic species                !
!!$    ! nAtoms           : total number of atoms                         !
!!$    ! atomType(i)      : atomic species of atom i                      !
!!$    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
!!$    !                    atom j                                        !
!!$    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
!!$    ! nlmax            : max number of angular momentum channels       !
!!$    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
!!$    ! norbs            : number of basis functions (matrix dimension)  !
!!$    ! H_0, S_0         : Hamilton and overlap matrices for T = (0,0,0) !
!!$    !------------------------------------------------------------------!
!!$
!!$    integer,                                    intent(in)  :: nTypes, nAtoms
!!$    integer,          dimension(nAtoms),        intent(in)  :: atomType
!!$    integer,          dimension(nAtoms),        intent(in)  :: aoNum
!!$    double precision, dimension(3,nAtoms),      intent(in)  :: cooLattice
!!$    double precision, dimension(3,3),           intent(in)  :: latticeVec
!!$    integer,          dimension(nTypes),        intent(in)  :: nl
!!$    integer,                                    intent(in)  :: nlmax
!!$    integer,          dimension(nlmax, nTypes), intent(in)  :: l_of_il
!!$    integer,                                    intent(in)  :: nloc
!!$    integer,                                    intent(in)  :: norbs
!!$    double precision, dimension(norbs,norbs),   intent(out) :: H_0, S_0
!!$
!!$    double precision, dimension(nloc,nloc) :: H_loc, S_loc
!!$    integer,          dimension(nlmax)     :: lvec1, lvec2
!!$    double precision, dimension(3)         :: coo1, coo2
!!$
!!$    integer          :: iat1, itype1, nl1, norbs1, iao1, iao1f
!!$    integer          :: iat2, itype2, nl2, norbs2, iao2, iao2f
!!$    integer          :: iorb, il1, m
!!$
!!$    H_0(:,:) = 0.0d0
!!$    S_0(:,:) = 0.0d0
!!$
!!$    atom1 : do iat1 = 1, nAtoms
!!$       itype1       = atomType(iat1)
!!$       coo1         = cooLattice(1:3,iat1)
!!$       nl1          = nl(itype1)
!!$       lvec1(1:nl1) = l_of_il(1:nl1, itype1)
!!$       norbs1       = sum(2*lvec1(1:nl1) + 1)
!!$       iao1         = aoNum(iat1)
!!$       iao1f        = iao1 + norbs1 - 1
!!$
!!$       H_0(iao1:iao1f,iao1:iao1f) = 0.0d0
!!$       S_0(iao1:iao1f,iao1:iao1f) = 0.0d0
!!$
!!$       iorb = iao1
!!$       do il1 = 1, nl1
!!$          do m = -lvec1(il1), lvec1(il1)
!!$             ! the factors 1/2 are for easier symmetrization:
!!$             H_0(iorb,iorb) = 0.5d0*bi_onsite(lvec1(il1),itype1)
!!$             S_0(iorb,iorb) = 0.5d0
!!$             iorb = iorb + 1
!!$          end do
!!$       end do
!!$
!!$       atom2 : do iat2 = iat1 + 1, nAtoms
!!$          itype2       = atomType(iat2)
!!$          coo2         = cooLattice(1:3,iat2)
!!$          nl2          = nl(itype2)
!!$          lvec2(1:nl2) = l_of_il(1:nl2, itype2)
!!$          norbs2       = sum(2*lvec2(1:nl2) + 1)
!!$          iao2         = aoNum(iat2)
!!$          iao2f        = iao2 + norbs2 - 1
!!$       
!!$          call mat_elem_0(itype1, coo1, nl1, lvec1, &
!!$                          itype2, coo2, nl2, lvec2, &
!!$                          latticeVec, nloc, H_loc, S_loc)
!!$
!!$          H_0(iao2:iao2f,iao1:iao1f) = H_loc(1:norbs2,1:norbs1)
!!$          S_0(iao2:iao2f,iao1:iao1f) = S_loc(1:norbs2,1:norbs1)
!!$                
!!$       end do atom2
!!$    end do atom1
!!$
!!$    ! symmetrize matrix:
!!$    H_0 = H_0 + transpose(H_0)
!!$    S_0 = S_0 + transpose(S_0)
!!$
!!$  end subroutine mat_setup_0
!!$
!!$  !--------------------------------------------------------------------!
!!$  !            add distance dependent on-site contributions            !
!!$  !--------------------------------------------------------------------!
!!$
!!$  subroutine mat_setup_onsite(nTypes, nAtoms, atomType, aoNum, cooLattice,  &
!!$                  latticeVec, nTvecs, Tvec, nl, nlmax, l_of_il, nloc, norbs, H)
!!$
!!$    implicit none
!!$
!!$    !------------------------------------------------------------------!
!!$    ! nTypes           : total number of atomic species                !
!!$    ! nAtoms           : total number of atoms                         !
!!$    ! atomType(i)      : atomic species of atom i                      !
!!$    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
!!$    !                    atom j                                        !
!!$    ! nTvecs           : number of real space translation vectors T    !
!!$    ! Tvec(i,j)        : i-th component of the j-th T vector           !
!!$    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
!!$    ! nlmax            : max number of angular momentum channels       !
!!$    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
!!$    ! norbs            : number of basis functions (matrix dimension)  !
!!$    ! H                : real valued Hamilton Matrix                   !
!!$    !------------------------------------------------------------------!
!!$
!!$    integer,                                    intent(in)    :: nTypes, nAtoms
!!$    integer,          dimension(nAtoms),        intent(in)    :: atomType
!!$    integer,          dimension(nAtoms),        intent(in)    :: aoNum
!!$    double precision, dimension(3,nAtoms),      intent(in)    :: cooLattice
!!$    double precision, dimension(3,3),           intent(in)    :: latticeVec
!!$    integer,                                    intent(in)    :: nTvecs
!!$    integer,          dimension(nTvecs),        intent(in)    :: Tvec
!!$    integer,          dimension(nTypes),        intent(in)    :: nl
!!$    integer,                                    intent(in)    :: nlmax
!!$    integer,          dimension(nlmax, nTypes), intent(in)    :: l_of_il
!!$    integer,                                    intent(in)    :: nloc
!!$    integer,                                    intent(in)    :: norbs
!!$    double precision, dimension(norbs,norbs),   intent(inout) :: H
!!$
!!$    double precision, dimension(nloc,nloc) :: H_loc
!!$    integer,          dimension(nlmax)     :: lvec1
!!$    double precision, dimension(3)         :: coo1, coo2
!!$
!!$    integer          :: iat1, itype1, nl1, norbs1, iao1, iao1f
!!$    integer          :: iat2, itype2
!!$
!!$    atom1 : do iat1 = 1, nAtoms
!!$       itype1        = atomType(iat1)
!!$       coo1          = cooLattice(1:3,iat1)
!!$       nl1           = nl(itype1)
!!$       lvec1(1:nl1)  = l_of_il(1:nl1, itype1)
!!$       norbs1        = sum(2*lvec1(1:nl1) + 1)
!!$       iao1          = aoNum(iat1)
!!$       iao1f         = iao1 + norbs1 - 1
!!$
!!$       atom2 : do iat2 = 1, nAtoms
!!$          itype2 = atomType(iat2)
!!$          coo2   = cooLattice(1:3,iat2)
!!$       
!!$          call mat_elem_ons(itype1, coo1, nl1, lvec1, itype2, coo2, &
!!$                            latticeVec, nTvecs, Tvec, nloc, H_loc)
!!$
!!$          H(iao1:iao1f,iao1:iao1f) = H(iao1:iao1f,iao1:iao1f) &
!!$                                   + H_loc(1:norbs1,1:norbs1)
!!$                
!!$       end do atom2
!!$    end do atom1
!!$
!!$  end subroutine mat_setup_onsite





!!$  subroutine old_mat_setup_0(nTypes, nAtoms, atomType, cooLattice,   &
!!$                         latticeVec, nl, nlmax, l_of_il,  norbs, &
!!$                         H_0, S_0)
!!$
!!$    implicit none
!!$
!!$    !------------------------------------------------------------------!
!!$    ! nTypes           : total number of atomic species                !
!!$    ! nAtoms           : total number of atoms                         !
!!$    ! atomType(i)      : atomic species of atom i                      !
!!$    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
!!$    !                    atom j                                        !
!!$    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
!!$    ! nlmax            : max number of angular momentum channels       !
!!$    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
!!$    ! norbs            : number of basis functions (matrix dimension)  !
!!$    ! H_0, S_0       : Hamilton and overlap matrices                 !
!!$    !------------------------------------------------------------------!
!!$
!!$    integer,                                    intent(in)  :: nTypes, nAtoms
!!$    integer,          dimension(nAtoms),        intent(in)  :: atomType
!!$    double precision, dimension(3,nAtoms),      intent(in)  :: cooLattice
!!$    double precision, dimension(3,3),           intent(in)  :: latticeVec
!!$    integer,          dimension(nTypes),        intent(in)  :: nl
!!$    integer,                                    intent(in)  :: nlmax
!!$    integer,          dimension(nlmax, nTypes), intent(in)  :: l_of_il
!!$    integer,                                    intent(in)  :: norbs
!!$    double precision, dimension(norbs,norbs),   intent(out) :: H_0, S_0
!!$
!!$    double precision               :: Hval, Sval
!!$    double precision, dimension(3) :: coo1, coo2
!!$    double precision               :: alpha, cosbeta, r
!!$
!!$    integer          :: iatom1, iatom2, itype1, itype2
!!$    integer          :: il1, il2, iorb1, iorb2
!!$    integer          :: l1, l2, m1, m2
!!$
!!$    if ((.not. mat_isInit)) then
!!$       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in `mat_setup_full()'."
!!$       return
!!$    end if
!!$
!!$    ! calculate real-space Hamilton and Overlap matrices:
!!$    iorb1 = 0
!!$    atom1 : do iatom1 = 1, nAtoms
!!$       itype1 = atomType(iatom1)
!!$       coo1   = cooLattice(1:3, iatom1)
!!$       angmom1 : do il1 = 1, nl(itype1)
!!$          l1 = l_of_il(il1, itype1)
!!$          mquantnum1 : do m1 = -l1, l1
!!$             iorb1 = iorb1 + 1
!!$  
!!$             ! second matrix dimension:
!!$             iorb2 = iorb1
!!$
!!$             ! on-site levels / diagonal elements:
!!$             H_0(iorb1, iorb1) = bi_onsite(l1, itype1)
!!$             S_0(iorb1, iorb1) = 1.0d0
!!$  
!!$             ! same atom, same angula momentum, different magn. quantum number:
!!$             do m2 = m1 + 1, l1
!!$                iorb2 = iorb2 + 1
!!$                H_0(iorb2, iorb1) = 0.0d0
!!$                S_0(iorb2, iorb1) = 0.0d0
!!$             end do
!!$  
!!$             ! same atom, different angula momentum:
!!$             do il2 = il1 + 1, nl(itype1)
!!$                l2 = l_of_il(il2, itype1)
!!$                do m2 = -l2, l2
!!$                   iorb2 = iorb2 + 1
!!$                   H_0(iorb2, iorb1) = 0.0d0
!!$                   S_0(iorb2, iorb1) = 0.0d0
!!$                end do
!!$             end do
!!$     
!!$             ! other atoms in the unit cell:
!!$             do iatom2 = iatom1 + 1, nAtoms
!!$                coo2 = cooLattice(1:3, iatom2)
!!$                call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)
!!$                itype2 = atomType(iatom2)
!!$                do il2 = 1, nl(itype2)
!!$                   l2 = l_of_il(il2, itype2)
!!$                   do m2 = -l2, l2
!!$                      iorb2 = iorb2 + 1
!!$
!!$                      call mat_bondint(itype1, itype2, l1, l2, m1, m2, &
!!$                           alpha, cosbeta, r, H=Hval, S=Sval)
!!$                      H_0(iorb2, iorb1) = Hval
!!$                      S_0(iorb2, iorb1) = Sval
!!$
!!$                   end do
!!$                end do
!!$             end do
!!$ 
!!$          end do mquantnum1
!!$       end do angmom1
!!$    end do atom1
!!$
!!$    ! symmetrization of the matrices:
!!$    do iorb1 = 1, norbs
!!$    do iorb2 = iorb1 + 1, norbs
!!$       ! H matrix:
!!$       H_0(iorb1, iorb2) = H_0(iorb2, iorb1)
!!$       ! S matrix:
!!$       S_0(iorb1, iorb2) = S_0(iorb2, iorb1)
!!$    end do
!!$    end do
!!$
!!$  end subroutine old_mat_setup_0
!!$
!!$  !--------------------------------------------------------------------!
!!$
!!$  subroutine old_mat_setup_z(nTypes, nAtoms, atomType, cooLattice, &
!!$    latticeVec, nTvecs, Tvec, nl, nlmax, l_of_il, rkpt, norbs, &
!!$    H_0, S_0, H_bloch, S_bloch)
!!$
!!$    implicit none
!!$
!!$    !------------------------------------------------------------------!
!!$    ! nTypes           : total number of atomic species                !
!!$    ! nAtoms           : total number of atoms                         !
!!$    ! atomType(i)      : atomic species of atom i                      !
!!$    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
!!$    !                    atom j                                        !
!!$    ! nTvecs           : number of real space translation vectors T    !
!!$    ! Tvec(i,j)        : i-th component of the j-th T vector           !
!!$    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
!!$    ! nlmax            : max number of angular momentum channels       !
!!$    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
!!$    ! rkpt(1:3)        : k-point in reciprocal coordinates             !
!!$    ! norbs            : matrix dimension                              !
!!$    ! H_0, S_0         : real-space matrices                           !
!!$    ! H_bloch, S_bloch : Bloch matrices (set-up here)                  !
!!$    !------------------------------------------------------------------!
!!$
!!$    integer,                                    intent(in)  :: nTypes, nAtoms
!!$    integer,          dimension(nAtoms),        intent(in)  :: atomType
!!$    double precision, dimension(3,nAtoms),      intent(in)  :: cooLattice
!!$    double precision, dimension(3,3),           intent(in)  :: latticeVec
!!$    integer,                                    intent(in)  :: nTvecs
!!$    integer,          dimension(3,nTvecs),      intent(in)  :: Tvec
!!$    integer,          dimension(nTypes),        intent(in)  :: nl
!!$    integer,                                    intent(in)  :: nlmax
!!$    integer,          dimension(nlmax, nTypes), intent(in)  :: l_of_il
!!$    double precision, dimension(3),             intent(in)  :: rkpt
!!$    integer,                                    intent(in)  :: norbs
!!$    double precision, dimension(norbs,norbs),   intent(in)  :: H_0, S_0
!!$    complex(kind=DP), dimension(norbs,norbs),   intent(out) :: H_bloch, S_bloch
!!$
!!$    double precision               :: KxT
!!$    complex(kind=DP)               :: sfac, val_z
!!$    double precision, dimension(3) :: coo1, coo2
!!$    double precision               :: alpha, cosbeta, r
!!$    double precision               :: Hval, Sval
!!$
!!$    integer          :: iatom1, iatom2, itype1, itype2
!!$    integer          :: il1, il2, iorb1, iorb2
!!$    integer          :: ivec
!!$    integer          :: l1, l2, m1, m2
!!$    integer          :: sgn
!!$
!!$    if ((.not. mat_isInit)) then
!!$       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in `mat_setup_full()'."
!!$       return
!!$    end if
!!$
!!$    H_bloch(:,:) = (0.0d0, 0.0d0)
!!$    S_bloch(:,:) = (0.0d0, 0.0d0)
!!$
!!$    ! calculate Bloch Hamiltonian and Overlap
!!$    iorb1 = 0
!!$    atom1 : do iatom1 = 1, nAtoms
!!$       itype1 = atomType(iatom1)
!!$       coo1   = cooLattice(1:3, iatom1)
!!$       angmom1 : do il1 = 1, nl(itype1)
!!$          l1 = l_of_il(il1, itype1)
!!$          mquantnum1 : do m1 = -l1, l1
!!$             iorb1 = iorb1 + 1 
!!$    
!!$             ! atoms in range, but outside of unit cell:
!!$             translations : do ivec = 1, ntvecs
!!$             signum       : do sgn = -1, 1, 2
!!$
!!$                iorb2 = iorb1
!!$
!!$                KxT  = sum(rkpt(1:3)*dble(sgn*Tvec(1:3,ivec)))
!!$                sfac = exp(dcmplx(0.0d0, PI2*kxT))
!!$
!!$                coo2 = cooLattice(1:3, iatom1) + dble(sgn*Tvec(1:3,ivec))
!!$                call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)
!!$
!!$                ! same atom, same l, m:
!!$                call mat_bondint(itype1, itype1, l1, l1, m1, m1, alpha, &
!!$                     cosbeta, r, H=Hval, S=Sval)
!!$                H_bloch(iorb1, iorb1) = H_bloch(iorb1, iorb1) + sfac*dcmplx(Hval, 0.0d0)
!!$                S_bloch(iorb1, iorb1) = S_bloch(iorb1, iorb1) + sfac*dcmplx(Sval, 0.0d0)
!!$  
!!$                ! same atom, same angular momentum, different magn. quantum number:
!!$                do m2 = m1 + 1, l1
!!$                   iorb2 = iorb2 + 1
!!$                   call mat_bondint(itype1, itype1, l1, l1, m1, m2, alpha, &
!!$                        cosbeta, r, H=Hval, S=Sval)
!!$                   H_bloch(iorb2, iorb1) = H_bloch(iorb2, iorb1) + sfac*dcmplx(Hval, 0.0d0)
!!$                   S_bloch(iorb2, iorb1) = S_bloch(iorb2, iorb1) + sfac*dcmplx(Sval, 0.0d0)
!!$                end do
!!$  
!!$                ! same atom, different angular momentum:
!!$                do il2 = il1 + 1, nl(itype1)
!!$                   l2 = l_of_il(il2, itype1)
!!$                   do m2 = -l2, l2
!!$                      iorb2 = iorb2 + 1
!!$                      call mat_bondint(itype1, itype1, l1, l2, m1, m2, alpha, &
!!$                           cosbeta, r, H=Hval, S=Sval)
!!$                      H_bloch(iorb2, iorb1) = H_bloch(iorb2, iorb1) + sfac*dcmplx(Hval, 0.0d0)
!!$                      S_bloch(iorb2, iorb1) = S_bloch(iorb2, iorb1) + sfac*dcmplx(Sval, 0.0d0)
!!$                   end do
!!$                end do
!!$                
!!$                ! other atoms:
!!$                do iatom2 = iatom1 + 1, nAtoms
!!$                   coo2 = cooLattice(1:3, iatom2) + dble(sgn*Tvec(1:3,ivec))
!!$                   call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)
!!$                   itype2 = atomType(iatom2)
!!$                   do il2 = 1, nl(itype2)
!!$                      l2 = l_of_il(il2, itype2)
!!$                      do m2 = -l2, l2
!!$                         iorb2 = iorb2 + 1
!!$                         call mat_bondint(itype1, itype2, l1, l2, m1, m2, alpha, &
!!$                              cosbeta, r, H=Hval, S=Sval)
!!$                         H_bloch(iorb2, iorb1) = H_bloch(iorb2, iorb1) + sfac*dcmplx(Hval, 0.0d0)
!!$                         S_bloch(iorb2, iorb1) = S_bloch(iorb2, iorb1) + sfac*dcmplx(Sval, 0.0d0)
!!$                      end do
!!$                   end do
!!$                end do
!!$
!!$             end do signum
!!$             end do translations
!!$  
!!$          end do mquantnum1
!!$       end do angmom1
!!$    end do atom1
!!$   
!!$    ! hermitization of the matrices:
!!$    do iorb1 = 1, norbs
!!$       H_bloch(iorb1, iorb1) = H_bloch(iorb1, iorb1) &
!!$                             + dcmplx(H_0(iorb1,iorb1), 0.0d0)
!!$       S_bloch(iorb1, iorb1) = S_bloch(iorb1, iorb1) &
!!$                             + dcmplx(S_0(iorb1,iorb1), 0.0d0)
!!$    do iorb2 = iorb1 + 1, norbs
!!$       ! H matrix:
!!$       val_z = H_bloch(iorb2, iorb1) + conjg(H_bloch(iorb1, iorb2)) &
!!$             + dcmplx(H_0(iorb2,iorb1), 0.0d0)
!!$       if (abs(val_z) < EPS) val_z = (0.0d0, 0.0d0)
!!$       H_bloch(iorb2, iorb1) = val_z 
!!$       H_bloch(iorb1, iorb2) = conjg(val_z)
!!$       ! S matrix:
!!$       val_z = S_bloch(iorb2, iorb1) + conjg(S_bloch(iorb1, iorb2)) &
!!$             + dcmplx(S_0(iorb2,iorb1), 0.0d0)
!!$       if (abs(val_z) < EPS) val_z = (0.0d0, 0.0d0)
!!$       S_bloch(iorb2, iorb1) = val_z
!!$       S_bloch(iorb1, iorb2) = conjg(val_z)
!!$    end do
!!$    end do
!!$
!!$  end subroutine old_mat_setup_z
!!$
!!$  !--------------------------------------------------------------------!
!!$
!!$  subroutine old_mat_setup_d(nTypes, nAtoms, atomType, cooLattice, &
!!$    latticeVec, nTvecs, Tvec, nl, nlmax, l_of_il, norbs, &
!!$    H_0, S_0, H_bloch, S_bloch)
!!$
!!$    implicit none
!!$
!!$    !------------------------------------------------------------------!
!!$    ! nTypes           : total number of atomic species                !
!!$    ! nAtoms           : total number of atoms                         !
!!$    ! atomType(i)      : atomic species of atom i                      !
!!$    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
!!$    !                    atom j                                        !
!!$    ! nTvecs           : number of real space translation vectors T    !
!!$    ! Tvec(i,j)        : i-th component of the j-th T vector           !
!!$    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
!!$    ! nlmax            : max number of angular momentum channels       !
!!$    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
!!$    ! norbs            : matrix dimension                              !
!!$    ! H_0, S_0       : real-space matrices                           !
!!$    ! H_bloch, S_bloch : Bloch matrices (set-up here)                  !
!!$    !------------------------------------------------------------------!
!!$
!!$    integer,                                    intent(in)  :: nTypes, nAtoms
!!$    integer,          dimension(nAtoms),        intent(in)  :: atomType
!!$    double precision, dimension(3,nAtoms),      intent(in)  :: cooLattice
!!$    double precision, dimension(3,3),           intent(in)  :: latticeVec
!!$    integer,                                    intent(in)  :: nTvecs
!!$    integer,          dimension(3,nTvecs),      intent(in)  :: Tvec
!!$    integer,          dimension(nTypes),        intent(in)  :: nl
!!$    integer,                                    intent(in)  :: nlmax
!!$    integer,          dimension(nlmax, nTypes), intent(in)  :: l_of_il
!!$    integer,                                    intent(in)  :: norbs
!!$    double precision, dimension(norbs,norbs),   intent(in)  :: H_0, S_0
!!$    double precision, dimension(norbs,norbs),   intent(out) :: H_bloch, S_bloch
!!$
!!$    double precision               :: val_d
!!$    double precision, dimension(3) :: coo1, coo2
!!$    double precision               :: alpha, cosbeta, r
!!$    double precision               :: Hval, Sval
!!$
!!$    integer          :: iatom1, iatom2, itype1, itype2
!!$    integer          :: il1, il2, iorb1, iorb2
!!$    integer          :: ivec
!!$    integer          :: l1, l2, m1, m2
!!$    integer          :: sgn
!!$
!!$    if ((.not. mat_isInit)) then
!!$       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in `mat_setup_full()'."
!!$       return
!!$    end if
!!$
!!$    H_bloch(:,:) = 0.0d0
!!$    S_bloch(:,:) = 0.0d0
!!$
!!$    ! calculate Bloch Hamiltonian and Overlap
!!$    iorb1 = 0
!!$    atom1 : do iatom1 = 1, nAtoms
!!$       itype1 = atomType(iatom1)
!!$       coo1   = cooLattice(1:3, iatom1)
!!$       angmom1 : do il1 = 1, nl(itype1)
!!$          l1 = l_of_il(il1, itype1)
!!$          mquantnum1 : do m1 = -l1, l1
!!$             iorb1 = iorb1 + 1 
!!$    
!!$             ! atoms in range, but outside of unit cell:
!!$             translations : do ivec = 1, ntvecs
!!$             signum       : do sgn = -1, 1, 2
!!$
!!$                iorb2 = iorb1
!!$
!!$                coo2 = cooLattice(1:3, iatom1) + dble(sgn*Tvec(1:3,ivec))
!!$                call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)
!!$
!!$                ! same atom, same l, m:
!!$                call mat_bondint(itype1, itype1, l1, l1, m1, m1, alpha, &
!!$                     cosbeta, r, H=Hval, S=Sval)
!!$                H_bloch(iorb1, iorb1) = H_bloch(iorb1, iorb1) + Hval
!!$                S_bloch(iorb1, iorb1) = S_bloch(iorb1, iorb1) + Sval
!!$  
!!$                ! same atom, same angula momentum, different magn. quantum number:
!!$                do m2 = m1 + 1, l1
!!$                   iorb2 = iorb2 + 1
!!$                   call mat_bondint(itype1, itype1, l1, l1, m1, m2, alpha, &
!!$                        cosbeta, r, H=Hval, S=Sval)
!!$                   H_bloch(iorb2, iorb1) = H_bloch(iorb2, iorb1) + Hval
!!$                   S_bloch(iorb2, iorb1) = S_bloch(iorb2, iorb1) + Sval
!!$                end do
!!$  
!!$                ! same atom, different angula momentum:
!!$                do il2 = il1 + 1, nl(itype1)
!!$                   l2 = l_of_il(il2, itype1)
!!$                   do m2 = -l2, l2
!!$                      iorb2 = iorb2 + 1
!!$                      call mat_bondint(itype1, itype1, l1, l2, m1, m2, alpha, &
!!$                           cosbeta, r, H=Hval, S=Sval)
!!$                      H_bloch(iorb2, iorb1) = H_bloch(iorb2, iorb1) + Hval
!!$                      S_bloch(iorb2, iorb1) = S_bloch(iorb2, iorb1) + Sval
!!$                   end do
!!$                end do
!!$                
!!$                ! other atoms:
!!$                do iatom2 = iatom1 + 1, nAtoms
!!$                   coo2 = cooLattice(1:3, iatom2) + dble(sgn*Tvec(1:3,ivec))
!!$                   call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)
!!$                   itype2 = atomType(iatom2)
!!$                   do il2 = 1, nl(itype2)
!!$                      l2 = l_of_il(il2, itype2)
!!$                      do m2 = -l2, l2
!!$                         iorb2 = iorb2 + 1
!!$                         call mat_bondint(itype1, itype2, l1, l2, m1, m2, alpha, &
!!$                              cosbeta, r, H=Hval, S=Sval)
!!$                         H_bloch(iorb2, iorb1) = H_bloch(iorb2, iorb1) + Hval
!!$                         S_bloch(iorb2, iorb1) = S_bloch(iorb2, iorb1) + Sval
!!$                      end do
!!$                   end do
!!$                end do
!!$
!!$             end do signum
!!$             end do translations
!!$  
!!$          end do mquantnum1
!!$       end do angmom1
!!$    end do atom1
!!$
!!$    ! symmetrization of the matrices:
!!$    do iorb1 = 1, norbs
!!$       H_bloch(iorb1,iorb1) = H_bloch(iorb1,iorb1) + H_0(iorb1,iorb1)
!!$       S_bloch(iorb1,iorb1) = S_bloch(iorb1,iorb1) + S_0(iorb1,iorb1)
!!$    do iorb2 = iorb1 + 1, norbs
!!$       ! H matrix:
!!$       val_d = H_bloch(iorb2,iorb1) + H_bloch(iorb1,iorb2) + H_0(iorb2,iorb1)
!!$       if (abs(val_d) < EPS) val_d = 0.0d0
!!$       H_bloch(iorb2, iorb1) = val_d 
!!$       H_bloch(iorb1, iorb2) = val_d
!!$       ! S matrix:
!!$       val_d = S_bloch(iorb2,iorb1) + S_bloch(iorb1,iorb2) + S_0(iorb2,iorb1)
!!$       if (abs(val_d) < EPS) val_d = 0.0d0
!!$       S_bloch(iorb2, iorb1) = val_d
!!$       S_bloch(iorb1, iorb2) = val_d
!!$    end do
!!$    end do
!!$
!!$  end subroutine old_mat_setup_d
!!$
!!$  !--------------------------------------------------------------------!
!!$
!!$  subroutine old_mat_setup_onsite(nTypes, nAtoms, atomType, cooLattice, &
!!$                              latticeVec, nTvecs, Tvec, nl, nlmax,  &
!!$                              l_of_il, norbs, H)
!!$
!!$    implicit none
!!$
!!$    !------------------------------------------------------------------!
!!$    ! nTypes           : total number of atomic species                !
!!$    ! nAtoms           : total number of atoms                         !
!!$    ! atomType(i)      : atomic species of atom i                      !
!!$    ! cooLattice(i,j)  : i-th component of the lattice coordinates of  !
!!$    !                    atom j                                        !
!!$    ! nTvecs           : number of real space translation vectors T    !
!!$    ! Tvec(i,j)        : i-th component of the j-th T vector           !
!!$    ! nl(il,itype)     : il-th angular momentum channel of atom type j !
!!$    ! nlmax            : max number of angular momentum channels       !
!!$    ! l_of_il(il,itype): angular momentum l of channel il of type itype!
!!$    ! norbs            : dimension of the H matrix                     !
!!$    ! H                : Hamilton matrix                               !
!!$    !------------------------------------------------------------------!
!!$
!!$    integer,                                    intent(in)    :: nTypes, nAtoms
!!$    integer,          dimension(nAtoms),        intent(in)    :: atomType
!!$    double precision, dimension(3,nAtoms),      intent(in)    :: cooLattice
!!$    double precision, dimension(3,3),           intent(in)    :: latticeVec
!!$    integer,                                    intent(in)    :: nTvecs
!!$    integer,          dimension(3,nTvecs),      intent(in)    :: Tvec
!!$    integer,          dimension(nTypes),        intent(in)    :: nl
!!$    integer,                                    intent(in)    :: nlmax
!!$    integer,          dimension(nlmax, nTypes), intent(in)    :: l_of_il
!!$    integer,                                    intent(in)    :: norbs
!!$    double precision, dimension(norbs,norbs),   intent(inout) :: H
!!$
!!$    double precision, dimension(3) :: coo1, coo2
!!$    double precision               :: alpha, cosbeta, r
!!$    double precision               :: Oval
!!$
!!$    integer :: iatom1, iatom2, itype1, itype2, il1, il2, iorb1, iorb2
!!$    integer :: ivec, sgn
!!$    integer :: l1, l2, m1, m2
!!$
!!$    if (.not. mat_isInit) then
!!$       write(0,*) "Error: `tbmatrix' module NOT INITIALIZED in `mat_setup_onsite()'."
!!$       return
!!$    end if
!!$
!!$    iorb1 = 0
!!$    atom1 : do iatom1 = 1, nAtoms
!!$       itype1 = atomType(iatom1)
!!$       coo1   = cooLattice(1:3, iatom1)
!!$       angmom1 : do il1 = 1, nl(itype1)
!!$          l1 = l_of_il(il1, itype1)
!!$          mquantnum1 : do m1 = -l1, l1
!!$             iorb1 = iorb1 + 1
!!$
!!$             ! atoms in the unit cell, e.g. T0 = (0,0,0)
!!$             T0 : do iatom2 = 1, nAtoms
!!$
!!$                if (iatom2 == iatom1) cycle T0
!!$
!!$                itype2 = atomType(iatom2)
!!$                coo2   = cooLattice(1:3, iatom2)
!!$                call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)
!!$
!!$                ! diagonal on-site element:
!!$                iorb2  = iorb1
!!$                call mat_bondint(itype1, itype2, l1, l1, m1, m1, &
!!$                     alpha, cosbeta, r, O=Oval)
!!$                H(iorb1,iorb1) = H(iorb1,iorb1) + Oval
!!$
!!$                ! same angula momentum, different magn. quantum number:
!!$                do m2 = m1 + 1, l1
!!$                   iorb2 = iorb2 + 1
!!$                   call mat_bondint(itype1, itype2, l1, l1, m1, m2, &
!!$                        alpha, cosbeta, r, O=Oval)
!!$                   H(iorb2,iorb1) = H(iorb2,iorb1) + Oval
!!$                   H(iorb1,iorb2) = H(iorb1,iorb2) + Oval
!!$                end do
!!$
!!$                ! different angula momentum:
!!$                do il2 = il1 + 1, nl(itype1)
!!$                l2 = l_of_il(il2, itype1)
!!$                do m2 = -l2, l2
!!$                   iorb2 = iorb2 + 1
!!$                   call mat_bondint(itype1, itype1, l2, l1, m2, m1, &
!!$                        alpha, cosbeta, r, O=Oval)
!!$                   H(iorb2,iorb1) = H(iorb2,iorb1) + Oval
!!$                   H(iorb1,iorb2) = H(iorb1,iorb2) + Oval
!!$                end do
!!$                end do
!!$                
!!$             end do T0
!!$
!!$     
!!$             ! atoms in range, but outside of unit cell, e.g. T /= (0,0,0)
!!$             Ti : do ivec = 1, ntvecs
!!$             vz : do sgn = -1, 1, 2
!!$
!!$                do iatom2 = 1, nAtoms
!!$
!!$                   itype2 = atomType(iatom2)
!!$                   coo2   = cooLattice(1:3, iatom2) + dble(sgn*Tvec(1:3,ivec))
!!$                   call mat_bondgeom(coo1, coo2, latticeVec, alpha, cosbeta, r)
!!$                   if (r > mat_r_max) cycle
!!$                   
!!$                   ! diagonal on-site element:
!!$                   iorb2  = iorb1
!!$                   call mat_bondint(itype1, itype2, l1, l1, m1, m1, &
!!$                        alpha, cosbeta, r, O=Oval)
!!$                   H(iorb1, iorb1) = H(iorb1,iorb1) + Oval
!!$                   
!!$                   ! same angula momentum, different magn. quantum number:
!!$                   do m2 = m1 + 1, l1
!!$                      iorb2 = iorb2 + 1
!!$                      call mat_bondint(itype1, itype2, l1, l1, m1, m2, &
!!$                           alpha, cosbeta, r, O=Oval)
!!$                      H(iorb2, iorb1) = H(iorb2,iorb1) + Oval
!!$                      H(iorb1, iorb2) = H(iorb1,iorb2) + Oval
!!$                   end do
!!$                   
!!$                   ! different angula momentum:
!!$                   do il2 = il1 + 1, nl(itype1)
!!$                   l2 = l_of_il(il2, itype1)
!!$                   do m2 = -l2, l2
!!$                      iorb2 = iorb2 + 1
!!$                      call mat_bondint(itype1, itype2, l1, l2, m1, m2, &
!!$                           alpha, cosbeta, r, O=Oval)
!!$                      H(iorb2, iorb1) = H(iorb2,iorb1) + Oval
!!$                      H(iorb1, iorb2) = H(iorb1,iorb2) + Oval
!!$                   end do
!!$                   end do
!!$
!!$                end do               
!!$
!!$             end do vz
!!$             end do Ti
!!$
!!$  
!!$          end do mquantnum1
!!$       end do angmom1
!!$    end do atom1
!!$
!!$  end subroutine old_mat_setup_onsite




!!$  !--------------------------------------------------------------------!
!!$  !        fill matrices with values for one set of iorb1/iorb2        !
!!$  !--------------------------------------------------------------------!
!!$
!!$  subroutine mat_element(iorb1, iorb2, itype1, itype2, l1, l2, m1, m2, &
!!$                         alpha, cosbeta, r, sfac)
!!$
!!$    implicit none
!!$
!!$    integer,                          intent(in)    :: iorb1, iorb2
!!$    integer,                          intent(in)    :: itype1, itype2
!!$    integer,                          intent(in)    :: l1, l2, m1, m2
!!$    double precision,                 intent(in)    :: alpha, cosbeta, r
!!$    complex(kind=DP),                 intent(in)    :: sfac
!!$
!!$    double precision :: Hval, Sval
!!$  
!!$
!!$    if (mat_has_S) then
!!$       call mat_bondint(itype1, itype2, l1, l2, m1, m2, alpha, cosbeta, r, &
!!$                        H=Hval, S=Sval)
!!$       select case(mat_mat)
!!$       case(M_COMPLEX_FULL)
!!$          mat_H_z(iorb2, iorb1) = mat_H_z(iorb2, iorb1) + sfac*dcmplx(Hval, 0.0d0)
!!$          mat_S_z(iorb2, iorb1) = mat_S_z(iorb2, iorb1) + sfac*dcmplx(Sval, 0.0d0)
!!$       case(M_REAL_FULL)
!!$          mat_H_d(iorb2, iorb1) = mat_H_d(iorb2, iorb1) + Hval
!!$          mat_S_d(iorb2, iorb1) = mat_S_d(iorb2, iorb1) + Sval
!!$       end select
!!$    else
!!$       call mat_bondint(itype1, itype2, l1, l2, m1, m2, alpha, cosbeta, r, &
!!$                        H=Hval)                         
!!$       select case(mat_mat)
!!$       case(M_COMPLEX_FULL)
!!$          mat_H_z(iorb2, iorb1) = mat_H_z(iorb2, iorb1) + sfac*dcmplx(Hval, 0.0d0)
!!$       case(M_REAL_FULL)
!!$          mat_H_d(iorb2, iorb1) = mat_H_d(iorb2, iorb1) + Hval
!!$       end select
!!$    end if
!!$    
!!$  end subroutine mat_element
!!$
!!$  !--------------------------------------------------------------------!
!!$  ! Look up distance dependent on-site contributions and update the    !
!!$  ! corresponding matrix elements.                                     !
!!$  !--------------------------------------------------------------------!
!!$
!!$  subroutine mat_update_onsite(iorb1, iorb2, itype1, itype2, l1, l2, &
!!$                               m1, m2, alpha, cosbeta, r)
!!$
!!$    implicit none
!!$
!!$    integer,          intent(in)    :: iorb1, iorb2
!!$    integer,          intent(in)    :: itype1, itype2
!!$    integer,          intent(in)    :: l1, l2, m1, m2
!!$    double precision, intent(in)    :: alpha, cosbeta, r
!!$
!!$    double precision :: Oval
!!$
!!$    call mat_bondint(itype1, itype2, l1, l2, m1, m2, alpha, cosbeta, r, O=Oval)
!!$
!!$    select case(mat_mat)
!!$    case(M_COMPLEX_FULL)
!!$       mat_H_z(iorb2, iorb1) = mat_H_z(iorb2,iorb1) + dcmplx(Oval, 0.0d0)
!!$    case(M_REAL_FULL)
!!$       mat_H_d(iorb2, iorb1) = mat_H_d(iorb2,iorb1) + Oval
!!$    end select
!!$
!!$    if (iorb1 /= iorb2) then
!!$       call mat_bondint(itype1, itype2, l2, l1, m2, m1, alpha, cosbeta, r, O=Oval)
!!$       select case(mat_mat)
!!$       case(M_COMPLEX_FULL)
!!$          mat_H_z(iorb1, iorb2) = mat_H_z(iorb1,iorb2) + dcmplx(Oval, 0.0d0)
!!$       case(M_REAL_FULL)
!!$          mat_H_d(iorb1, iorb2) = mat_H_d(iorb1,iorb2) + Oval
!!$       end select
!!$    end if
!!$
!!$  end subroutine mat_update_onsite

