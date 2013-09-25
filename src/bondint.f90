module bondint

  !--------------------------------------------------------------------!
  ! This module contains the paramterization of the bond integrals (BI)!
  ! for tight-binding (TB) calculations.  BIs are integrals over cubic !
  ! harmonics (i.e. atomic-like orbitals) that are functions of the    !
  ! angular momentum l and the magnetic quantum number m:              !
  !                                                                    !
  ! cubic harmonics K(l,m)          bondint(l1,l2,|m|)                 !
  ! ======================          = int{ K(l1,m) Oper. K(l2,m) }     !
  !                                                                    !
  !   l   m   orbital               l1   l2   |m|   bond               !
  !   -------------------           ------------------------           !
  !   0   0   s                     0    0     0    ss sigma           !
  !   1   0   p_z                   0    1     0    sp sigma           !
  !   1   1   p_x                   0    2     0    sd sigma           !
  !   1  -1   p_y                   1    0     0    ps sigma           !
  !   2   0   d_{z^2}               1    1     0    pp sigma           !
  !   2   1   d_{xz}                1    1     1    pp pi              !
  !   2  -1   d_{yz}                1    2     0    pd sigma           !
  !   2   2   d_{x^2-y^2}           1    2     1    pd pi              !
  !   2  -2   d_{xy}                2    0     0    ds sigma           !
  !                                 2    1     0    dp sigma           !
  !                                 2    1     1    dp pi              !
  !                                 2    2     0    dd sigma           !
  !                                 2    2     1    dd pi              !
  !                                 2    2     2    dd delta           !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! 2010-10-17 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use interpol, only: deriv2,     &
                      cs, cs_d1, cs_d2

  use io,       only: io_adjustl, &
                      io_lower,   &
                      io_upper

  implicit none
  save

  public  :: bi_init,              &
             bi_final,             &
             bi_bondint,           &
             bi_onsite,            &
             bi_pairpot,           &
             bi_tabulate,          &
             bi_set_onsitelevel,   &
             bi_get_isInit,        &
             bi_get_hasHmat,       &
             bi_get_hasSmat,       &
             bi_get_hasOmat,       &
             bi_get_hasV,          &
             bi_get_r_min,         &
             bi_get_r_max,         &
             bi_get_r_min_pot,     &
             bi_get_r_max_pot,     &
             bi_get_l_min,         &
             bi_get_l_max,         &
             bi_get_nbi,           &
             bi_get_bi,            &
             polynomial,           &
             chebypol,             &
             dchebypol,            &
             polyexp,              &
             polyexp_cut
             
             

  private :: bi_memory_index,      &
             bi_pot_idx,           &
             bi_tabulate_mat,      &
             bi_tabulate_pot,      &
             generate_bondint_IDs, &
             tab_chebychev,        &
             tab_polyexp
             
  !---------------------------- constants -----------------------------!
  ! NPNTS   : number of nodes for the spline interpolation             !
  !--------------------------------------------------------------------!

  integer,          parameter, private :: NPNTS = 50
  double precision, parameter, private :: EPS = 1.0d-10

  !--------------------------------------------------------------------!
  ! bi_isInit     : .true., if the module has been initialized         !
  !                                                                    !
  ! ibMat         : next index for matrix elements initialization      !
  ! nBonds        : number of different interactions                   !
  ! memSize       : dimension of the matrix elements arrays            !
  ! memSizePot    : dimension of the pair potential array              !
  !                                                                    !
  ! bi(l1,l2,|m|) : unique bond integral ID for angular momenta l1, l2 !
  !                 and unsigned value of the magn. quantum number m   !
  !                 note: |m| <= min(l1,l2)                            !
  ! nbi           : total number of those bond integral IDs            !
  ! biIdx(ibi, itype1, itype2)                                         !
  !               : matrix element index for the ibi-th bond between   !
  !                 types itype1 and itype2; ibi = bi(l1,l2,|m|)       !
  !                                                                    !
  ! bi_ntypes     : total number of atomic species                     !
  ! bi_l_min      : lowest appearing angular momentum                  !
  ! bi_l_max      : highest appearing angular momentum                 !
  ! bi_r_min      : lower radial cut-off for bond integral curves      !
  ! bi_r_max      : upper radial cut-off                               !
  !                                                                    !
  ! bi_hasHmat    : .true. if Hamilton matrix parameters are provided  !
  ! bi_hasSmat    : .true. if overlap matrix parameters are provided   !
  ! bi_hasOmat    : .true. if distance dependent on-site parameters    !
  !                 are provided                                       !
  ! bi_hasV       : .true., if pair potential parameters are available !
  !                                                                    !
  ! onsite_level(l, itype) : on-site level for angular momentum l and  !
  !                          atomic species number itype               !
  !                                                                    !
  ! R_mat(i)      : r values for the tabulated range for the matrix    !
  !                 elements                                           !
  ! H_mat(idx)    : Hamilton matrix elements parameters                !
  ! S_mat(idx)    : overlap/screening matrix elements parameters       !
  ! d2?_mat(idx)  : second derivatives of the matrix elements for      !
  !                 ? in {H, S, O}                                     !
  !                                                                    !
  ! R_pot(i)      : r values for the tabulated range of the pair pot.  !
  ! V_pot(idxv)   : pair potential parameters                          !
  !                                                                    !
  ! Storage of the parameter cvurves                                   !
  ! ================================                                   !
  !                                                                    !
  ! The parameter curves for the matrices (H, S or O) are interpolated !
  ! with cubic splines on NPNTS nodes.  Values for the NPNTS nodes are !
  ! stored in the 1-d arrays ?_mat and the 2nd derivatives (necessary  !
  ! for the calculation of the cubic splines) are stored in d2?_mat.   !
  !   To access parameters for the bond between elements A and B for   !
  ! angular momenta l1 and l2 and a given unsigned value of the magn.  !
  ! quantum number m:                                                  !
  !                                                                    !
  !   (1)  get the ID for this bond kind:                              !
  !        --> ibi = bi(l1, l2, |m|)                                   !
  !   (2)  get index for that bond kind between elements A and B       !
  !        --> idx = biIdx(ibi, itypeA, itypeB)                        !
  !        where itypeA and itypeB are the atomic species indices of   !
  !        A and B.                                                    !
  !   (3)  access the parameters:                                      !
  !        --> value = ?_mat(idx:idx+NPNTS-1)                          !
  !                                                                    !
  !--------------------------------------------------------------------!
  
  logical,                                         private :: bi_isInit = .false.

  integer,                                         private :: ibMat
  integer,                                         private :: ibOns
  integer,                                         private :: nBonds
  integer,                                         private :: nOns
  integer,                                         private :: memSize
  integer,                                         private :: memSizeOns
  integer,                                         private :: memSizePot

  integer,          dimension(:,:,:), allocatable, private :: bi
  integer,                                         private :: nbi
  integer,          dimension(:,:,:), allocatable, private :: biIdx
  integer,          dimension(:,:,:), allocatable, private :: onsIdx

  integer,                                         private :: bi_ntypes
  integer,                                         private :: bi_l_min
  integer,                                         private :: bi_l_max
  double precision,                                private :: bi_r_min
  double precision,                                private :: bi_r_max
  double precision,                                private :: bi_r_min_pot
  double precision,                                private :: bi_r_max_pot

  logical,                                         private :: bi_hasHmat
  logical,                                         private :: bi_hasSmat
  logical,                                         private :: bi_hasOmat
  logical,                                         private :: bi_hasV

  double precision, dimension(:,:),   allocatable, private :: onsite_level

  double precision, dimension(:), allocatable, target, private :: R_mat
  double precision, dimension(:), allocatable, target, private :: H_mat, d2H_mat
  double precision, dimension(:), allocatable, target, private :: S_mat, d2S_mat

  double precision, dimension(:), allocatable, target, private :: R_ons
  double precision, dimension(:), allocatable, target, private :: H_ons, d2H_ons 

  double precision, dimension(:), allocatable,         private :: R_pot
  double precision, dimension(:), allocatable, target, private :: V_pot, d2V_pot

  interface bi_tabulate
     module procedure bi_tabulate_mat, bi_tabulate_pot
  end interface

contains

  !--------------------------------------------------------------------!
  !               module initialization and finalization               !
  !--------------------------------------------------------------------!

  subroutine bi_init(ntypes, nlmax, nl, l_of_il, r0, r1, hamilton, &
                     notb, onsite, r0pot, r1pot, pairpot)
    
    implicit none

    !------------------------------------------------------------------!
    ! ntypes            : number of atomic species                     !
    ! l_max             : highest angular momentum                     !
    ! nl(itype)         : num. angular momentum channels of type itype !
    ! l_of_il(il,itype) : il-th angular momentum of type itype         !
    ! r0, r1            : lower and upper cut-off for r                !
    ! notb              : .true., if overlap is parametrized           !
    ! onsite            : .true., if onsite elements are parametrized  !
    ! r0pot, r1pot      : r range for the pair potentials              !
    ! pairpot           : .true., if pair pot. parameters are available!
    !------------------------------------------------------------------!

    integer,                             intent(in) :: ntypes
    integer,                             intent(in) :: nlmax
    integer, dimension(ntypes),          intent(in) :: nl
    integer, dimension(nlmax, ntypes),   intent(in) :: l_of_il
    double precision,                    intent(in) :: r0, r1
    logical,                             intent(in) :: hamilton
    logical,                             intent(in) :: notb
    logical,                             intent(in) :: onsite
    double precision,                    intent(in) :: r0pot, r1pot
    logical,                             intent(in) :: pairpot

    integer          :: itype1, itype2, il1, il2, ipt
    integer          :: l1, l2, l_min, l_max
    double precision :: r, dr

    ! paramterizations provided:
    bi_hasHmat = hamilton
    bi_hasSmat = hamilton
    bi_hasOmat = onsite
    bi_hasV    = pairpot

    ibMat = 1
    ibOns = 1

    !-----------------------------------------!
    ! loop over all angular momentum channel  !
    ! combinations of all atom types          !
    ! --> count all possible bonds = nbonds   !
    !     determine min. angular mom. = l_min !
    !-----------------------------------------!

    nBonds = 0
    nOns   = 0

    l_min  = l_of_il(1, 1)
    l_max  = l_min

    do itype1 = 1, ntypes
    do il1 = 1, nl(itype1)
       l1     = l_of_il(il1, itype1)
       l_min  = min(l_min, l1)
       l_max  = max(l_max, l1)
       ! same type, same angular momentum:
       nBonds = nBonds + l1 + 1
       nOns   = nOns   + ntypes*(l1 + 1)
       ! same type, different angular momentum:
       do il2 = il1 + 1, nl(itype1)
          l2     = l_of_il(il2, itype1)
          nBonds = nBonds + min(l1,l2) + 1
          nOns   = nOns   + ntypes*(min(l1,l2) + 1)
       end do
       ! different type:
       do itype2 = itype1 + 1, ntypes
          ! bonds:
          do il2 = 1, nl(itype2)
             l2     = l_of_il(il2, itype2)
             nBonds = nBonds + min(l1,l2) + 1
          end do
       end do
    end do
    end do

    
    !----------------------------------------!
    ! number of different bond integrals     !
    ! and bond IDs:                          !
    !                                        !
    ! defined afterwards: nbi, bi(l1, l2, m) !
    !----------------------------------------!

    call generate_bondint_IDs(l_min, l_max)


    !-----------------!
    ! allocate memory !
    !-----------------!

    memSize    = NPNTS*nBonds
    memSizeOns = NPNTS*nOns
    memSizePot = (ntypes*(ntypes + 1))/2*NPNTS

    allocate( biIdx(nbi, ntypes, ntypes), onsIdx(nbi,ntypes,ntypes), &
              R_mat(NPNTS), onsite_level(l_min:l_max, ntypes) )

    if (hamilton) then
       allocate(H_mat(memSize), d2H_mat(memSize))
       H_mat(:)   = 0.0d0
       d2H_mat(:) = 0.0d0
    end if
    if (notb) then
       allocate(S_mat(memSize), d2S_mat(memSize))
       S_mat(:)   = 0.0d0
       d2S_mat(:) = 0.0d0
    end if
    if (onsite) then
       allocate(R_ons(NPNTS), H_ons(memSizeOns), d2H_ons(memSizeOns))
       H_ons(:)   = 0.0d0
       d2H_ons(:) = 0.0d0
    end if
    if (pairpot) then
       allocate(R_pot(NPNTS), V_pot(memSizePot), d2V_pot(memSizePot))
       V_pot(:)   = 0.0d0
       d2V_pot(:) = 0.0d0
    end if

    biIdx(1:nbi, 1:ntypes, 1:ntypes)  = 0
    onsIdx(1:nbi, 1:ntypes, 1:ntypes) = 0
    
    
    !------------------!
    ! list of r values !
    !------------------!

    ! matrix elements:
    dr = (r1 - r0)/dble(NPNTS - 1)
    r  = r0
    do ipt = 1, NPNTS
       R_mat(ipt) = r
       r = r + dr
    end do

    ! on-site contributions:
    if (onsite) R_ons(:) = R_mat(:)

    ! pair potentials:
    if (pairpot) then
       dr = (r1pot - r0pot)/dble(NPNTS - 1)
       r  = r0
       do ipt = 1, NPNTS
          R_pot(ipt) = r
          r = r + dr
       end do
    end if


    !--------------------!
    ! things to remember !
    !--------------------!

    bi_l_min     = l_min
    bi_l_max     = l_max
    bi_r_min     = r0
    bi_r_max     = r1
    bi_r_min_pot = r0pot
    bi_r_max_pot = r1pot
    bi_ntypes    = ntypes

    bi_isInit = .true.

  end subroutine bi_init

  !--------------------------------------------------------------------!

  subroutine bi_final()

    implicit none

    if (allocated(bi))           deallocate(bi)
    if (allocated(biIdx))        deallocate(biIdx)
    if (allocated(onsIdx))       deallocate(onsIdx)
    if (allocated(R_mat))        deallocate(R_mat)
    if (allocated(H_mat))        deallocate(H_mat, d2H_mat)
    if (allocated(S_mat))        deallocate(S_mat, d2S_mat)
    if (allocated(R_ons))        deallocate(R_ons)
    if (allocated(H_ons))        deallocate(H_ons, d2H_ons)
    if (allocated(R_pot))        deallocate(R_pot)
    if (allocated(V_pot))        deallocate(V_pot, d2V_pot)
    if (allocated(onsite_level)) deallocate(onsite_level)

    bi_isInit = .false.

  end subroutine bi_final

  !--------------------------------------------------------------------!
  !                      unique bond integral ID                       !
  !--------------------------------------------------------------------!

  subroutine generate_bondint_IDs(l_min, l_max)

    implicit none

    integer, intent(in) :: l_min, l_max

    integer :: l1, l2, m

    allocate(bi(l_min:l_max, l_min:l_max, 0:l_max))

    nbi = 0
    do l1 = l_min, l_max
    do l2 = l_min, l_max
    do m  = 0, min(l1, l2)
       nbi         = nbi + 1
       bi(l1,l2,m) = nbi
    end do
    end do
    end do

  end subroutine generate_bondint_IDs

  !--------------------------------------------------------------------!
  !       unique index for each combination of itype1 >= itype2        !
  !--------------------------------------------------------------------!

  function bi_pot_idx(itype1, itype2) result(idx)

    implicit none

    integer, intent(in) :: itype1, itype2
    integer             :: idx

    integer             :: it1, it2

    ! assert that type1 >= type2
    if (itype1 >= itype2) then
       it1 = itype1
       it2 = itype2
    else
       it1 = itype2
       it2 = itype1
    end if
    
    idx = (it1*(it1 - 1))/2 + it2 - 1
    idx = idx*NPNTS + 1

  end function bi_pot_idx

  !--------------------------------------------------------------------!
  !                   look up matrix element values                    !
  !--------------------------------------------------------------------!

  subroutine bi_bondint(x, itype1, itype2, l1, l2, m, H, S, O, dH, dS, DO)

    implicit none

    double precision,               intent(in)  :: x
    integer,                        intent(in)  :: itype1, itype2
    integer,                        intent(in)  :: l1, l2, m
    double precision, optional,     intent(out) :: H, dH, S, dS, O, dO

    integer          :: idx1H, idx2H, idx1S, idx2S, idx1O, idx2O
    double precision :: facH, facS, facO

    if (x < bi_r_min) then
       write(0,*) "Error: out of parametrized region for r = ", x
       stop
    end if
    
    ! retrieve memory indices:
    if (present(H) .or. present(dH)) then
       facH = 1.0d0
       call bi_memory_index('hamilton', itype1, itype2, l1, l2, m, &
                            idx1H, idx2H, facH)       
    end if
    if (present(S) .or. present(dS)) then
       facS = 1.0d0
       call bi_memory_index('overlap', itype1, itype2, l1, l2, m, &
                            idx1S, idx2S, facS)       
    end if
    if (present(O) .or. present(dO)) then
       facO = 1.0d0
       call bi_memory_index('onsite', itype1, itype2, l1, l2, m, &
                            idx1O, idx2O, facO)
    end if

    ! interpolate parameter curves:
    if (x > bi_r_max) then
       if (present(H))  H  = 0.0d0
       if (present(dH)) dH = 0.0d0
       if (present(S))  S  = 0.0d0
       if (present(dS)) dS = 0.0d0
       if (present(O))  O  = 0.0d0
       if (present(dO)) dO = 0.0d0
    else
       if (present(H)) then
          call cs(R_mat(1:NPNTS), H_mat(idx1H:idx2H), &
                  d2H_mat(idx1H:idx2H), NPNTS, x, H)
          H = facH*H
       end if
       if (present(dH)) then
          call cs_d1(R_mat(1:NPNTS), H_mat(idx1H:idx2H), &
                     d2H_mat(idx1H:idx2H), NPNTS, x, dH)
          dH = facH*dH
       end if
       if (present(S)) then
          call cs(R_mat(1:NPNTS), S_mat(idx1S:idx2S), &
                  d2S_mat(idx1S:idx2S), NPNTS, x, S)
          S = facS*S
       end if
       if (present(dS)) then
          call cs_d1(R_mat(1:NPNTS), S_mat(idx1S:idx2S), &
                     d2S_mat(idx1S:idx2S), NPNTS, x, dS)
          dS = facS*dS
       end if
       if (present(O)) then
          call cs(R_ons(1:NPNTS), H_ons(idx1O:idx2O), &
                  d2H_ons(idx1O:idx2O), NPNTS, x, O)
          O = facO*O
       end if
       if (present(dO)) then
          call cs_d1(R_ons(1:NPNTS), H_ons(idx1O:idx2O), &
                  d2H_ons(idx1O:idx2O), NPNTS, x, dO)
          dO = facO*dO
       end if
    end if

  end subroutine bi_bondint

  !--------------------------------------------------------------------!
  !                    look up pair potential value                    !
  !--------------------------------------------------------------------!

  subroutine bi_pairpot(x, itype1, itype2, V, dV)

    implicit none

    double precision,               intent(in)  :: x
    integer,                        intent(in)  :: itype1, itype2
    double precision, optional,     intent(out) :: V, dV

    integer          :: idx1, idx2

    idx1 = bi_pot_idx(itype1, itype2)
    idx2 = idx1 + NPNTS - 1

    if (x < bi_r_min_pot) then
       write(0,*) "Error: out of parametrized region for r = ", x
       stop
    end if
    
    if (x > bi_r_max_pot) then
       if (present(V))  V  = 0.0d0
       if (present(dV)) dV = 0.0d0
    else
       if (present(V)) then
          call cs(R_pot(1:NPNTS), V_pot(idx1:idx2), &
                  d2V_pot(idx1:idx2), NPNTS, x, V)
       end if
       if (present(dV)) then
          call cs_d1(R_pot(1:NPNTS), V_pot(idx1:idx2), &
                  d2V_pot(idx1:idx2), NPNTS, x, dV)
       end if
    end if

  end subroutine bi_pairpot

  !--------------------------------------------------------------------!
  !                       look up on-site level                        !
  !--------------------------------------------------------------------!

  function bi_onsite(l, itype) result(level)

    implicit none

    integer, intent(in) :: l
    integer, intent(in) :: itype
    double precision    :: level

    if ((l >= bi_l_min) .and. (l <= bi_l_max) .and. (itype <= bi_ntypes)) then
       level = onsite_level(l, itype)
    else
       write(0,*) "Warning: invalid request for onsite level in `bi_onsite'."
       level = 0.0d0
    end if
    
  end function bi_onsite

  !--------------------------------------------------------------------!
  !                         property requests                          !
  !                                                                    !
  ! Values of some variables can be accessed read-only via these       !
  ! functions.  See comments for variables (above) for explanations.   !
  !--------------------------------------------------------------------!

  function bi_get_isInit() result(val)
    implicit none
    logical :: val
    val = bi_isInit
  end function bi_get_isInit

  !--------------------------------------------------------------------!

  function bi_get_hasHmat() result(val)
    implicit none
    logical :: val
    val = bi_hasHmat
  end function bi_get_hasHmat
  
  !--------------------------------------------------------------------!

  function bi_get_hasSmat() result(val)
    implicit none
    logical :: val
    val = bi_hasSmat
  end function bi_get_hasSmat
  
  !--------------------------------------------------------------------!

  function bi_get_hasOmat() result(val)
    implicit none
    logical :: val
    val = bi_hasOmat
  end function bi_get_hasOmat
  
  !--------------------------------------------------------------------!

  function bi_get_hasV() result(val)
    implicit none
    logical :: val
    val = bi_hasV
  end function bi_get_hasV
  
  !--------------------------------------------------------------------!

  function bi_get_r_min() result(val)
    implicit none
    double precision :: val
    val = bi_r_min
  end function bi_get_r_min

  !--------------------------------------------------------------------!

  function bi_get_r_max() result(val)
    implicit none
    double precision :: val
    val = bi_r_max
  end function bi_get_r_max

  !--------------------------------------------------------------------!

  function bi_get_r_min_pot() result(val)
    implicit none
    double precision :: val
    val = bi_r_min_pot
  end function bi_get_r_min_pot

  !--------------------------------------------------------------------!

  function bi_get_r_max_pot() result(val)
    implicit none
    double precision :: val
    val = bi_r_max_pot
  end function bi_get_r_max_pot

  !--------------------------------------------------------------------!

  function bi_get_l_min() result(val)
    implicit none
    integer :: val
    val = bi_l_min
  end function bi_get_l_min

  !--------------------------------------------------------------------!

  function bi_get_l_max() result(val)
    implicit none
    integer :: val
    val = bi_l_max
  end function bi_get_l_max

  !--------------------------------------------------------------------!

  function bi_get_nbi() result(val)
    implicit none
    integer :: val
    val = nbi
  end function bi_get_nbi

  !--------------------------------------------------------------------!

  function bi_get_bi(l1, l2, m) result(val)
    implicit none
    integer, intent(in) :: l1, l2, m
    integer             :: val
    val = bi(l1,l2,m)
  end function bi_get_bi

  !--------------------------------------------------------------------!
  !                       fill parameter tables                        !
  !--------------------------------------------------------------------!

  subroutine bi_tabulate_mat(scheme, itype1, itype2, matrix, l1, l2, m, &
                             r0, r1, params, xscale, yscale)

    implicit none

    character(len=*),               intent(in) :: scheme
    integer,                        intent(in) :: itype1, itype2
    character(len=*),               intent(in) :: matrix
    integer,                        intent(in) :: l1, l2, m
    double precision,               intent(in) :: r0, r1
    double precision, dimension(:), intent(in) :: params
    double precision, optional,     intent(in) :: xscale, yscale

    double precision, dimension(:), pointer    :: A, d2A, R

    double precision :: xfac, yfac
    integer          :: idx1, idx2
  
    ! scaling factors (e.g. for unit conversion):
    if (present(xscale)) then
       xfac = xscale
    else
       xfac = 1.0d0
    end if
    if (present(yscale)) then
       yfac = yscale
    else
       yfac = 1.0d0
    end if

    ! check boundaries of r:
    if ((xfac*r0 < bi_r_min) .or. (xfac*r1 > bi_r_max)) then
       write(0,*) itype1, itype2, matrix, l1, l2, m, xfac*r0, xfac*r1
       write(0,*) "Error: out of distance interval in `bi_tabulate()'"
       stop
    end if

    ! retrieve memory positions:
    call bi_memory_index(matrix, itype1, itype2, l1, l2, m, idx1, idx2, yfac)
    
    ! select matrix elements' memory region:
    sel_matrix : select case(trim(io_upper(matrix)))
    case('H', 'HAMILTON', 'HOPPING')
       A   => H_mat(idx1:idx2)
       d2A => d2H_mat(idx1:idx2)
       R   => R_mat(:)
    case('S', 'OVERLAP', 'SCREENING')
       A   => S_mat(idx1:idx2)
       d2A => d2S_mat(idx1:idx2)
       R   => R_mat(:)
    case('O', 'ONSITE')
       A   => H_ons(idx1:idx2)
       d2A => d2H_ons(idx1:idx2)
       R   => R_ons(:)
    case default      
       write(0,*) "Error: unknown matrix name in `bi_tabulate()': ", &
            trim(adjustl(matrix))
       return
    end select sel_matrix

    ! evaluate and store parameter functions:
    sel_scheme : select case(trim(io_lower(scheme)))
    case('chebychev')
       call tab_chebychev(params, r0, r1, xfac, yfac, R, A, d2A)
    case('csplines')
       call tab_csplines(params, r1, xfac, yfac, R, A, d2A)
    case('polyexp')
       call tab_polyexp(params, r1, xfac, yfac, R, A, d2A)
    case default
       write(0,*) "Error: unknown parameter scheme in `bi_tabulate()': ", &
            trim(adjustl(scheme))
       return
    end select sel_scheme
  
  end subroutine bi_tabulate_mat

  !--------------------------------------------------------------------!

  subroutine bi_tabulate_pot(scheme, itype1, itype2, pottype, &
                             r0, r1, params, xscale, yscale)

    implicit none

    character(len=*),               intent(in) :: scheme
    integer,                        intent(in) :: itype1, itype2
    character(len=*),               intent(in) :: pottype
    double precision,               intent(in) :: r0, r1
    double precision, dimension(:), intent(in) :: params
    double precision, optional,     intent(in) :: xscale, yscale

    double precision, dimension(:), pointer    :: A, d2A

    double precision :: xfac, yfac
    integer          :: idx
  
    ! scaling factors (e.g. for unit conversion):
    if (present(xscale)) then
       xfac = xscale
    else
       xfac = 1.0d0
    end if
    if (present(yscale)) then
       yfac = yscale
    else
       yfac = 1.0d0
    end if

    ! check boundaries of r:
    if ((xfac*r0 < bi_r_min_pot) .or. (xfac*r1 > bi_r_max_pot)) then
       write(0,*) xfac*r0, bi_r_min_pot, xfac*r1, bi_r_max_pot
       write(0,*) itype1, itype2, pottype
       write(0,*) "Error: out of distance interval in `bi_tabulate()'"
       stop
    end if

    ! index for the parameters entry:
    idx = bi_pot_idx(itype1, itype2)

    ! select potential type:
    select  case(trim(io_lower(pottype)))
    case('pairpot')
       A   => V_pot(idx:idx+NPNTS-1)
       d2A => d2V_pot(idx:idx+NPNTS-1)
    case default
       write(0,*) "Error: unknown potential type name in `bi_tabulate()': ", &
            trim(adjustl(pottype))
       return       
    end select

    ! select evaluation routine:
    sel_scheme : select case(trim(io_lower(scheme)))
    case('chebychev')
       call tab_chebychev(params, r0, r1, xfac, yfac, R_pot, A, d2A)
    case('csplines')
       call tab_csplines(params, r1, xfac, yfac, R_pot, A, d2A)
    case('polyexp')
       call tab_polyexp(params, r1, xfac, yfac, R_pot, A, d2A)
    case default
       write(0,*) "Error: unknown parameter scheme in `bi_tabulate()': ", &
            trim(adjustl(scheme))
       return
    end select sel_scheme
    
  end subroutine bi_tabulate_pot

  !--------------------------------------------------------------------!
  !   standardize parameter order for itype1, itype2, l1, l2, m        !
  !                                                                    !
  ! Every bond parametrization is stored just once, but the symmetry   !
  ! with respect to the exchange of the atom types must be obeyed.     !
  !                                                                    !
  ! The following rules are used for (itype1, itype2, l1, l2, m)       !
  ! - itype1 <= itype2                                                 !
  ! - if (itype1 == itype2) ---> l1 <= l2                              !
  !                                                                    !
  ! For distance dependent on-site contributions there is NO SYMMETRY  !
  ! with respect to exchange of the atom type, because both orbitals   !
  ! are on the same atom on either one of the two atom types.  On the  !
  ! other hand only one l1/l2 combination has to be stored.            !
  !--------------------------------------------------------------------!

  subroutine bi_memory_index(dataType, itype1, itype2, l1,  l2, m, &
                             idx1, idx2, fac )

    implicit none

    !------------------------------------------------------------------!
    ! dataType       : currently either "hamilton", "overlap"          !
    !                  or "onsite" ("h", "s", "o")                     !
    ! itype1, itype2 : input atom types                                !
    ! l1, l2         : input angular momenta                           !
    ! m              : common magnetic quantum number                  !
    ! idx1, idx2     : memory range indices for the selected data type !
    ! fac            : scaling factor, that might change the sign      !
    !------------------------------------------------------------------!

    character(len=*), intent(in)    :: dataType
    integer,          intent(in)    :: itype1, itype2, l1, l2, m
    integer,          intent(out)   :: idx1, idx2
    double precision, intent(inout) :: fac

    integer :: it1, it2, il1, il2
    integer :: ibond

    select case(io_lower(trim(dataType)))

       !----------------------- bond integrals ------------------------!

    case('s','overlap','screening','h','hamilton','hopping')
       if (itype1 == itype2) then
          it1 = itype1
          it2 = itype2
          if (l1 <= l2) then
             il1 = l1
             il2 = l2
          else
             il1 = l2
             il2 = l1
             if (mod(il2-il1,2) == 1) fac = -fac
          end if
       else if (itype1 <= itype2) then
          it1 = itype1
          it2 = itype2
          il1 = l1
          il2 = l2
       else
          it1 = itype2
          it2 = itype1
          il1 = l2
          il2 = l1
          if (mod(abs(il2-il1),2) == 1) fac = -fac
       end if
       ibond = bi(il1, il2, m)

       if (biIdx(ibond, it1, it2) == 0) then
          if (ibMat > (memSize - NPNTS + 1)) then
             write(*,*) itype1, itype2
             write(*,*) trim(dataType), l1, l2, m
             write(0,*) "Error: allocated memory exceeded in `bondint' module."
             stop
          else
             biIdx(ibond, it1, it2) = ibMat
             ibMat = ibMat + NPNTS
          end if
       end if

       idx1 = biIdx(ibond, it1, it2)
       idx2 = biIdx(ibond, it1, it2) + NPNTS - 1

       !--------------------- on-site parameters ----------------------!

    case('o','onsite')
       it1 = itype1
       it2 = itype2
       if (l1 <= l2) then
          il1 = l1
          il2 = l2
       else
          il1 = l2
          il2 = l1
       end if
       ibond = bi(il1, il2, m)

       if (onsIdx(ibond, it1, it2) == 0) then
          if (ibOns > (memSizeOns - NPNTS + 1)) then
             write(*,*) itype1, itype2
             write(*,*) trim(dataType), l1, l2, m
             write(0,*) "Error: allocated memory exceeded in `bondint' module."
             stop
          else
             onsIdx(ibond, it1, it2) = ibOns
             ibOns = ibOns + NPNTS
          end if
       end if

       idx1 = onsIdx(ibond, it1, it2)
       idx2 = onsIdx(ibond, it1, it2) + NPNTS - 1

       !---------------------------------------------------------------!

    case default
       write(0,*) "Error: unknown data type in `bi_standard_index()': ", &
            trim(adjustl(dataType))
       stop
    end select
       
 
  end subroutine bi_memory_index

  !--------------------------------------------------------------------!
  !                        store on-site level                         !
  !--------------------------------------------------------------------!

  subroutine bi_set_onsitelevel(l, itype, level)

    implicit none

    integer,          intent(in) :: l
    integer,          intent(in) :: itype
    double precision, intent(in) :: level

    if ((l >= bi_l_min) .and. (l <= bi_l_max) .and. (itype <= bi_ntypes)) then
       onsite_level(l, itype) = level
    end if

  end subroutine bi_set_onsitelevel

  !--------------------------------------------------------------------!
  !                   parameter function evaluation                    !
  !--------------------------------------------------------------------!

  subroutine tab_chebychev(params, rmin, rcut, xfac, yfac, X, A, d2A)
    
    implicit none

    double precision, dimension(:),     intent(in)  :: params
    double precision,                   intent(in)  :: rmin, rcut
    double precision,                   intent(in)  :: xfac, yfac
    double precision, dimension(NPNTS), intent(in)  :: X
    double precision, dimension(NPNTS), intent(out) :: A, d2A

    integer          :: nparams, ncut
    integer          :: ipt
    double precision :: r, dA_1, dA_n

    nparams = size(params)

    ncut = NPNTS
    do ipt = 1, NPNTS
       r = X(ipt)/xfac
       if (r < rcut) then
          A(ipt)  = yfac*chebypol(r, rmin, rcut, nparams, params)
       else
          ncut = ipt - 1
          exit
       end if
    end do

    dA_1 = yfac*dchebypol(X(1)/xfac, rmin, rcut, nparams, params)
    dA_n = yfac*dchebypol(X(ncut)/xfac, rmin, rcut, nparams, params)

    call deriv2(X(1:ncut), A(1:ncut), ncut, dA_1, dA_n, d2A(1:ncut))

    if (ncut < NPNTS) then
       A(ncut+1:NPNTS)   = 0.0d0
       d2A(ncut+1:NPNTS) = 0.0d0
    end if

  end subroutine tab_chebychev

  !--------------------------------------------------------------------!

  subroutine tab_csplines(params, rcut, xfac, yfac, R, A, d2A)
    
    implicit none

    double precision, dimension(:),     intent(in)  :: params
    double precision,                   intent(in)  :: rcut
    double precision,                   intent(in)  :: xfac, yfac
    double precision, dimension(NPNTS), intent(in)  :: R
    double precision, dimension(NPNTS), intent(out) :: A, d2A

    integer          :: i, j, n, ipt, ncut
    double precision :: y1_1, y1_n, x0
    double precision, dimension(:), allocatable :: x, y, y2

    n = (size(params) - 2)/2
    y1_1 = params(1)
    y1_n = params(2)

    allocate(x(n), y(n), y2(n))

    j = 0
    do i = 3, size(params), 2
       j = j + 1
       x(j) = params(i)
       y(j) = params(i+1)
    end do
    call deriv2(x, y, n, y1_1, y1_n, y2)

    ncut = NPNTS
    do ipt = 1, NPNTS
       x0 = R(ipt)/xfac
       if (x0 < rcut) then
          call cs(x, y, y2, n, x0, A(ipt))
          A(ipt) = yfac*A(ipt)
       else
          ncut = ipt - 1
          exit
       end if
    end do
    call cs_d1(x, y, y2, n, R(1)/xfac, y1_1)
    call cs_d1(x, y, y2, n, R(ncut)/xfac, y1_n)
    y1_1 = yfac*y1_1
    y1_n = yfac*y1_n

    call deriv2(R(1:ncut), A(1:ncut), ncut, y1_1, y1_n, d2A(1:ncut))

    if (ncut < NPNTS) then
       A(ncut+1:NPNTS)   = 0.0d0
       d2A(ncut+1:NPNTS) = 0.0d0
    end if

    deallocate(x, y, y2)

  end subroutine tab_csplines

  !--------------------------------------------------------------------!

  subroutine tab_polyexp(params, rcut, xfac, yfac, X, A, d2A)

    implicit none

    double precision, dimension(:),     intent(in)  :: params
    double precision,                   intent(in)  :: rcut
    double precision,                   intent(in)  :: xfac, yfac
    double precision, dimension(NPNTS), intent(in)  :: X
    double precision, dimension(NPNTS), intent(out) :: A, d2A

    integer          :: order
    double precision :: r
    double precision :: y, dy, d2y, d3y
    double precision :: yc, dyc, d2yc, d3yc
    double precision :: fc, dfc, d2fc, d3fc
    double precision :: y1_1, y1_n
    integer          :: ipt, ncut

    order = size(params) - 2

    call polyexp(rcut, params(1), order, params(2:order+2), yc, dyc, d2yc, d3yc)

    ncut = NPNTS
    do ipt = 1, NPNTS
       r = X(ipt)/xfac
       if (r < rcut) then
          call polyexp(r, params(1), order, params(2:order+2), y, dy, d2y, d3y)
          call polyexp_cut(r, rcut, yc, dyc, d2yc, fc, dfc, d2fc, d3fc)
          A(ipt)   = yfac*(y   - fc)
       else
          ncut = ipt - 1
          exit
       end if
    end do

    call polyexp(X(1)/xfac, params(1), order, params(2:order+2), y, dy, d2y, d3y)
    call polyexp_cut(X(1)/xfac, rcut, yc, dyc, d2yc, fc, dfc, d2fc, d3fc)
    y1_1 = yfac*(dy - dfc)
    call polyexp(X(ncut)/xfac, params(1), order, params(2:order+2), y, dy, d2y, d3y)
    call polyexp_cut(X(ncut)/xfac, rcut, yc, dyc, d2yc, fc, dfc, d2fc, d3fc)
    y1_n = yfac*(dy - dfc)

    call deriv2(X(1:ncut), A(1:ncut), ncut, y1_1, y1_n, d2A(1:ncut))

    if (ncut < NPNTS) then
       A(ncut+1:NPNTS)     = 0.0d0
       d2A(ncut+1:NPNTS)   = 0.0d0
    end if

  end subroutine tab_polyexp

  !--------------------------------------------------------------------!

!!$  subroutine tab_polyexp(params, rcut, fac, X, A, d2A)
!!$
!!$    implicit none
!!$
!!$    double precision, dimension(:),     intent(in)  :: params
!!$    double precision,                   intent(in)  :: rcut
!!$    double precision,                   intent(in)  :: fac
!!$    double precision, dimension(NPNTS), intent(in)  :: X
!!$    double precision, dimension(NPNTS), intent(out) :: A, d2A
!!$
!!$    integer          :: order
!!$    double precision :: r
!!$    double precision :: y, dy, d2y, d3y
!!$    double precision :: yc, dyc, d2yc, d3yc
!!$    double precision :: fc, dfc, d2fc, d3fc
!!$    integer          :: ipt
!!$
!!$    order = size(params) - 2
!!$
!!$    call polyexp(rcut, params(1), order, params(2:order+2), yc, dyc, d2yc, d3yc)
!!$
!!$    do ipt = 1, NPNTS
!!$       r = X(ipt)
!!$       if (r < rcut) then
!!$          call polyexp(r, params(1), order, params(2:order+2), y, dy, d2y, d3y)
!!$          call polyexp_cut(r, rcut, yc, dyc, d2yc, fc, dfc, d2fc, d3fc)
!!$          A(ipt)   = fac*(y   - fc)
!!$          d2A(ipt) = fac*(d2y - d2fc)
!!$       else
!!$          A(ipt)   = 0.0d0
!!$          d2A(ipt) = 0.0d0
!!$       end if
!!$    end do
!!$
!!$  end subroutine tab_polyexp

  !--------------------------------------------------------------------!
  !                              polynomial                            !
  !--------------------------------------------------------------------!

  subroutine polynomial(x, order, coeff, y, dy, d2y, d3y)

    implicit none

    double precision,                     intent(in)  :: x
    integer,                              intent(in)  :: order
    double precision, dimension(0:order), intent(in)  :: coeff
    double precision,                     intent(out) :: y, dy, d2y, d3y

    integer          :: i

    ! polynomial:
    y = 0.0d0
    do i = order, 1, -1
       y = x*(coeff(i) + y)
    end do
    y = y + coeff(0)

    ! 1st derivative:
    dy = 0.0d0
    do i = order, 2, -1
       dy = x*(dble(i)*coeff(i) + dy)
    end do
    if (order >= 1) dy = dy + coeff(1)

    ! 2nd derivative:
    d2y = 0.0d0
    do i = order, 3, -1
       d2y = x*(dble(i*(i-1))*coeff(i) + d2y)
    end do
    if (order >= 2) d2y = d2y + 2.0d0*coeff(2)

    ! 3rd derivative:
    d3y = 0.0d0
    do i = order, 4, -1
       d3y = x*(dble(i*(i-1)*(i-2))*coeff(i) + d3y)
    end do
    if (order >= 3) d3y = d3y + 6.0d0*coeff(3)

  end subroutine polynomial


  !--------------------------------------------------------------------!
  !      Evaluate Chebyshev polynomial with coefficients c(1:n).       !
  !                 [see Numerical Recipes, Chapter 5.8]               !
  !--------------------------------------------------------------------!

  function chebypol(r, r0, r1, n, para) result(val)

    implicit none

    !------------------------------------------------------------------!
    ! r            : distance; para(1) <= r <= para(2)                 !
    ! r0, r1       : lower and upper boundaries                        !
    ! n            : number of coefficients                            !
    ! para(1...n)  : coefficients                                      !
    ! val          : value of the Chebychev polynomial                 !
    !------------------------------------------------------------------!

    double precision,               intent(in) :: r, r0, r1
    integer,                        intent(in) :: n
    double precision, dimension(n), intent(in) :: para
    double precision                           :: val

    integer                      :: i
    double precision             :: d1, d2, sv, y, y2

    val = 0.d0

    ! Change of variable (rescale to interval [-1,1]):
    y = (2.d0*r - r0 - r1)/(r1 - r0)
    y2 = 2.d0*y

    ! Evaluate Chebychev polynomial (Clenshaw's recurrence):
    d1 = 0.d0
    d2 = 0.d0
    do i = n, 2, -1
       sv = d1
       d1 = y2*d1 - d2 + para(i)
       d2 = sv
    end do
    val = y*d1 - d2 + 0.5d0*para(1)

  end function chebypol

  !--------------------------------------------------------------------!
  !         Evaluate the derivative of a Chebyshev polynomial          !
  !                     with coefficients c(1:n).                      !
  !             [see Numerical Recipes, Chapter 5.5 + 5.8]             !
  !--------------------------------------------------------------------!

  function dchebypol(r, r0, r1, n, para) result(dval)
    
    implicit none

    !------------------------------------------------------------------!
    ! r            : distance; para(1) <= r <= para(2)                 !
    ! r0, r1       : lower and upper boundaries                        !
    ! n            : number of parameters                              !
    ! para(1...n)  : coefficients                                      !
    ! dval         : derivative of the Chebychev polynomial            !
    !------------------------------------------------------------------!

    double precision,               intent(in) :: r, r0, r1
    integer,                        intent(in) :: n
    double precision, dimension(n), intent(in) :: para
    double precision                           :: dval

    integer                       :: i
    double precision              :: d1, d2, denom, sv, y2

    dval = 0.d0

    ! Change of variable (rescale to interval [-1,1]):
    denom = 2.D0/(r1 - r0)
    y2 = (2.d0*r - r0 - r1)*denom
    
    ! Evaluate derivative of Chebychev polynomial (Clenshaw's recurrence):
    d1 = 0.d0
    d2 = 0.d0
    do i = n, 3, -1
       sv = d1
       d1 = y2*d1 - d2 + dble(i-1)*para(i)
       d2 = sv
    enddo
    dval = y2*d1 - d2 + para(2)

    ! Renormalize (due to rescaling):
    dval = denom*dval

  end function dchebypol

  !--------------------------------------------------------------------!
  !                 polyexp = polynomial x exponential                 !
  !--------------------------------------------------------------------!

  subroutine polyexp(x, a, order, coeff, y, dy, d2y, d3y)

    implicit none

    double precision,                     intent(in)  :: x
    double precision,                     intent(in)  :: a
    integer,                              intent(in)  :: order
    double precision, dimension(0:order), intent(in)  :: coeff
    double precision,                     intent(out) :: y, dy, d2y, d3y
    
    double precision :: p, dp, d2p, d3p
    double precision :: expval

    call polynomial(x, order, coeff, p, dp, d2p, d3p)
    
    expval = exp(-a*x)

    y   = p*expval
    dy  = (dp - a*p)*expval
    d2y = (d2p - a*(2.0d0*dp - a*p))*expval
    d3y = (d3p - a*(3.0d0*d2p - a*(3.0d0*dp - a*p)))*expval

  end subroutine polyexp

  !--------------------------------------------------------------------!

  subroutine polyexp_cut(x, xcut, fc, dfc, d2fc, y, dy, d2y, d3y)

    implicit none

    double precision,                     intent(in)  :: x, xcut
    double precision,                     intent(in)  :: fc, dfc, d2fc
    double precision,                     intent(out) :: y, dy, d2y, d3y

    y   = fc + (x - xcut)*dfc + 0.5d0*(x - xcut)*(x - xcut)*d2fc
    dy  = dfc + (x - xcut)*d2fc
    d2y = d2fc
    d3y = 0

  end subroutine polyexp_cut


end module bondint
