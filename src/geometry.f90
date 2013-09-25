module geometry

  !--------------------------------------------------------------------!
  ! Structural and related information provided by one or more input   !
  ! files are handled in this module.                                  !
  !   In principle the idea is to provide a unified interface for      !
  ! different input file formats.  Nonetheless, currently only a single!
  ! format is supported; namely the format of Bernd Meyer's Mixed Basis!
  ! Pseudo-Potential (MBPP) code.                                      !
  !--------------------------------------------------------------------!
  ! 2010-10-29 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use constants, only: PI2, eV2Ha, Ang2Bohr

  use io,     only: io_adjustl,              &
                    io_lower,                &
                    io_upper

  implicit none

  public  :: geo_init,              &
             geo_final,             &
             geo_print_status,      &
             geo_print_lattice,     &
             geo_print_coordinates

  private :: calc_recLattVec,       &
             calc_metric,           &
             geo_init_mbpp,         &
             geo_init_xsf

  !----------- simulation cell, atom types and coordinates ------------!
  ! nAtoms              : total number of atoms                        !
  ! nTypes              : total number of atomic species               !
  ! latticeConstant     : the lattice constant in Bohr units (a.u.)    !
  ! unitCellVolume      : volume of the unit cell in Bohr^3            !
  ! recUnitCellVol      : volume of the reciprocal unit cell           !
  ! lattiveVec(i,ivec)  : i-th component of the ivec-th lattice vector !
  ! nAtomsOfType(itype) : number of atoms of type itype                !
  ! atomType(iatom)     : atomic species number of atom iatom          !
  ! atomTypeName(itype) : name of atomic species itype                 !
  ! cooLattice(i, iat)  : i-th component of the lattice coordinates of !
  !                       the iat-th atom                              !
  ! recLattVec(i,j)     : i-th component of the j-th rec. latt. vector !
  ! latticeVecMet(i,j)  : = sum(latticeVec(1:3,i)*latticeVec(i:3,j))   !
  ! recLattVecMet(i,j)  : = sum(recLattVec(1:3,i)*recLattVec(i:3,j))   !
  !                                                                    !
  !------------------------- k-point sampling -------------------------!
  ! kptMultiple(i)      : number of k-points in dimension i            !
  ! kptShift(1:3)       : shift vecor (rec. coo.) of the k-point mesh  !
  ! useSymmetry         : .true. if only irreducible k-points shall be !
  !                       used                                         !
  !                                                                    !
  !----------------------------- smearing -----------------------------!
  ! metallic            : .true., if the system is metallic            !
  ! smearWidth          : width [Ha] for Gaussian smearing             !
  !                                                                    !
  !----------------------- forces / relaxation ------------------------!
  ! doForces            : if .true. --> calculate forces               !
  ! doGeoOpt            : if .true. --> optimize the geometry          !
  ! optIter             : number of iterations for the optimization    !
  ! optTol              : optimization threshold for the max. force    !
  ! optAlpha            : initial shift proportional to forces         !
  ! optAll              : if .true. --> optimize all atoms             !
  ! optNumActive        : number of active atoms for optimization      !
  ! optActive(i)        : i-th active atom                             !
  ! optNumDoF           : total number of active degrees of freedom    !
  ! optDoF(i,iat)       : status of the i-th dimension of the iat-th   !
  !                       active atom; 1 --> on / 0 --> off            !
  !                                                                    !
  !------------------------ band plot options -------------------------!
  ! doBandPlot          : .true., if band structure plot is desired    !
  ! bandPlotBands       : number of bands in the band structure        !
  ! bandPlotPoints      : number of k-points in the band structure     !
  ! bandPlotPoint(1:3,i): the i-th k-point for the band structure      !
  ! bandPlotStepSize    : step size for interpolation (rez. units)     !
  !                                                                    !
  !------------------------- DOS plot options -------------------------!
  ! doDOSPlot           : .true., if DOS plot shall be calculated      !
  ! DOSPoints           : number of energy points in DOS plot          !
  ! DOSWidth            : smearing width for DOS calculation           !
  ! DOSEmin             : minimum energy, if given                     !
  ! DOSEmax             : maximum energy, if givenm                    !
  !--------------------------------------------------------------------!

  integer,                                       public :: nAtoms
  integer,                                       public :: nTypes
  double precision,                              public :: latticeConstant
  double precision,                              public :: unitCellVolume
  double precision,                              public :: recUnitCellVol
  double precision, dimension(3,3),              public :: latticeVec
  integer,          dimension(:),   allocatable, public :: nAtomsOfType
  integer,          dimension(:),   allocatable, public :: atomType
  character(len=2), dimension(:),   allocatable, public :: atomTypeName
  double precision, dimension(:,:), allocatable, public :: cooLattice
  double precision, dimension(3,3),              public :: recLattVec
  double precision, dimension(3,3),              public :: latticeVecMet
  double precision, dimension(3,3),              public :: recLattVecMet

  integer,          dimension(3),                public :: kptMultiple
  double precision, dimension(3),                public :: kptShift
  logical,                                       public :: useSymmetry

  logical,                                       public :: metallic
  double precision,                              public :: smearWidth

  logical,                                       public :: doForces
  logical,                                       public :: doGeoOpt
  integer,                                       public :: optIter
  double precision,                              public :: optTol
  double precision,                              public :: optAlpha
  logical,                                       public :: optAll
  integer,                                       public :: optNumActive
  integer,          dimension(:),   allocatable, public :: optActive
  integer,                                       public :: optNumDoF
  integer,          dimension(:,:), allocatable, public :: optDoF

  logical,                                       public :: doDOSPlot
  integer,                                       public :: DOSPoints
  double precision,                              public :: DOSWidth
  double precision,                              public :: DOSEmin
  double precision,                              public :: DOSEmax

  logical,                                       public :: doBandPlot
  integer,                                       public :: bandPlotBands
  integer,                                       public :: bandPlotPoints
  double precision, dimension(:,:), allocatable, public :: bandPlotPoint
  double precision,                              public :: bandPlotStepSize

contains

  subroutine geo_init(infile, filetype)

    implicit none

    character(len=*), intent(in) :: infile, filetype

    ! defaults:
    doBandPlot  = .false.
    metallic    = .false.
    useSymmetry = .true.

    select case(trim(io_lower(filetype)))
    case('mbpp')
       call geo_init_mbpp(infile)
    case('xsf')
       call geo_init_xsf(infile)
    case default
       write(0,*) "Error: unknown input file format: ", &
                  trim(adjustl(filetype))
       stop
    end select

  end subroutine geo_init

  !--------------------------------------------------------------------!

  subroutine geo_final()

    implicit none

    if (allocated(atomType))        deallocate(atomType)
    if (allocated(atomTypeName))    deallocate(atomTypeName)
    if (allocated(bandPlotPoint))   deallocate(bandPlotPoint)
    if (allocated(cooLattice))      deallocate(cooLattice)
    if (allocated(nAtomsOfType))    deallocate(nAtomsOfType)
    if (allocated(optActive))       deallocate(optActive, optDoF)
    
  end subroutine geo_final





  !=============================== I/O ================================!

  subroutine geo_print_status()

    implicit none

    write(*,*) 'Input data summary:'
    write(*,*) '==================='
    write(*,*)
    write(*,'(1x,"Number of atoms    : ",A)') trim(io_adjustl(nAtoms))
    write(*,'(1x,"Number of species  : ",A)') trim(io_adjustl(nTypes))
    write(*,'(1x,"Metallic           : ",L1)') metallic
    write(*,'(1x,"Smearing width     : ",A," Ha")') trim(io_adjustl(smearWidth,8))
    write(*,'(1x,"Unit cell volume   : ",A, " Bohr^3")') trim(io_adjustl(unitCellVolume,3))
    write(*,'(1x,"k-point multiples  : ",A,1x,A,1x,A)') trim(io_adjustl(kptMultiple(1))), &
         trim(io_adjustl(kptMultiple(2))), trim(io_adjustl(kptMultiple(3)))
    write(*,'(1x,"k-point grid shift : ",A,1x,A,1x,A)') trim(io_adjustl(kptShift(1),2)), &
         trim(io_adjustl(kptShift(2),2)), trim(io_adjustl(kptShift(3),2))
    write(*,*)
    write(*,*)

  end subroutine geo_print_status

  !--------------------------------------------------------------------!

  subroutine geo_print_lattice()

    implicit none

    write(*,*) 'Lattice vectors:'
    write(*,*)
    write(*,'(4x,3(F15.8,1x))') latticeVec(1:3,1)
    write(*,'(4x,3(F15.8,1x))') latticeVec(1:3,2)
    write(*,'(4x,3(F15.8,1x))') latticeVec(1:3,3)
    write(*,*)

    write(*,*) 'Reciprocal lattice vectors:'
    write(*,*)
    write(*,'(4x,3(F15.8,1x))') recLattVec(1:3,1)
    write(*,'(4x,3(F15.8,1x))') recLattVec(1:3,2)
    write(*,'(4x,3(F15.8,1x))') recLattVec(1:3,3)
    write(*,*)
    write(*,*)

  end subroutine geo_print_lattice

  !--------------------------------------------------------------------!

  subroutine geo_print_coordinates()

    implicit none

    integer :: iatom, itype, iat

    write(*,*) 'Lattice ccoordinates:'
    write(*,*)
    iatom = 0
    do itype = 1, nTypes
    do iat   = 1, nAtomsOfType(itype)
       iatom = iatom + 1
       write(*,'(1x,A2,1x,3(F15.8,1x))') atomTypeName(itype), cooLattice(1:3,iatom)
    end do
    write(*,*)
    end do
    write(*,*)

  end subroutine geo_print_coordinates



  !======================= auxiliary procedures =======================!





  !--------------------------------------------------------------------!
  !           calculation of the reciprocal lattice vectors            !
  !--------------------------------------------------------------------!

  subroutine calc_recLattVec(avec, avol, bvec, bvol)

    implicit none

    double precision, dimension(3,3), intent(in)  :: avec
    double precision,                 intent(in)  :: avol
    double precision, dimension(3,3), intent(out) :: bvec
    double precision,                 intent(out) :: bvol

    bvec(1,1) =  avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3)
    bvec(2,1) =  avec(3,2)*avec(1,3) - avec(1,2)*avec(3,3)
    bvec(3,1) =  avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3)

    bvec(1,2) =  avec(2,3)*avec(3,1) - avec(3,3)*avec(2,1)
    bvec(2,2) =  avec(3,3)*avec(1,1) - avec(1,3)*avec(3,1)
    bvec(3,2) =  avec(1,3)*avec(2,1) - avec(2,3)*avec(1,1)

    bvec(1,3) =  avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)
    bvec(2,3) =  avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)
    bvec(3,3) =  avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)

    bvec(:,:) = bvec(:,:)*PI2/avol

    bvol = bvec(1,1)*bvec(2,2)*bvec(3,3) &
         + bvec(2,1)*bvec(3,2)*bvec(1,3) &
         + bvec(3,1)*bvec(1,2)*bvec(2,3) &
         - bvec(3,1)*bvec(2,2)*bvec(1,3) &
         - bvec(1,1)*bvec(3,2)*bvec(2,3) &
         - bvec(2,1)*bvec(1,2)*bvec(3,3)

  end subroutine calc_recLattVec

  !--------------------------------------------------------------------!
  !       calculation of real space and reciprocal space metrics       !
  !--------------------------------------------------------------------!

  subroutine calc_metric(avec, amet)

    implicit none

    double precision, dimension(3,3), intent(in)  :: avec
    double precision, dimension(3,3), intent(out) :: amet

    integer :: i, j

     do j = 1, 3
     do i = 1, 3
        amet(i,j) = sum(avec(:,i)*avec(:,j))
     end do
     end do

  end subroutine calc_metric




  !======================== input file formats ========================!





  !--------------------------------------------------------------------!
  !            Format: Mixed Basis Pseudopotential Program             !
  !--------------------------------------------------------------------!

  subroutine geo_init_mbpp(infile)
    
  use mbppio, only: mbpp_io_init,            &
                    mbpp_io_final,           &
                    mbpp_get_alat,           &
                    mbpp_get_avec,           &
                    mbpp_get_celvol,         &
                    mbpp_get_coorat,         &
                    mbpp_get_kpmeth,         &
                    mbpp_get_nameat,         &
                    mbpp_get_natoms,         &
                    mbpp_get_symop,          &
                    mbpp_get_nkxyz,          &
                    mbpp_get_ntype,          &
                    mbpp_get_shift,          &
                    mbpp_get_niter,          &
                    mbpp_get_ifmet,          &
                    mbpp_get_intmeth,        &
                    mbpp_get_smearwidth,     &
                    mbpp_get_doforces,       &
                    mbpp_get_geo_niter,      &
                    mbpp_get_geo_relmeth,    &
                    mbpp_get_geo_relspec,    &
                    mbpp_get_geo_alph,       &
                    mbpp_get_geo_ftol,       &
                    mbpp_get_geo_natrel,     &
                    mbpp_get_geo_atoms,      &
                    mbpp_get_dos_nplt,       &
                    mbpp_get_dos_epltmin,    &
                    mbpp_get_dos_epltmax,    &
                    mbpp_get_dos_smearwidth, &
                    mbpp_get_bandplot_nline, &
                    mbpp_get_bandplot_nband, &
                    mbpp_get_bandplot_start, &
                    mbpp_get_bandplot_kpt,   &
                    mbpp_get_bandplot_nstep 

    implicit none

    character(len=*), intent(in) :: infile

    double precision, dimension(3) :: rvec
    double precision               :: rnorm
    character(len=10)              :: relspec

    integer :: itype, iatom, iat, ivec, ipt, i, j

    call mbpp_io_init(file=trim(adjustl(infile)))

    nTypes = mbpp_get_ntype()
    allocate( nAtomsOfType(nTypes), &
              atomTypeName(nTypes)  )

    ! count number of atoms:
    nAtoms = 0
    do itype = 1, nTypes
       nAtomsOfType(itype) = mbpp_get_natoms(itype)
       atomTypeName(itype) = mbpp_get_nameat(itype)
       nAtoms              = nAtoms + nAtomsOfType(itype)
    end do

    allocate(cooLattice(3, nAtoms), &
             atomType(nAtoms)       )
   
    ! atomic coordinates and atom types:
    iatom = 0
    do itype = 1, nTypes
    do iat   = 1, nAtomsOfType(itype)
       iatom = iatom + 1
       atomType(iatom) = itype
       do i     = 1, 3
          cooLattice(i, iatom) = mbpp_get_coorat(i,iat,itype)
       end do
    end do
    end do

    ! unit cell:
    latticeConstant = mbpp_get_alat()
    unitCellVolume  = mbpp_get_celvol()
    do ivec = 1, 3
    do i    = 1, 3
       latticeVec(i,ivec) = mbpp_get_avec(i,ivec)
    end do
    end do

    ! k-point sampling:
    if (trim(io_lower(mbpp_get_kpmeth())) == 'mp') then
       kptMultiple(1:3) = mbpp_get_nkxyz()
       kptShift(1:3)    = mbpp_get_shift()
    else
       write(0,*) 'Error: k-point set method not implemented.'
       kptMultiple(1:3) = (/ 1, 1, 1 /)
       kptShift(1:3)    = (/ 0.0d0, 0.0d0, 0.0d0 /)
    end if

    ! use symmetry for k-points?
    if (trim(mbpp_get_symop()) == 'gen') then
       useSymmetry = .true.
    else
       useSymmetry = .false.
    end if

    ! reciprocal lattice vectors:
    call calc_recLattVec(latticeVec, unitCellVolume, &
                         recLattVec, recUnitCellVol  )

    ! metrics:
    call calc_metric(latticeVec, latticeVecMet)
    call calc_metric(recLattVec, recLattVecMet)

    ! smearing options:
    metallic   = mbpp_get_ifmet()
    smearWidth = mbpp_get_smearwidth()*eV2Ha

    ! forces / geometry optimization:
    doForces = mbpp_get_doforces()
    optIter  = mbpp_get_geo_niter()
    if (optIter > 0) then
       doGeoOpt = .true.
       optTol   = mbpp_get_geo_ftol()
       optAlpha = mbpp_get_geo_alph()
       relspec  = trim(mbpp_get_geo_relspec())
       if (trim(relspec) == 'full') then
          optAll    = .true.
          optNumDoF = 3*nAtoms
       else
          optAll = .false.
          optNumActive = mbpp_get_geo_natrel()
          allocate(optActive(optNumActive), optDoF(3,optNumActive))
          iatom = 0
          j     = 0
          types : do itype = 1, nTypes
          do iat   = 1, nAtomsOfType(itype)
             if (j >= optNumActive) exit types
             iatom = iatom + 1
             do i     = 1, 3
                optDoF(i,j+1) = mbpp_get_geo_atoms(i,iat,itype)
             end do
             if (sum(optDoF(:,j+1))>0) then
                j = j + 1
                optActive(j) = iatom
             end if
          end do
          end do types
          optNumDoF = sum(optDoF(:,:))
          if (j /= optNumActive) then
             write(0,*) "Error: invalid input for geometry optimization."
             write(0,*) optNumActive, j
             stop
          end if
       end if
    end if

    ! DOS plot options:
    DOSPoints = mbpp_get_dos_nplt()
    if (DOSPoints > 0) then
       doDOSPlot = .true.
       DOSEmin   = mbpp_get_dos_epltmin()*eV2Ha
       DOSEmax   = mbpp_get_dos_epltmax()*eV2Ha
       DOSWidth  = mbpp_get_dos_smearwidth()*eV2Ha
    else
       doDOSPlot = .false.
    end if

    ! band plot options:
    bandPlotPoints = mbpp_get_bandplot_nline()
    if (bandPlotPoints > 0) then
       doBandPlot = .true.
       bandPlotPoints = bandPlotPoints + 1
       bandPlotBands  = mbpp_get_bandplot_nband()
       allocate(bandPlotPoint(3,bandPlotPoints))
       bandPlotPoint(1:3,1) = mbpp_get_bandplot_start()
       do ipt = 2, bandPlotPoints
          bandPlotPoint(1:3,ipt) = mbpp_get_bandplot_kpt(ipt-1)
       end do
       rvec(1:3) = bandPlotPoint(1:3,2) - bandPlotPoint(1:3,1)
       rnorm = 0.0d0
       do i = 1, 3
       do j = 1, 3
             rnorm = rnorm + rvec(j)*recLattVecMet(j,i)*rvec(i)
       end do
       end do
       rnorm = sqrt(rnorm)
       bandPlotStepSize = rnorm/dble(mbpp_get_bandplot_nstep(1))
    else
       doBandPlot = .false.
    end if
    
    call mbpp_io_final()
    
  end subroutine geo_init_mbpp

  !--------------------------------------------------------------------!
  !              Format: XCrysDen Structure Format - XSF               !
  !--------------------------------------------------------------------!

  subroutine geo_init_xsf(infile)
    
  use xsflib, only: xsf_init,                &
                    xsf_final,               &
                    xsf_get_alat,            &
                    xsf_get_avec,            &
                    xsf_get_celvol,          &
                    xsf_get_coorat,          &
                    xsf_get_nameat,          &
                    xsf_get_natoms,          &
                    xsf_get_nkxyz,           &
                    xsf_get_ntype,           &
                    xsf_get_shift,           &
                    xsf_get_niter,           &
                    xsf_get_ifmet,           &
                    xsf_get_smearwidth
  
    implicit none

    character(len=*), intent(in)   :: infile

    integer :: itype, iatom, iat, ivec, i

    call xsf_init(file=trim(adjustl(infile)))

    nTypes = xsf_get_ntype()
    allocate( nAtomsOfType(nTypes), &
              atomTypeName(nTypes)  )

    ! count number of atoms:
    nAtoms = 0
    do itype = 1, nTypes
       nAtomsOfType(itype) = xsf_get_natoms(itype)
       atomTypeName(itype) = xsf_get_nameat(itype)
       nAtoms              = nAtoms + nAtomsOfType(itype)
    end do

    allocate(cooLattice(3, nAtoms), atomType(nAtoms) )
   
    ! atomic coordinates and atom types:
    iatom = 0
    do itype = 1, nTypes
    do iat   = 1, nAtomsOfType(itype)
       iatom = iatom + 1
       atomType(iatom) = itype
       do i     = 1, 3
          cooLattice(i, iatom) = xsf_get_coorat(i,iat,itype)*Ang2Bohr
       end do
    end do
    end do

    ! unit cell:
    latticeConstant = 1.0d0
    unitCellVolume  = xsf_get_celvol()*Ang2Bohr*Ang2Bohr*Ang2Bohr
    do ivec = 1, 3
    do i    = 1, 3
       latticeVec(i,ivec) = xsf_get_avec(i,ivec)*Ang2Bohr
    end do
    end do

    ! k-point sampling:
    kptMultiple(1:3) = xsf_get_nkxyz()
    kptShift(1:3)    = xsf_get_shift()

    ! reciprocal lattice vectors:
    call calc_recLattVec(latticeVec, unitCellVolume, &
                         recLattVec, recUnitCellVol  )

    ! convert coordinates to lattice vector units:
    do iatom = 1, nAtoms
       cooLattice(1:3,iatom) = cooLattice(1,iatom)*recLattVec(1,1:3) &
                             + cooLattice(2,iatom)*recLattVec(2,1:3) &
                             + cooLattice(3,iatom)*recLattVec(3,1:3)
       cooLattice(1:3,iatom) = cooLattice(1:3,iatom)/PI2
    end do

    ! metrics:
    call calc_metric(latticeVec, latticeVecMet)
    call calc_metric(recLattVec, recLattVecMet)

    ! smearing options:
    metallic   = xsf_get_ifmet()
    smearWidth = xsf_get_smearwidth()*eV2Ha

    ! band plots not supported:
    doBandPlot = .false.
    
    call xsf_final()
    
  end subroutine geo_init_xsf










end module geometry
