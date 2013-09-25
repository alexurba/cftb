module mbppio

  !--------------------------------------------------------------------!
  ! A library for parsing MBPP file formats                            !
  !--------------------------------------------------------------------!
  ! 2011-01-09 Alexander Urban (AU)                                    !
  ! 2011-04-09 AU -- sections 50 and 55 added (forces/relaxation)      !
  !--------------------------------------------------------------------!

  use io, only: io_adjustl,  &
                io_isfigure, &
                io_isletter, &
                io_readval,  &
                io_lower,    &
                io_upper      

  implicit none

  public  :: defined,                 &
             mbpp_io_init,            &
             mbpp_io_final,           &
             mbpp_parse_output,       &
             mbpp_get_alat,           &
             mbpp_get_AtomsVersion,   &
             mbpp_get_avec,           &
             mbpp_get_BondsVersion,   &
             mbpp_get_celvol,         &
             mbpp_get_coorat,         &
             mbpp_get_nameat,         &
             mbpp_get_natomax,        &
             mbpp_get_natoms,         &
             mbpp_get_symop,          &
             mbpp_get_ntype,          &
             mbpp_get_Rcmax,          &
             mbpp_get_title,          &
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
             mbpp_get_dos_intmeth,    &
             mbpp_get_dos_smearwidth, &
             mbpp_get_dos_nplot,      &
             mbpp_get_dos_totdos,     &
             mbpp_get_dos_filename,   &
             mbpp_get_dos_nrtyp,      &
             mbpp_get_dos_nrat,       &
             mbpp_get_bandplot_nline, &
             mbpp_get_bandplot_nband, &
             mbpp_get_bandplot_start, &
             mbpp_get_bandplot_kpt,   &
             mbpp_get_bandplot_nstep 

  private :: mbpp_readline,           &
             mbpp_read_bandplot,      &
             mbpp_read_coordinates,   &
             mbpp_read_tboptions,     &
             mbpp_parse_coo


  !--------------------------------------------------------------------!
  ! lenline : max. length of lines (number of char's) in input file    !
  ! lentype : max. length of type names (number of char's)             !
  !--------------------------------------------------------------------!
  ! u_in    : unit number of the input file (INP)                      !
  ! u_co    : unit number of the coordinates file (COORAT)             !
  !--------------------------------------------------------------------!

  integer, parameter, private :: lenline   = 132
  integer, parameter, private :: lentype   = 2
  integer, parameter, private :: lenbopfox = 50

  integer, parameter, private :: u_in    = 20
  integer, parameter, private :: u_co    = 30
  integer, parameter, private :: u_out   = 40


  character(len=100), private :: f_co    = 'COORAT'


  !--------------------------------------------------------------------!
  ! title              : title (1st line) of the input file            !
  !--------------------------------------------------------------------!
  ! ntype              : number of atom types / atomic species         !
  ! natomax            : max. number of atoms per atom type            !
  ! alat               : lattice constant (as read from the input file !
  ! celvol             : unit cell volume                              !
  ! avec(i,j)          : i-th component of the j-th lattice vector     !
  ! nameat(ityp)       : name of atomic species ityp; ityp = 1, ntype  !
  ! natoms(ityp)       : number of atoms of atomic species ityp        !
  ! coorat(i,iat,ityp) : i-th component of the lattice coordinates of  !
  !                      the iat-th atom of atomic species ityp        !
  ! forlat(i,iat,ityp) : atomic forces (lattice units), same logic     !
  ! symop              : value of the `symop' flag.                    !
  !                      'gen' --> generate symmetry operations        !
  !                                                                    !
  !------------------------- k-point sampling -------------------------!
  ! kpmeth             : method for generating the k-point set         !
  ! nband              : number of bands to be calculated              !
  ! ifmin              : lowest occupied band                          !
  ! ifmax              : highest possibly occupied band                !
  ! nkxyz(1:3)         : numbers of k-points per dimension             !
  ! shift(1:3)         : shift vector of the k-point mesh              !
  !                                                                    !
  !--------------------------- SCF settings ---------------------------!
  ! niter              : max. number of SCF iterations                 !
  ! ifmet              : .true., if the system is metallic             !
  ! intmeth            : integration method for the DOS integration    !
  ! smearwidth         : width [eV] for smearing of the Fermi edge     !
  !                                                                    !
  !---------------------- geometry optimization -----------------------!
  ! doforces           : if .true. --> calculate forces                !
  ! geo_niter          : number of geometry optiumization iterations   !
  ! geo_relmeth        : relaxation method for the geometry opt.       !
  ! geo_relspec        : specification of the degrees of freedom       !
  ! geo_alph           : initial shift, proportional to forces         !
  ! geo_ftol           : convergence threshold for the max. force      !
  ! geo_natrel         : number of atoms to be optimized               !
  ! geo_atoms(i,iat)   : for the iat-th atom of those geo_natrel atoms !
  !                      i = 1 --> type index; i = 2 --> atom index    !
  !                      i = 3..4 --> (de)activate dimensions          !
  !                                                                    !
  !---------------------------- band plots ----------------------------!
  ! bandplot_nkpts      : number of k-points for a band plot           !
  ! bandplot_nband      : number of bands for the band plot            !
  ! bandplot_start(1:3) : starting point for a band plot               !
  ! bandplot_kpt(1:3,i) : i-th k-point in a band plot                  !
  ! bandplot_nstep(i)   : number of sampling points to k-point i       !
  !                                                                    !
  !---------------------------- DOS plots -----------------------------!
  ! dos_nplt          : number of data points for the DOS plot         !
  ! dos_epltmin       : initial energy [eV] relative to Fermi level    !
  ! dos_epltmax       : final energy [eV] relative to Fermi level      !
  ! dos_intmeth       : integration method for DOS integratio          !
  ! dos_smearwidth    : broadening width [eV] for smearing             !
  ! dos_nplot         : number of pDOS plots                           !
  ! dod_totdos        : .true., if the total DOS shall be computed     !
  ! dos_filename(i)   : file name for the i-th pDOS plot               !
  ! dos_nrtyp(i)      : atom type number for the i-th pDOS plot        !
  ! dos_nrat(i)       : atom index for the i-th pDOS plot              !
  !                                                                    !
  !-------------------------- tight binding ---------------------------!
  ! AtomsVersion       : BOPfox tight-binding atoms version            !
  ! BondsVersion       : BOPfox tight-binding bonds version            !
  ! Rcmax              : tight-binding max. radial cut-off             !
  !                                                                    !
  !----------------------------- energies -----------------------------!
  ! E_tot              : total energy                                  !
  ! E_band             : band sum energy                               !
  ! E_vac              : average vacuum energy level                   !
  !--------------------------------------------------------------------!

  character(len=lenline),                                public  :: title

  integer,                                               public  :: ntype
  integer,                                               public  :: natomax
  double precision,                                      public  :: alat
  double precision,                                      public  :: celvol
  double precision,       dimension(3,3),                public  :: avec
  character(len=lentype), dimension(:),     allocatable, public  :: nameat
  integer,                dimension(:),     allocatable, public  :: natoms
  double precision,       dimension(:,:,:), allocatable, public  :: coorat
  double precision,       dimension(:,:,:), allocatable, public  :: forlat
  character(len=lenline),                                public  :: symop

  character(len=lenline),                                public  :: kpmeth
  integer,                                               public  :: nband 
  integer,                                               public  :: ifmin 
  integer,                                               public  :: ifmax 
  integer,                dimension(3),                  public  :: nkxyz 
  double precision,       dimension(3),                  public  :: shift 

  integer,                                               public  :: niter     
  logical,                                               public  :: ifmet     
  character(len=lenline),                                public  :: intmeth   
  double precision,                                      public  :: smearwidth

  logical,                                               public  :: doforces   
  integer,                                               public  :: geo_niter  
  character(len=lenline),                                public  :: geo_relmeth
  character(len=lenline),                                public  :: geo_relspec
  double precision,                                      public  :: geo_alph   
  double precision,                                      public  :: geo_ftol
  integer,                                               public  :: geo_natrel
  integer,                dimension(:,:),   allocatable, public  :: geo_atoms
  
  integer,                                               public  :: bandplot_nline
  integer,                                               public  :: bandplot_nband
  double precision,       dimension(3),                  public  :: bandplot_start
  double precision,       dimension(:,:),   allocatable, public  :: bandplot_kpt
  integer,                dimension(:),     allocatable, public  :: bandplot_nstep

  integer,                                               public  :: dos_nplt
  double precision,                                      public  :: dos_epltmin
  double precision,                                      public  :: dos_epltmax
  character(len=lenline),                                public  :: dos_intmeth
  double precision,                                      public  :: dos_smearwidth
  integer,                                               public  :: dos_nplot
  logical,                                               public  :: dos_totdos
  character(len=lenline), dimension(:),     allocatable, public  :: dos_filename
  integer,                dimension(:),     allocatable, public  :: dos_nrtyp
  integer,                dimension(:),     allocatable, public  :: dos_nrat

  character(len=lenbopfox),                              public  :: AtomsVersion
  character(len=lenbopfox),                              public  :: BondsVersion
  double precision,                                      public  :: Rcmax

  double precision,                                      public  :: E_tot
  double precision,                                      public  :: E_band
  double precision,                                      public  :: E_vac

  !--------------------------------------------------------------------!
  !                              PRIVATE                               !
  !                                                                    !
  ! lenelem    : max. length of `#' statements                         !
  ! nelsmax    : max. number of `#' definitions                        !
  ! nels       : number of `#' definitions set                         !
  ! element(i) : i-th `#' definition from input file                   !
  !--------------------------------------------------------------------!
  ! iline   : number of the current line of the input file             !
  !--------------------------------------------------------------------!

  integer,                                    parameter, private :: lenelem = 30
  integer,                                    parameter, private :: nelsmax = 50
  integer,                                               private :: nels    =  0
  character(len=lenelem), dimension(nelsmax),            private :: element

  integer,                                               private :: iline

contains

  !--------------------------------------------------------------------!
  !                     parse the MBPP input file                      !
  !--------------------------------------------------------------------!

  subroutine mbpp_io_init(file, coofile)

    implicit none

    character(len=*), optional, intent(in) :: file
    character(len=*), optional, intent(in) :: coofile

    !------------------------------------------------------------------!
    ! fname    : file name of the input file                           !
    ! line     : line buffer (single line from input file)             !
    ! fexist   : .true. if the input file was found                    !
    ! iblock   : last read section number                              !
    !------------------------------------------------------------------!

    character(len=50)                      :: fname
    character(len=lenline)                 :: line
    logical                                :: fexists
    integer                                :: iblock

    ! set input file name:
    if (present(file)) then
       fname = trim(adjustl(file))
    else
       fname = 'INP'
    end if

    ! set coordinates file name:
    if (present(coofile)) then
       f_co = trim(adjustl(coofile))
    end if

    ! stop, if input file does not exist:
    inquire(file=trim(fname), exist=fexists)
    if (.not. fexists) then
       write(0,*) 'Error: File not found (MBPP input): ', trim(fname)
       stop
    end if

    ! defaults:
    ifmet     = .false.
    niter     = 0
    doforces  = .false.
    dos_nplt  = 0
    geo_niter = 0

    open(u_in, file=trim(fname), status='old', action='read')
    rewind(u_in)

    ! input file title/comment:
    read(u_in, '(A)') title

    iline = 0
    readINP : do

       call mbpp_readline(line)

       ! new block ID found?
       if ((len_trim(line)==2) .and. io_isfigure(line(1:1)) &
                               .and. io_isfigure(line(2:2)) ) then
          read(line, *) iblock

          block : select case(iblock)
          case(00) ! end of input file
             exit readINP
          case(10) ! coordinates section
             call mbpp_read_coordinates()
          case(25) ! tight-binding options
             call mbpp_read_tboptions()
          case(30) ! k-point sampling
             call mbpp_read_kpoint_options()
          case(40,45) ! scf and smearing
             call mbpp_read_scf()
          case(50) ! calculate forces
             doforces = .true.
          case(55) ! geometry optimization
             doforces = .true.
             call mbpp_read_geoopt()
          case(70) ! DOS plot
             call mbpp_read_dosplot()
          case(72) ! band plot information
             call mbpp_read_bandplot()
          case default ! unknown section -> read on.
             continue
          end select block

       end if

    end do readINP

    close(u_in)

  end subroutine mbpp_io_init

  !--------------------------------------------------------------------!
  !                         deallocate memory                          !
  !--------------------------------------------------------------------!

  subroutine mbpp_io_final()

    implicit none

    if (allocated(coorat)) then
       deallocate(coorat, nameat, natoms)
    end if

    if (allocated(forlat))    deallocate(forlat)
    if (allocated(geo_atoms)) deallocate(geo_atoms)

    if (allocated(bandplot_kpt)) then
       deallocate(bandplot_kpt, bandplot_nstep)
    end if

    if (allocated(dos_filename)) then
       deallocate(dos_filename, dos_nrtyp, dos_nrat)
    end if

  end subroutine mbpp_io_final

  !--------------------------------------------------------------------!
  !                  check if a string was `#define'd                  !
  !--------------------------------------------------------------------!

  function defined(str) result(isdef)

    implicit none

    character(len=*), intent(in) :: str

    logical :: isdef
    integer :: iels

    isdef = .false.
    checkdef : do iels = 1, nels
       if ( trim(adjustl(str)) == trim(adjustl(element(iels))) ) then
          isdef = .true.
          exit checkdef
       end if
    end do checkdef
    
  end function defined

  !--------------------------------------------------------------------!
  !                         output file parser                         !
  !--------------------------------------------------------------------!

  subroutine mbpp_parse_output(outfile)

    implicit none
  
    character(len=*), intent(in) :: outfile

    logical                      :: fexists
    character(len=256)           :: line
    integer                      :: i, iat, itype, iline
    integer                      :: status
    character(len=50)            :: str1, str2

    inquire(file=trim(outfile), exist=fexists)
    if (.not. fexists) then
       write(0,*) 'Error: File not found: ', trim(outfile)
       return
    end if

    E_tot  = 0.0d0
    E_band = 0.0d0
    E_vac  = 0.0d0

    open(u_out, file=trim(outfile), status='old', action='read')

    read(u_out, *)
    read(u_out, '(A)') title
    read(u_out, *)
   
    nels  = 0
    iline = 3
    parse : do
       read(u_out, '(A)', iostat=status) line
       if (status /= 0) exit parse
       iline = iline + 1
       line  = trim(adjustl(line))

       ! define statements:
       i = index(line, '<---- defining')
       if (i /= 0) then         
          if (nels >= nelsmax) then
             write(0,*) "Warning: Too many `#define' statements."
          else
             nels = nels + 1
             read(line(i+14:256), *) element(nels)
          end if
          cycle parse
       end if

       ! cell volume:
       i = index(line, 'Cell volume =')
       if (i /= 0) then         
          read(line(i+13:256), *) celvol
          cycle parse
       end if

       ! lattice vectors:
       i = index(line, 'lattice vectors (a.u.)')
       if (i /= 0) then
          read(u_out,*)
          read(u_out, *) str1, avec(1:3,1)
          read(u_out, *) str1, avec(1:3,2)
          read(u_out, *) str1, avec(1:3,3)
          iline = iline + 4
          cycle parse
       end if      

       ! atomic coordinates:
       i = index(line, 'Atomic Positions:')
       if (i /= 0) then
          read(u_out,*)
          iline = iline + 1
          ! count atoms
          ntype   = 0
          natomax = 0
          str1    = ''
          i       = 0
          countat : do
             read(u_out, '(A)') line
             if (len_trim(line) == 0) cycle countat
             if (scan(line,'#') == 0) exit countat
             read(line, *) str2
             if (trim(str2) /= trim(str1)) then
                i     = 0
                str1  = str2
                ntype = ntype + 1
             end if
             i       = i + 1
             natomax = max(natomax, i)
          end do countat
          allocate(nameat(ntype), natoms(ntype), coorat(3,natomax,ntype), &
                   forlat(3,natomax,ntype))
          rewind(u_out)
          do i = 1, iline
             read(u_out, *)
          end do

          read(u_out, '(A)') line
          iline = iline + 1
          itype = 1
          natoms(itype) = 1
          read(line, *) nameat(itype), str1, str2, coorat(1:3, natoms(itype), itype)
          i = scan(line, ':') + 1
          readcoo : do
             read(u_out, '(A)') line
             iline = iline + 1
             if (len_trim(line) == 0) cycle readcoo
             if (scan(line,'#') == 0) exit readcoo
             read(line, *) str1
             if (trim(str1) == trim(nameat(itype))) then
                natoms(itype) = natoms(itype) + 1
             else
                itype         = itype + 1
                natoms(itype) = 1
                nameat(itype) = trim(str1)
             end if
             read(line(i:256), *) coorat(1:3, natoms(itype), itype)
          end do readcoo

          cycle parse
       end if      

       ! k-points:
       i = index(line, 'nkabc =')
       if (i /= 0) then
          call io_readval(line, name='nkabc', value=nkxyz, n=3)
          call io_readval(line, name='shift', value=shift, n=3)
          cycle parse
       end if      

       ! number of iterations:
       i = index(line, 'max number of iterations:')
       if (i /= 0) then
          read(line(i+26:256),*) niter
          cycle parse
       end if      

       ! metal / insulator:
       i = index(line, 'insulator/metal:')
       if (i /= 0) then
          read(line(i+17:256),*) str1
          i = scan(str1,'m')
          if (i > 0) then
             ifmet = .true.
          else
             ifmet = .false.
          end if
          cycle parse
       end if      

       ! smearing width:
       i = index(line, 'BZ integration method:')
       if (i /= 0) then
          read(line(i+23:256),*) str1
          i = scan(str1,'-')
          read(str1(1:i-1), *) smearwidth
          cycle parse
       end if      

       ! band energy:
       i = index(line, 'band energy          =')
       if (i /= 0) then
          read(line(i+23:256),*) E_band
          cycle parse
       end if      

       ! total energy:
       i = index(line, 'total energy         =')
       if (i /= 0) then
          read(line(i+23:256),*) E_tot
          cycle parse
       end if      

       ! vacuum energy:
       i = index(line, 'Vacuum level (Ry)')
       if (i /= 0) then
          i = index(line, 'av:')
          if (i /= 0) then
             read(line(i+4:256),*) E_vac
          else
             i = index(line, '0.5 :')
             read(line(i+6:256),*) E_vac
          end if
          cycle parse
       end if      

       ! atomic forces (lattice units):
       i = index(line, 'Force/Stress analysis :')
       if (i /= 0) then
          do
             read(u_out, '(A)') line
             iline = iline + 1
             i = index(line, 'force(latt.)')
             if (i>0) exit
          end do
          read(u_out, '(A)')
          iline = iline + 1
          do itype = 1, ntype
             do iat = 1, natoms(itype)
                read(u_out, *) str1, str1, str1, str1, str1, str1, &
                               forlat(1:3, iat, itype)
             end do
          end do
          cycle parse
       end if

    end do parse

    close(u_out)

    ! calculate lattice constant:
    alat = avec(1,1)*avec(2,2)*avec(3,3) &
         + avec(2,1)*avec(3,2)*avec(1,3) &
         + avec(3,1)*avec(1,2)*avec(2,3) &
         - avec(1,3)*avec(2,2)*avec(3,1) &
         - avec(2,3)*avec(3,2)*avec(1,1) &
         - avec(3,3)*avec(1,2)*avec(2,1)
    alat = (celvol/alat)**(1.0d0/3.0d0)
  
  end subroutine mbpp_parse_output

  !--------------------------------------------------------------------!
  !                data access interfaces / properties                 !
  !--------------------------------------------------------------------!

  function mbpp_get_title() result(out_title)
    implicit none
    character(len=lenline) :: out_title
    out_title = title
  end function mbpp_get_title

  !--------------------------------------------------------------------!

  function mbpp_get_ntype() result(out_ntype)
    implicit none
    integer :: out_ntype
    out_ntype = ntype
  end function mbpp_get_ntype

  !--------------------------------------------------------------------!

  function mbpp_get_natomax() result(out_natomax)
    implicit none
    integer :: out_natomax
    out_natomax = natomax
  end function mbpp_get_natomax

  !--------------------------------------------------------------------!

  function mbpp_get_alat() result(out_alat)
    implicit none
    double precision :: out_alat
    out_alat = alat
  end function mbpp_get_alat

  !--------------------------------------------------------------------!

  function mbpp_get_celvol() result(out_celvol)
    implicit none
    double precision :: out_celvol
    out_celvol = celvol
  end function mbpp_get_celvol

  !--------------------------------------------------------------------!

  function mbpp_get_avec(i, ivec) result(out_avec)
    implicit none
    integer, intent(in) :: i, ivec
    double precision    :: out_avec
    out_avec = avec(i,ivec)
  end function mbpp_get_avec

  !--------------------------------------------------------------------!

  function mbpp_get_natoms(itype) result(out_natoms)
    implicit none
    integer, intent(in) :: itype
    integer             :: out_natoms
    if (allocated(natoms)) then
       out_natoms = natoms(itype)
    else
       out_natoms = 0
    end if
  end function mbpp_get_natoms

  !--------------------------------------------------------------------!

  function mbpp_get_nameat(itype) result(out_nameat)
    implicit none
    integer,    intent(in) :: itype
    character(len=lentype) :: out_nameat
    if (allocated(nameat)) then
       out_nameat = nameat(itype)
    else
       out_nameat = ' '
    end if
  end function mbpp_get_nameat

  !--------------------------------------------------------------------!

  function mbpp_get_coorat(i, iatom, itype) result(out_coorat)
    implicit none
    integer, intent(in) :: i, iatom, itype
    double precision    :: out_coorat
    if (allocated(coorat)) then
       out_coorat = coorat(i, iatom, itype)
    else
       out_coorat = 0.0d0
    end if
  end function mbpp_get_coorat

  !--------------------------------------------------------------------!

  function mbpp_get_symop() result(out_symop)
    implicit none
    character(len=lenline) :: out_symop
    if (allocated(coorat)) then
       out_symop = trim(symop)
    else 
       out_symop = 'none'
    end if
  end function mbpp_get_symop

  !--------------------------------------------------------------------!

  function mbpp_get_kpmeth() result(out_kpmeth)
    implicit none
    character(len=lenline) :: out_kpmeth
    out_kpmeth = kpmeth
  end function mbpp_get_kpmeth

  !--------------------------------------------------------------------!

  function mbpp_get_nband() result(out_nband)
    implicit none
    integer :: out_nband
    out_nband = nband
  end function mbpp_get_nband

  !--------------------------------------------------------------------!
  
  function mbpp_get_ifmin() result(out_ifmin)
    implicit none
    integer :: out_ifmin
    out_ifmin = ifmin
  end function mbpp_get_ifmin

  !--------------------------------------------------------------------!

  function mbpp_get_ifmax() result(out_ifmax)
    implicit none
    integer :: out_ifmax
    out_ifmax = ifmax
  end function mbpp_get_ifmax

  !--------------------------------------------------------------------!

  function mbpp_get_nkxyz() result(out_nkxyz)
    implicit none
    integer, dimension(3) :: out_nkxyz
    out_nkxyz = nkxyz
  end function mbpp_get_nkxyz

  !--------------------------------------------------------------------!

  function mbpp_get_shift() result(out_shift)
    implicit none
    double precision, dimension(3) :: out_shift
    out_shift(1:3) = shift
  end function mbpp_get_shift
  
  !--------------------------------------------------------------------!

  function mbpp_get_AtomsVersion() result(out_AtomsVersion)
    implicit none
    character(len=lenbopfox) :: out_AtomsVersion
    out_AtomsVersion = AtomsVersion
  end function mbpp_get_AtomsVersion

  !--------------------------------------------------------------------!

  function mbpp_get_BondsVersion() result(out_BondsVersion)
    implicit none
    character(len=lenbopfox) :: out_BondsVersion
    out_BondsVersion = BondsVersion
  end function mbpp_get_BondsVersion

  !--------------------------------------------------------------------!

  function mbpp_get_Rcmax() result(out_Rcmax)
    implicit none
    double precision :: out_Rcmax
    out_Rcmax = Rcmax
  end function mbpp_get_Rcmax

  !--------------------------------------------------------------------!

  function mbpp_get_niter() result(out_niter)
    implicit none
    integer :: out_niter
    out_niter = niter
  end function mbpp_get_niter

  !--------------------------------------------------------------------!

  function mbpp_get_ifmet() result(out_ifmet)
    implicit none
    logical :: out_ifmet
    out_ifmet = ifmet
  end function mbpp_get_ifmet

  !--------------------------------------------------------------------!

  function mbpp_get_intmeth() result(out_intmeth)
    implicit none
    character(len=lenline) :: out_intmeth
    out_intmeth = trim(intmeth)
  end function mbpp_get_intmeth

  !--------------------------------------------------------------------!

  function mbpp_get_smearwidth() result(out_width)
    implicit none
    double precision :: out_width
    out_width = smearwidth
  end function mbpp_get_smearwidth

  !--------------------------------------------------------------------!

  function mbpp_get_doforces() result(out_doforces)
    implicit none
    logical :: out_doforces
    out_doforces = doforces
  end function mbpp_get_doforces

  !--------------------------------------------------------------------!

  function mbpp_get_geo_niter() result(out)
    implicit none
    integer :: out
    out = geo_niter
  end function mbpp_get_geo_niter

  !--------------------------------------------------------------------!

  function mbpp_get_geo_relmeth() result(out)
    implicit none
    character(len=lenline) :: out
    out = geo_relmeth
  end function mbpp_get_geo_relmeth

  !--------------------------------------------------------------------!

  function mbpp_get_geo_relspec() result(out)
    implicit none
    character(len=lenline) :: out
    out = geo_relspec
  end function mbpp_get_geo_relspec

  !--------------------------------------------------------------------!

  function mbpp_get_geo_alph() result(out)
    implicit none
    double precision :: out
    out = geo_alph
  end function mbpp_get_geo_alph

  !--------------------------------------------------------------------!

  function mbpp_get_geo_ftol() result(out)
    implicit none
    double precision :: out
    out = geo_ftol
  end function mbpp_get_geo_ftol

  !--------------------------------------------------------------------!

  function mbpp_get_geo_natrel() result(out)
    implicit none
    integer :: out
    out = geo_natrel
  end function mbpp_get_geo_natrel

  !--------------------------------------------------------------------!

  function mbpp_get_geo_atoms(idim, iatom, itype) result(out)
    implicit none
    integer, intent(in) :: idim, iatom, itype
    integer :: out
    integer :: i
    if (.not. allocated(geo_atoms)) then
       out = 1
       return
    else
       out = 0
       do i = 1, geo_natrel
          if ((geo_atoms(1,i)==itype) .and. (geo_atoms(2,i)==iatom)) then
             out = geo_atoms(idim+2, i)
             exit
          end if
       end do
    end if
  end function mbpp_get_geo_atoms

  !--------------------------------------------------------------------!

  function mbpp_get_dos_nplt() result(nplt)
    implicit none
    integer :: nplt
    nplt = dos_nplt
  end function mbpp_get_dos_nplt

  !--------------------------------------------------------------------!

  function mbpp_get_dos_epltmin() result(epltmin)
    implicit none
    double precision :: epltmin
    if (dos_nplt > 0) then
       epltmin = dos_epltmin
    else
       epltmin = 0.0d0
    end if
  end function mbpp_get_dos_epltmin

  !--------------------------------------------------------------------!

  function mbpp_get_dos_epltmax() result(epltmax)
    implicit none
    double precision :: epltmax
    if (dos_nplt > 0) then
       epltmax = dos_epltmax
    else
       epltmax = 0.0d0
    end if
  end function mbpp_get_dos_epltmax

  !--------------------------------------------------------------------!

  function mbpp_get_dos_intmeth() result(intmeth)
    implicit none
    character(len=lenline) :: intmeth
    if (dos_nplt > 0) then
       intmeth = dos_intmeth
    else
       intmeth = 'gauss'
    end if
  end function mbpp_get_dos_intmeth

  !--------------------------------------------------------------------!

  function mbpp_get_dos_smearwidth() result(smearwidth)
    implicit none
    double precision :: smearwidth
    if (dos_nplt > 0) then
       smearwidth = dos_smearwidth
    else
       smearwidth = 0.0d0
    end if
  end function mbpp_get_dos_smearwidth

  !--------------------------------------------------------------------!

  function mbpp_get_dos_nplot() result(nplot)
    implicit none
    integer :: nplot
    if (dos_nplt > 0) then
       nplot = dos_nplot
    else
       nplot = 0
    end if
  end function mbpp_get_dos_nplot

  !--------------------------------------------------------------------!

  function mbpp_get_dos_totdos() result(totdos)
    implicit none
    logical :: totdos
    if (dos_nplt > 0) then
       totdos = dos_totdos
    else
       totdos = .true.
    end if
  end function mbpp_get_dos_totdos

  !--------------------------------------------------------------------!

  function mbpp_get_dos_filename(iplot) result(filename)
    implicit none
    integer, intent(in)    :: iplot
    character(len=lenline) :: filename
    if ((allocated(dos_filename)) .and. (iplot <= dos_nplot)) then
       filename = dos_filename(iplot)
    else
       filename = ''
    end if
  end function mbpp_get_dos_filename

  !--------------------------------------------------------------------!

  function mbpp_get_dos_nrtyp(iplot) result(nrtyp)
    implicit none
    integer, intent(in)    :: iplot
    integer                :: nrtyp
    if ((allocated(dos_nrtyp)) .and. (iplot <= dos_nplot)) then
       nrtyp = dos_nrtyp(iplot)
    else
       nrtyp = 0
    end if
  end function mbpp_get_dos_nrtyp

  !--------------------------------------------------------------------!

  function mbpp_get_dos_nrat(iplot) result(nrat)
    implicit none
    integer, intent(in)    :: iplot
    integer                :: nrat
    if ((allocated(dos_nrat)) .and. (iplot <= dos_nplot)) then
       nrat = dos_nrat(iplot)
    else
       nrat = 0
    end if
  end function mbpp_get_dos_nrat

  !--------------------------------------------------------------------!

  function mbpp_get_bandplot_nline() result(out_nline)
    implicit none
    integer :: out_nline
    if (allocated(bandplot_kpt)) then
       out_nline = bandplot_nline
    else
       out_nline = 0
    end if
  end function mbpp_get_bandplot_nline

  !--------------------------------------------------------------------!

  function mbpp_get_bandplot_nband() result(out_nband)
    implicit none
    integer :: out_nband
    if (allocated(bandplot_kpt)) then
       out_nband = bandplot_nband
    else
       out_nband = 0
    end if
  end function mbpp_get_bandplot_nband

  !--------------------------------------------------------------------!

  function mbpp_get_bandplot_start() result(out_start)
    implicit none
    double precision, dimension(3) :: out_start
    if (allocated(bandplot_kpt)) then
       out_start(1:3) = bandplot_start(1:3)
    else
       out_start(1:3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
    end if
  end function mbpp_get_bandplot_start

  !--------------------------------------------------------------------!

  function mbpp_get_bandplot_kpt(i) result(out_kpt)
    implicit none
    integer,            intent(in) :: i
    double precision, dimension(3) :: out_kpt
    if (allocated(bandplot_kpt)) then
       out_kpt(1:3) = bandplot_kpt(1:3, i)
    else
       out_kpt(1:3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
    end if
  end function mbpp_get_bandplot_kpt

  !--------------------------------------------------------------------!

  function mbpp_get_bandplot_nstep(i) result(out_nstep)
    implicit none
    integer, intent(in) :: i
    integer             :: out_nstep
    if (allocated(bandplot_kpt)) then
       out_nstep = bandplot_nstep(i)
    else
       out_nstep = 0
    end if
  end function mbpp_get_bandplot_nstep

  !============================= PRIVATE ==============================!



  !--------------------------------------------------------------------!
  !                   read next line from input file                   !
  !--------------------------------------------------------------------!

  subroutine mbpp_readline(line, unit)

    implicit none

    character(len=lenline), intent(out) :: line
    integer, optional,      intent(in)  :: unit

    logical                             :: fopened
    integer                             :: u_read
    integer                             :: status

    if (present(unit)) then
       u_read = unit
    else
       u_read = u_in
    end if

    inquire(unit=u_read, opened=fopened)
    if (.not. fopened) then
       write(0,*) 'Error: No input file open in MBPP_read_line().'
       stop
    end if

    readline : do

       iline = iline + 1
       read(u_read, '(A)', iostat=status) line
       if (status /= 0) then
          write(0,*) 'Error: End of file during read in MBPP_read_line().'
          stop
       end if

       line = trim(adjustl(line))

       ! comments and empty lines:
       if (line(1:1) == '!')    cycle readline
       if (len_trim(line) == 0) cycle readline

       ! `#' statements:
       if (line(1:1) == '#') then
          if (nels >= nelsmax) then
             write(0,*) "Warning: Too many `#define' statements. Ignoring line ", &
                        io_adjustl(iline)
             write(0,*) trim(line)
             cycle readline
          end if
          nels = nels + 1
          read(line(8:lenline), *) element(nels)
          cycle readline
       end if

       exit readline

    end do readline

  end subroutine mbpp_readline

  !--------------------------------------------------------------------!
  !               read coordinates section of input file               !
  !--------------------------------------------------------------------!

  subroutine mbpp_read_coordinates()

    implicit none

    character(len=lenline) :: line

    character(len=20)      :: source
    integer                :: coo_unit
    logical                :: fexists
    double precision       :: vol
    integer                :: ityp, iat, ipos, i

    call mbpp_readline(line)

    call io_readval(line, name='ntype',   value=ntype)
    call io_readval(line, name='natomax', value=natomax)

    allocate(nameat(ntype), natoms(ntype), coorat(3,natomax,ntype))

    ipos = index(line, 'source')
    if (ipos > 0) then
       call io_readval(line, name='source', value=source)
    else
       source='inp'
    end if

    ! open external coordinates file, if necessary:
    if (trim(source) == 'file') then
       inquire(file=trim(f_co), exist=fexists)
       if (.not. fexists) then
          write(0,*) 'Error: Coordinates file '// trim(f_co) //' not found.'
          call mbpp_io_final()
          stop
       end if
       open(u_co, file=trim(f_co), status='old', action='read')
       coo_unit = u_co
    else
       coo_unit = u_in
    end if

    ! read atomic coordinates:
    do ityp = 1, ntype
       call mbpp_readline(line, unit=coo_unit)
       call io_readval(line, name='natom', value=natoms(ityp))
       call io_readval(line, name='name',  value=nameat(ityp))
       do iat = 1, natoms(ityp)
          call mbpp_readline(line, unit=coo_unit)
          coorat(:,iat,ityp) = mbpp_parse_coo(line)
       end do
    end do
    
    ! close coordinates file, if it was opened:
    if (trim(source) == 'file') then
       close(u_co)
    end if

    ! lattice constant or unit cell volume:
    call mbpp_readline(line)
    ipos = index(line, 'alat')
    if (ipos > 0) then
       call io_readval(line, name='alat', value=alat)
       vol  = -1.0d0
    else
       call io_readval(line, name='vol',  value=celvol)
       alat = -1.0d0
    end if
    
    ! lattice vectors:
    do i = 1, 3
       call mbpp_readline(line)
       avec(:,i) = mbpp_parse_coo(line)
    end do

    vol = avec(1,1)*avec(2,2)*avec(3,3) &
        + avec(2,1)*avec(3,2)*avec(1,3) &
        + avec(3,1)*avec(1,2)*avec(2,3) &
        - avec(3,1)*avec(2,2)*avec(1,3) &
        - avec(2,1)*avec(1,2)*avec(3,3) &
        - avec(1,1)*avec(3,2)*avec(2,3)

    if (vol < 0) then
       write(0,*) 'Error: Lattice vectors define lefthanded coordinate system.'
       call mbpp_io_final()
       stop
    end if

    ! rescale lattice vectors with respect to the lattice constant:
    if (alat == -1.0d0) then
       alat = (celvol/vol)**(1.0d0/3.0d0)
    else
       celvol = vol*alat*alat*alat
    end if
    avec(:,:) = avec(:,:)*alat

    ! read first line of symmetry settings:
    call mbpp_readline(line)
    call io_readval(line, name='symop', value=symop)

    ! Any further lines in this section (e.g. symmetry operations) are
    ! ignored.

  end subroutine mbpp_read_coordinates

  !--------------------------------------------------------------------!
  !                   read k-point sampling options                    !
  !--------------------------------------------------------------------!

  subroutine mbpp_read_kpoint_options()

    implicit none

    character(len=lenline) :: line

    call mbpp_readline(line)
    call io_readval(line, name='kpmeth', value=kpmeth)
    call mbpp_readline(line)
    call io_readval(line, name='nband',  value=nband)
    call io_readval(line, name='ifmin',  value=ifmin)
    call io_readval(line, name='ifmax',  value=ifmax)
    call io_readval(line, name='nkxyz',  value=nkxyz, n=3)
    call io_readval(line, name='shift',  value=shift, n=3)

  end subroutine mbpp_read_kpoint_options

  !--------------------------------------------------------------------!
  !                  read SCF and smearing parameters                  !
  !--------------------------------------------------------------------!

  subroutine mbpp_read_scf()

    implicit none

    character(len=lenline) :: line, metal

    call mbpp_readline(line)
    call io_readval(line, name='niter',   value=niter)
    call io_readval(line, name='ifmet',   value=metal)
    call io_readval(line, name='intmeth', value=intmeth)
    call io_readval(line, name='width',   value=smearwidth)

    if (trim(io_lower(metal)) == 'yes') then
       ifmet = .true.
    else
       ifmet = .false.
    end if

  end subroutine mbpp_read_scf

  !--------------------------------------------------------------------!
  !                       geometry optimization                        !
  !--------------------------------------------------------------------!

  subroutine mbpp_read_geoopt()

    implicit none

    character(len=lenline) :: line
    integer                :: i
    character(len=3)       :: yn

    call mbpp_readline(line)
    geo_niter   = 0
    geo_relmeth = 'bfgs'
    geo_relspec = 'full'
    geo_alph    = 0.1d0
    geo_ftol    = 1.0d-4
    geo_natrel  = 0
    call io_readval(line, name='niter',   value=geo_niter)
    call io_readval(line, name='relmeth', value=geo_relmeth)
    call io_readval(line, name='relspec', value=geo_relspec)
    call io_readval(line, name='alph',    value=geo_alph)
    call io_readval(line, name='ftol',    value=geo_ftol)
    call io_readval(line, name='natrel',  value=geo_natrel)

    if ((geo_relspec /= 'full') .and. (geo_natrel > 0)) then
       allocate(geo_atoms(5,geo_natrel))
       do i = 1, geo_natrel
          geo_atoms(:,i) = 1
          call mbpp_readline(line)
          call io_readval(line, name='nrtyp',   value=geo_atoms(1,i))
          call io_readval(line, name='nrat',    value=geo_atoms(2,i))
          yn = 'yes'
          call io_readval(line, name='a1_axis', value=yn)
          if (yn == 'no') geo_atoms(3,i) = 0
          yn = 'yes'
          call io_readval(line, name='a2_axis', value=yn)
          if (yn == 'no') geo_atoms(4,i) = 0
          yn = 'yes'
          call io_readval(line, name='a3_axis', value=yn)
          if (yn == 'no') geo_atoms(5,i) = 0
       end do
    end if

  end subroutine mbpp_read_geoopt

  !--------------------------------------------------------------------!
  !                       read DOS plot options                        !
  !--------------------------------------------------------------------!

  subroutine mbpp_read_dosplot()

    implicit none

    character(len=lenline) :: line, str

    dos_nplt       = 5000
    dos_epltmin    = 0.0d0
    dos_epltmax    = 0.0d0
    dos_intmeth    = 'gauss'
    dos_smearwidth = 0.4d0
    dos_nplot      = 0
    dos_totdos     = .true.

    call mbpp_readline(line)
    call io_readval(line, name='nplt',    value=dos_nplt)
    call io_readval(line, name='epltmin', value=dos_epltmin)
    call io_readval(line, name='epltmax', value=dos_epltmax)
    call io_readval(line, name='intmeth', value=dos_intmeth)
    call io_readval(line, name='width',   value=dos_smearwidth)

    call mbpp_readline(line)
    call io_readval(line, name='nplot',   value=dos_nplot)
    call io_readval(line, name='totdos',  value=str)
    if (trim(str) /= 'yes') dos_totdos = .false.
    
    if (dos_nplot > 0) then
       allocate( dos_filename(dos_nplot), &
                 dos_nrtyp(dos_nplot),    &
                 dos_nrat(dos_nplot)      )
    end if

  end subroutine mbpp_read_dosplot

  !--------------------------------------------------------------------!
  !                       band plot information                        !
  !--------------------------------------------------------------------!

  subroutine mbpp_read_bandplot()

    implicit none

    character(len=lenline) :: line
    integer                :: ikpt

    call mbpp_readline(line)
    call io_readval(line, name='nline', value=bandplot_nline)
    call io_readval(line, name='nband', value=bandplot_nband)

    allocate(bandplot_kpt(3,bandplot_nline), &
             bandplot_nstep(bandplot_nline))
    
    call mbpp_readline(line)
    read(line, *) bandplot_start

    do ikpt = 1, bandplot_nline
       call mbpp_readline(line)
       read(line, *) bandplot_kpt(1:3, ikpt), bandplot_nstep(ikpt)
    end do

  end subroutine mbpp_read_bandplot

  !--------------------------------------------------------------------!
  !                 read tight-binding options section                 !
  !--------------------------------------------------------------------!

  subroutine mbpp_read_tboptions()

    implicit none

    character(len=lenline) :: line

    call mbpp_readline(line)
    call io_readval(line, name='AtomsVersion', value=AtomsVersion)
    call io_readval(line, name='BondsVersion', value=BondsVersion)
    call io_readval(line, name='Rcmax',        value=Rcmax)

  end subroutine mbpp_read_tboptions

  !--------------------------------------------------------------------!
  !         parse numeric vector containing literal constants          !
  !--------------------------------------------------------------------!

  function mbpp_parse_coo(str) result(coo)

    implicit none
    
    character(len=*),               intent(in) :: str
    double precision, dimension(3)             :: coo

    character(len=*), parameter :: signum    = '+-.'
    character(len=*), parameter :: operators = '*/'
    character(len=*), parameter :: constants = 'otsnrhlexy'

    double precision, dimension(10) :: value
    
    character         :: chr, op
    double precision  :: buff
    integer           :: i1, i2
    integer           :: icoo, ipos

    value( 1) = 1.D0             ! o
    value( 2) = 3.D0             ! t
    value( 3) = 7.D0             ! s
    value( 4) = 9.D0             ! n
    value( 5) = sqrt(2.D0)       ! r
    value( 6) = sqrt(3.D0)       ! h
    value( 7) = sqrt(5.D0)       ! l
    value( 8) = sqrt(10.D0)      ! e
    value( 9) = sqrt(8.D0/3.D0)  ! x
    value(10) = sqrt(6.D0)       ! y

    coo(1:3) = 1.0d0
    op    = "*"
    icoo  = 1
    i1 = 1
    i2 = 1
    parse : do while(i2 <= len_trim(str))

       chr = str(i2:i2)

       ! not part of a number:
       if ( (.not. io_isfigure(chr)) .and. (scan(signum, chr) == 0) ) then

          ! if anything was read so far, it was a number:
          if ( ((i2-1) >= i1) .and. (len_trim(str(i1:i2-1))>0) ) then
             read(str(i1:i2-1), *) buff
             ! if no operator is in memory, the last value is done:
             if (op == ' ') then
                icoo  = icoo + 1
                coo(icoo) = buff
             ! otherwise apply the operator:
             else if (op == "*") then
                coo(icoo) = coo(icoo)*buff
             else
                coo(icoo) = coo(icoo)/buff
             end if
             ! reset operator:
             op = ' '
          end if

          ! exit, if the rest of the line is commented:
          if (chr == '!') exit parse

          ! if a named constant was found, evaluate it:
          ipos = scan(constants, chr)
          if (ipos > 0) then
             buff = value(ipos)
             ! if no operator is in memory, the last value is done:
             if (op == ' ') then
                icoo  = icoo + 1
                coo(icoo) = buff
             ! otherwise apply the operator:
             else if (op == "*") then
                coo(icoo) = coo(icoo)*buff
             else
                coo(icoo) = coo(icoo)/buff
             end if
             ! reset operator:
             op = ' '
          end if

          ! if an operator was found, remember it:
          ipos = scan(operators, chr)
          if (ipos > 0) then
             op = chr
          end if         

          i1 = i2+1
          
       end if ! not part of a number:
          
       i2 = i2 + 1

    end do parse

    ! if the string was not parsed till the end:
    if ( ((i2-1) >= i1) .and. (len_trim(str(i1:i2-1))>0) ) then
       read(str(i1:i2-1), *) buff
       ! if no operator is in memory, the last value is done:
       if (op == ' ') then
          icoo  = icoo + 1
          coo(icoo) = buff
       ! otherwise apply the operator:
       else if (op == "*") then
          coo(icoo) = coo(icoo)*buff
       else
          coo(icoo) = coo(icoo)/buff
       end if
       ! reset operator:
       op = ' '
    end if

    if (icoo /= 3) then
       write(0,*) 'Error: Invalid string for numeric evaluation in MBPP_eval_coo():'
       write(0,*) trim(str)
       stop
    end if

  end function mbpp_parse_coo


end module mbppio
