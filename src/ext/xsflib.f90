module xsflib

  !--------------------------------------------------------------------!
  ! A library with I/O routines to access the XCrysDen Structure File  !
  ! format (XSF).                                                      !
  !--------------------------------------------------------------------!
  ! 2011-02-11 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use io, only: io_adjustl,  &
                io_isfigure, &
                io_isletter, &
                io_readnext, &
                io_readval,  &
                io_lower,    &
                io_upper      

  implicit none

  public :: xsf_init,               &
            xsf_final,              &
            xsf_get_alat,           &
            xsf_get_avec,           &
            xsf_get_celvol,         &
            xsf_get_coorat,         &
            xsf_get_forlat,         &
            xsf_get_nameat,         &
            xsf_get_natoms,         &
            xsf_get_nkxyz,          &
            xsf_get_ntype,          &
            xsf_get_shift,          &
            xsf_get_niter,          &
            xsf_get_ifmet,          &
            xsf_get_smearwidth

  !--------------------------------------------------------------------!
  ! lenline : max. length of lines (number of char's) in input file    !
  ! lentype : max. length of type names (number of char's)             !
  !--------------------------------------------------------------------!
  ! u_in    : unit number of the input file                            !
  !--------------------------------------------------------------------!

  integer, parameter, private :: lenline   = 132
  integer, parameter, private :: lentype   = 2

  integer, parameter, private :: u_in    = 20

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
  ! forlat(i,iat,ityp) : lattice forces in the same format             !
  !                                                                    !
  !------------------------- k-point sampling -------------------------!
  ! nband              : number of bands to be calculated              !
  ! nkxyz(1:3)         : numbers of k-points per dimension             !
  ! shift(1:3)         : shift vector of the k-point mesh              !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! niter              : max. number of SCF iterations                 !
  ! ifmet              : .true., if the system is metallic             !
  ! smearwidth         : width [eV] for smearing of the Fermi edge     !
  !                                                                    !
  !------------------------------ output ------------------------------!
  ! output             : .true., if file contains output information   !
  ! forces             : .true., if forces are present                 !
  ! E_tot              : total energy                                  !
  ! E_coh              : cohesive energy                               !
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

  character(len=lenline),                                public  :: kpmeth
  integer,                                               public  :: nband 
  integer,                dimension(3),                  public  :: nkxyz 
  double precision,       dimension(3),                  public  :: shift 

  integer,                                               public  :: niter     
  logical,                                               public  :: ifmet     
  double precision,                                      public  :: smearwidth

  logical,                                               public  :: output
  logical,                                               public  :: forces
  double precision,                                      public  :: E_tot
  double precision,                                      public  :: E_coh

  logical,                                               public  :: tbenergies
  double precision,                                      public  :: E_tbband
  double precision,                                      public  :: E_tbbond
  double precision,                                      public  :: E_pairpot


  !--------------------------------------------------------------------!
  !                              PRIVATE                               !
  !                                                                    !
  ! iline   : number of the current line of the input file             !
  ! line    : input buffer for single line from the input file         !
  !--------------------------------------------------------------------!

  integer,                                               private :: iline
  character(len=lenline),                                private :: line
  
contains

  subroutine xsf_init(file)

    implicit none

    character(len=*), intent(in) :: file

    logical                      :: fexists
    integer                      :: status
    character(len=lentype)       :: el

    character(len=lenline)       :: str
    integer                      :: i, nat, k, ityp

    inquire(file=trim(file), exist=fexists)
    if (.not. fexists) then
       write(0,*) 'Error: File not found: ', trim(file)
       return
    end if

    ! initial values:
    output     = .false.
    forces     = .false.
    E_tot      = 0.0d0
    E_coh      = 0.0d0
    niter      = 1
    nkxyz(1:3) = (/ 1, 1, 1 /)
    shift(1:3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
    ifmet      = .false.
    smearwidth = 0.0d0

    tbenergies = .false.
    E_tbband   = 0.0d0
    E_tbbond   = 0.0d0
    E_pairpot  = 0.0d0

    open(u_in, file=trim(file), status='old', action='read')

    ! title::
    read(u_in, '(A)') line
    if (len_trim(line(2:lenline)) > 0) then
       read(line(2:lenline), '(A)') title
       title = trim(adjustl(title))
    end if

    iline = 2
    read_input : do

       read(u_in, '(A)', iostat=status) line
       if (status /= 0) exit read_input
       iline = iline + 1

       if (len_trim(line) == 0) cycle read_input

       ! total energy:
       if (index(line, 'total energy') > 0) then
          call io_readval(line, name='total energy', value=E_tot)
          output = .true.
          cycle read_input
       end if

       ! cohesive energy:
       if (index(line, 'cohesive energy') > 0) then
          call io_readval(line, name='cohesive energy', value=E_coh)
          output = .true.
          cycle read_input
       end if

       ! SCF iterations:
       if (index(line, 'SCF iterations') > 0) then
          call io_readval(line, name='SCF iterations', value=niter)
          cycle read_input
       end if

       ! number of k-points:
       if (index(line, 'number of k-pts.') > 0) then
          call io_readval(line, name='number of k-pts.', value=nkxyz, n=3)
          cycle read_input
       end if

       ! shift of k-points grid:
       if (index(line, 'k-point shift') > 0) then
          call io_readval(line, name='k-point shift', value=shift, n=3)
          cycle read_input
       end if

       ! insulator/metal:
       if (index(line, 'metallic') > 0) then
          if (scan(line(10:lenline), 'T') > 0) then
             ifmet = .true.
          else
             ifmet = .false.
          end if
          cycle read_input
       end if

       ! smearing width:
       if (index(line, 'smearing width') > 0) then
          call io_readval(line, name='smearing width', value=smearwidth)
          cycle read_input
       end if

       !--------------------- TB energies        ----------------------!

       ! TB-Band cohesive energy:
       if (index(line, 'TB-Band coh. energy') > 0) then
          call io_readval(line, name='TB-Band coh. energy', value=E_tbband)
          tbenergies = .true.
          cycle read_input
       end if

       ! TB-Bond cohesive energy:
       if (index(line, 'TB-Bond coh. energy') > 0) then
          call io_readval(line, name='TB-Bond coh. energy', value=E_tbbond)
          tbenergies = .true.
          cycle read_input
       end if

       ! pair potential energy:
       if (index(line, 'Pair potential energy') > 0) then
          call io_readval(line, name='Pair potential energy', value=E_pairpot)
          tbenergies = .true.
          cycle read_input
       end if
       

       ! skip all remaining comments:
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle read_input

       ! lattice vectors:
       if (index(line, 'PRIMVEC') > 0) then
          read(u_in, *) avec(1:3, 1)
          read(u_in, *) avec(1:3, 2)
          read(u_in, *) avec(1:3, 3)
          iline = iline + 3
          cycle read_input
       end if

       ! lattice coordinates and forces:
       primcoord : if (index(line, 'PRIMCOORD') > 0) then

          read(u_in, *) nat

          ! count atom types and atoms:
          ntype   = 0
          str     = ''
          natomax = 0
          k       = 0
          do i = 1, nat
             read(u_in, '(A)') line
             read(line, *) el
             if (trim(el) /= trim(str)) then
                str     = el
                ntype   = ntype + 1
                natomax = max(natomax, k)
                k       = 1
             end if
             k = k + 1
          end do
          natomax = max(natomax, k)
          ! are forces present?
          k = 0
          i = 1
          do while (i > 0)
             k = k + 1
             call io_readnext(line, i, str)
          end do
          if (k >= 7) forces = .true.
          ! rewind to previous position:
          rewind(u_in)
          do i = 1, iline
             read(u_in, *)
          end do

          ! allocate memory:
          if (forces) then
             allocate(coorat(3,natomax,ntype), natoms(ntype), &
                      nameat(ntype), forlat(3,natomax,ntype))
          else
             allocate(forlat(3,natomax,ntype), natoms(ntype), &
                      nameat(ntype))
          end if

          ! read cartesian coordinates and forces:
          el   = ''
          ityp = 0
          do i = 1, nat
             read(u_in, '(A)') line
             iline = iline + 1
             read(line, *) str
             if (trim(str) /= trim(el)) then
                ityp         = ityp + 1
                natoms(ityp) = 0
                nameat(ityp) = trim(str)
             end if
             natoms(ityp) = natoms(ityp) + 1
             if (forces) then
                read(line, *) el, coorat(1:3,natoms(ityp),ityp), &
                                  forlat(1:3,natoms(ityp),ityp)
             else
                read(line, *) el, coorat(1:3,natoms(ityp),ityp)
             end if
          end do
          cycle read_input
       end if primcoord

    end do read_input

    close(u_in)

    ! unit cell volume:
    celvol = avec(1,1)*avec(2,2)*avec(3,3) &
           + avec(2,1)*avec(3,2)*avec(1,3) &
           + avec(3,1)*avec(1,2)*avec(2,3) &
           - avec(3,1)*avec(2,2)*avec(1,3) &
           - avec(1,1)*avec(3,2)*avec(2,3) &
           - avec(2,1)*avec(1,2)*avec(3,3)
   
  end subroutine xsf_init

  !--------------------------------------------------------------------!

  subroutine xsf_final()

    implicit none

    if (allocated(coorat)) then
       deallocate(coorat, nameat, natoms)
    end if
    if (allocated(forlat)) deallocate(forlat)
    
  end subroutine xsf_final

  !--------------------------------------------------------------------!
  !                data access interfaces / properties                 !
  !--------------------------------------------------------------------!

  function xsf_get_title() result(out_title)
    implicit none
    character(len=lenline) :: out_title
    out_title = title
  end function xsf_get_title

  !--------------------------------------------------------------------!

  function xsf_get_ntype() result(out_ntype)
    implicit none
    integer :: out_ntype
    out_ntype = ntype
  end function xsf_get_ntype

  !--------------------------------------------------------------------!

  function xsf_get_natomax() result(out_natomax)
    implicit none
    integer :: out_natomax
    out_natomax = natomax
  end function xsf_get_natomax

  !--------------------------------------------------------------------!

  function xsf_get_alat() result(out_alat)
    implicit none
    double precision :: out_alat
    out_alat = alat
  end function xsf_get_alat

  !--------------------------------------------------------------------!

  function xsf_get_celvol() result(out_celvol)
    implicit none
    double precision :: out_celvol
    out_celvol = celvol
  end function xsf_get_celvol

  !--------------------------------------------------------------------!

  function xsf_get_avec(i, ivec) result(out_avec)
    implicit none
    integer, intent(in) :: i, ivec
    double precision    :: out_avec
    out_avec = avec(i,ivec)
  end function xsf_get_avec

  !--------------------------------------------------------------------!

  function xsf_get_natoms(itype) result(out_natoms)
    implicit none
    integer, intent(in) :: itype
    integer             :: out_natoms
    if (allocated(natoms)) then
       out_natoms = natoms(itype)
    else
       out_natoms = 0
    end if
  end function xsf_get_natoms

  !--------------------------------------------------------------------!

  function xsf_get_nameat(itype) result(out_nameat)
    implicit none
    integer,    intent(in) :: itype
    character(len=lentype) :: out_nameat
    if (allocated(nameat)) then
       out_nameat = nameat(itype)
    else
       out_nameat = ' '
    end if
  end function xsf_get_nameat

  !--------------------------------------------------------------------!

  function xsf_get_coorat(i, iatom, itype) result(out_coorat)
    implicit none
    integer, intent(in) :: i, iatom, itype
    double precision    :: out_coorat
    if (allocated(coorat)) then
       out_coorat = coorat(i, iatom, itype)
    else
       out_coorat = 0.0d0
    end if
  end function xsf_get_coorat

  !--------------------------------------------------------------------!

  function xsf_get_forlat(i, iatom, itype) result(out_forlat)
    implicit none
    integer, intent(in) :: i, iatom, itype
    double precision    :: out_forlat
    if (allocated(forlat)) then
       out_forlat = forlat(i, iatom, itype)
    else
       out_forlat = 0.0d0
    end if
  end function xsf_get_forlat

  !--------------------------------------------------------------------!

  function xsf_get_nband() result(out_nband)
    implicit none
    integer :: out_nband
    out_nband = nband
  end function xsf_get_nband

  !--------------------------------------------------------------------!

  function xsf_get_nkxyz() result(out_nkxyz)
    implicit none
    integer, dimension(3) :: out_nkxyz
    out_nkxyz = nkxyz
  end function xsf_get_nkxyz

  !--------------------------------------------------------------------!

  function xsf_get_shift() result(out_shift)
    implicit none
    double precision, dimension(3) :: out_shift
    out_shift(1:3) = shift
  end function xsf_get_shift

  !--------------------------------------------------------------------!

  function xsf_get_ifmet() result(out_ifmet)
    implicit none
    logical :: out_ifmet
    out_ifmet = ifmet
  end function xsf_get_ifmet

  !--------------------------------------------------------------------!

  function xsf_get_niter() result(out_niter)
    implicit none
    integer :: out_niter
    out_niter = niter
  end function xsf_get_niter

  !--------------------------------------------------------------------!

  function xsf_get_smearwidth() result(out_smearwidth)
    implicit none
    double precision :: out_smearwidth
    out_smearwidth = smearwidth
  end function xsf_get_smearwidth

end module xsflib
