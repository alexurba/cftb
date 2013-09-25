module tbparam

  !--------------------------------------------------------------------!
  ! read the tight-binding parameters and tabulate them using          !
  ! bondint.f90                                                        !
  !--------------------------------------------------------------------!
  ! 2010-11-02 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use constants, only: DIRSEP

  use io,        only: io_lower,           &
                       io_readnext,        &
                       io_readval

  use bondint,   only: bi_init,            &
                       bi_final,           &
                       bi_set_onsitelevel, &
                       bi_tabulate

  implicit none
  save

  public  :: tbp_init,            &
             tbp_final

  private :: tbp_parse_atom,      &
             tbp_parse_info,      &
             tbp_parse_matparams, &
             tbp_parse_params,    &
             convert_energy,      &
             convert_length

  integer, parameter :: lenline = 2048

  !--------------------------------------------------------------------!
  ! nl(itype)      : number of angular momentum channels of type i     !
  ! l_of_il(il,it) : il-th angular momentum of atom type it            !
  !--------------------------------------------------------------------!
  ! tbp_nElec(it)  : number of valence electrons of type it            !
  ! tbp_mass(it)   : atomic mass of atom type it                       !
  !--------------------------------------------------------------------!
  ! tbp_*          : * in H,S,O --> .true., if * is parametrized       !
  ! tbp_pot        : .true., if a (pair-) potential is parametrized    !
  !--------------------------------------------------------------------!
  ! tbp_level(il,i): on-site level of il-th angular mom. of type i     !
  !--------------------------------------------------------------------!

  integer,          dimension(:),   allocatable, public  :: nl
  integer,          dimension(:,:), allocatable, public  :: l_of_il

  integer,          dimension(:),   allocatable, public  :: tbp_nElec
  double precision, dimension(:),   allocatable, public  :: tbp_mass

  logical,                                       public  :: tbp_H
  logical,                                       public  :: tbp_S
  logical,                                       public  :: tbp_O
  logical,                                       public  :: tbp_pot

  double precision, dimension(:,:), allocatable, public  :: tbp_level

contains

  subroutine tbp_init(nlmax, ntypes, type_name, param_name, root_dir)
    
    implicit none

    integer,                              intent(in) :: ntypes
    integer,                              intent(in) :: nlmax
    character(len=2),  dimension(ntypes), intent(in) :: type_name
    character(len=*),                     intent(in) :: param_name
    character(len=*),                     intent(in) :: root_dir

    integer, parameter :: u_bond = 20
    integer, parameter :: NPMAX  = 100

    character(len=100) :: fname, dir
    logical            :: fexists

    character(len=2)   :: name1, name2
    double precision   :: rmin, rmax, rmin_pot, rmax_pot
    integer            :: itype1, itype2, it1, it2, i, n
    character(len=lenline) :: line
    integer            :: status
    character(len=30)  :: keyword, funct, eunit, lunit
    integer            :: l, il
    integer            :: l1, l2, m
    double precision   :: r0, r1, efac, lfac
    logical            :: reverse
    integer                            :: np
    double precision, dimension(NPMAX) :: param

    allocate(nl(ntypes), l_of_il(nlmax, ntypes), tbp_level(nlmax, ntypes), &
             tbp_nElec(ntypes), tbp_mass(ntypes) )

    !---------------------------------------------!
    ! read information about the TB parameter set !
    !---------------------------------------------!

    dir = trim(root_dir)
    if (dir(len_trim(dir):len_trim(dir)) /= DIRSEP) then
       dir = trim(dir) // DIRSEP
    end if

    fname = trim(adjustl(dir)) // trim(adjustl(param_name)) // '.info'
    call tbp_parse_info(fname, rmin, rmax, rmin_pot, rmax_pot, tbp_H, &
                        tbp_S, tbp_O, tbp_pot)


    !-----------------------------------------------!
    ! read on-site levels and count angular momenta !
    ! --> nl(itype), l_of_il(il,itype), l_max       !
    !-----------------------------------------------!

    do itype1 = 1, ntypes
       fname = trim(adjustl(dir)) // trim(adjustl(type_name(itype1))) &
               // '.' // trim(adjustl(param_name)) // '.dat'
       call tbp_parse_atom(fname, itype1, ntypes, nlmax, nl, l_of_il, tbp_level)
    end do


    ! initialize bond integrals table:
    call bi_init(ntypes, nlmax, nl, l_of_il, rmin, rmax, tbp_H, tbp_S, &
                 tbp_O, rmin_pot, rmax_pot, tbp_pot)


    !--------------------!
    ! set on-site levels !
    !--------------------!

    do itype1 = 1, ntypes
       do il = 1, nl(itype1)
          l = l_of_il(il, itype1)
          call bi_set_onsitelevel(l, itype1, tbp_level(il, itype1))
       end do
    end do

    !-----------------------------------------------------------!
    ! loop over all parameter files and load the bond integrals !
    !-----------------------------------------------------------!

    do itype1 = 1, ntypes
    do itype2 = itype1, ntypes
       name1 = type_name(itype1)
       name2 = type_name(itype2)
       fname = trim(adjustl(dir)) // trim(name1) // '-' // trim(name2) &
               // '.' // trim(adjustl(param_name)) // '.dat'
       inquire(file=trim(fname), exist=fexists)
       if (fexists) then
          it1 = itype1
          it2 = itype2
       else
          fname = trim(adjustl(dir)) // trim(name2) // '-' // trim(name1) &
                  // '.' // trim(adjustl(param_name)) // '.dat'
          inquire(file=trim(fname), exist=fexists)
          it1 = itype2
          it2 = itype1          
       end if
       if (.not. fexists) then
          write(0,*) "Error: file not found in `tbp_init()': " // trim(fname)
          stop
       end if

       open(u_bond, file=trim(fname), status='old', action='read')
       readloop : do
          call tbp_readline(u_bond, line, status)
          if (status /= 0) exit readloop
          read(line, *) keyword
          select case(trim(io_lower(keyword)))
          ! Hamilton matrix:
          case('hamilton', 'hopping', 'h')
             call io_readval(line, name='funct',   value=funct)
             call io_readval(line, name='EUnit',   value=eunit)
             call io_readval(line, name='LenUnit', value=lunit)
             efac = convert_energy(eunit)
             lfac = convert_length(lunit)
             read(line, *) keyword, n
             do i = 1, n
                call tbp_readline(u_bond, line, status)
                if (tbp_H) then
                   call tbp_parse_matparams(line, NPMAX, l1, l2, &
                                            m, r0, r1, param, np)
                   call bi_tabulate(trim(funct), it1, it2, 'H',  &
                        l1, l2, m, r0, r1, param(1:np), xscale=lfac, yscale=efac)
                end if
             end do
          ! Overlap matrix:
          case('overlap', 'screening', 's')
             call io_readval(line, name='funct',   value=funct)
             call io_readval(line, name='LenUnit', value=lunit)
             lfac = convert_length(lunit)
             read(line, *) keyword, n
             do i = 1, n
                call tbp_readline(u_bond, line, status)
                if (tbp_S) then
                   call tbp_parse_matparams(line, NPMAX, l1, l2, &
                                            m, r0, r1, param, np)
                   call bi_tabulate(trim(funct), it1, it2, 'S',  &
                        l1, l2, m, r0, r1, param(1:np), xscale=lfac)
                end if
             end do
          ! On-site matrix:
          case('onsite', 'o')
             call io_readval(line, name='funct',   value=funct)
             call io_readval(line, name='EUnit',   value=eunit)
             call io_readval(line, name='LenUnit', value=lunit)
             call io_readval(line, name='reverse', value=reverse)
             efac = convert_energy(eunit)
             lfac = convert_length(lunit)
             read(line, *) keyword, n
             do i = 1, n
                call tbp_readline(u_bond, line, status)
                if (tbp_O) then
                   call tbp_parse_matparams(line, NPMAX, l1, l2, &
                                            m, r0, r1, param, np)
                   if (reverse) then
                      ! reverse direction parameters type2 --> type1:
                      call bi_tabulate(trim(funct), it2, it1, 'O',  &
                           l1, l2, m, r0, r1, param(1:np), xscale=lfac, yscale=efac)
                   else
                      ! parameters type1 --> type2:
                      call bi_tabulate(trim(funct), it1, it2, 'O',  &
                           l1, l2, m, r0, r1, param(1:np), xscale=lfac, yscale=efac)
                   end if
                end if
             end do
          ! pair potential matrix:
          case('potential', 'pot')
             call io_readval(line, name='funct',   value=funct)
             call io_readval(line, name='EUnit',   value=eunit)
             call io_readval(line, name='LenUnit', value=lunit)
             efac = convert_energy(eunit)
             lfac = convert_length(lunit)
             call tbp_readline(u_bond, line, status)
             if (tbp_pot) then
                call tbp_parse_params(line, NPMAX, r0, r1, param, np)
                call bi_tabulate(trim(funct), it1, it2, 'pairpot', &
                                 r0, r1, param(1:np), xscale=lfac, yscale=efac)
             end if
          case default
             write(0,*) "Warning: ignoring unknown keyword: ", trim(keyword)
          end select
       end do readloop
       close(u_bond)

    end do
    end do

  end subroutine tbp_init

  !--------------------------------------------------------------------!

  subroutine tbp_final()

    implicit none

    call bi_final()
    if (allocated(nl)) deallocate(nl, l_of_il, tbp_level, tbp_nElec, tbp_mass)

  end subroutine tbp_final

  !--------------------------------------------------------------------!
  !             parser routines for different input files              !
  !--------------------------------------------------------------------!

  subroutine tbp_parse_atom(fname, itype, ntypes, nlmax, nl, l_of_il, level)

    implicit none 

    character(len=*),                           intent(in)    :: fname
    integer,                                    intent(in)    :: ntypes, itype
    integer,                                    intent(in)    :: nlmax
    integer,          dimension(ntypes),        intent(inout) :: nl
    integer,          dimension(nlmax, ntypes), intent(inout) :: l_of_il
    double precision, dimension(nlmax, ntypes), intent(inout) :: level

    integer, parameter :: u_atom = 20

    character(len=lenline) :: line
    integer            :: status
    integer            :: l, il
    character(len=30)  :: keyword, eunit, munit
    double precision   :: fac, val
    logical            :: fexists

    ! does the file exists ?
    inquire(file=trim(fname), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: file not found in `tbp_parse_atom()': " // trim(fname)
       stop
    end if

    level(:,itype)   = 0.0d0
    tbp_mass(itype)  = 1.0d0
    tbp_nElec(itype) = 0

    open(u_atom, file=trim(fname), status='old', action='read')
    readloop : do

       call tbp_readline(u_atom, line, status)
       if (status /= 0) exit readloop

       read(line, *) keyword
       select case(trim(io_lower(keyword)))
       case('electrons')
          read(line, *) keyword, tbp_nElec(itype)
       case('levels')
          call io_readval(line, name="EUnit", value=eunit)
          fac = convert_energy(eunit)
          read(line, *) keyword, nl(itype)
          if (nl(itype) > nlmax) then
             write(0,*) 'Error: number of different angular momenta too high!'
             close(u_atom)
             stop
          end if
          do il = 1, nl(itype)
             read(u_atom, *) l, val
             l_of_il(il, itype) = l
             level(il,itype)    = fac*val
          end do
       case('mass')
          call io_readval(line, name="MUnit", value=munit)
          fac = convert_mass(munit)
          read(line, *) keyword, tbp_mass(itype)
          tbp_mass(itype) = fac*tbp_mass(itype)
       case('name') 
          continue
       case default
          write(0,*) 'Warning: ignoring unknown keyword: ' // trim(keyword)
          write(0,*) '         in file: ' // trim(fname)
       end select

    end do readloop
    close(u_atom)

  end subroutine tbp_parse_atom

  !--------------------------------------------------------------------!

  subroutine tbp_parse_info(fname, rmin, rmax, rmin_pot, rmax_pot, &
                            Hmat, Smat, Omat, pot)

    implicit none

    character(len=*), intent(in)  :: fname
    double precision, intent(out) :: rmin, rmax
    double precision, intent(out) :: rmin_pot, rmax_pot
    logical,          intent(out) :: Hmat, Smat, Omat, pot

    integer, parameter :: u_info = 20

    character(len=lenline) :: line
    integer            :: status, n, i
    character(len=30)  :: keyword, unit, str
    double precision   :: fac
    logical            :: fexists

    ! file exists ?
    inquire(file=trim(fname), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: file not found in `tbp_parse_info()': " // trim(fname)
       stop
    end if

    ! default values:
    rmin     = 0.0d0
    rmax     = 0.0d0
    rmin_pot = 0.0d0
    rmax_pot = 0.0d0
    Hmat = .false.
    Smat = .false.
    Omat = .false.
    pot  = .false.

    open(u_info, file=trim(fname), status='old', action='read')

    readloop : do

       call tbp_readline(u_info, line, status)
       if (status /= 0) exit readloop

       read(line, *) keyword
       select case(trim(io_lower(keyword)))
       case('range')
          call io_readval(line, name='rmin',    value=rmin)
          call io_readval(line, name='rmax',    value=rmax)
          call io_readval(line, name='LenUnit', value=unit)
          fac  = convert_length(unit)
          rmin = fac*rmin
          rmax = fac*rmax
       case('potrange', 'vrange')
          call io_readval(line, name='rmin',    value=rmin_pot)
          call io_readval(line, name='rmax',    value=rmax_pot)
          call io_readval(line, name='LenUnit', value=unit)
          fac  = convert_length(unit)
          rmin_pot = fac*rmin_pot
          rmax_pot = fac*rmax_pot
       case('hamilton', 'hopping', 'h')
          read(line,*) keyword, str
          select case(trim(io_lower(str)))
          case('true', 't', '1')
             Hmat = .true.
          case default
             Hmat = .false.
          end select
       case('overlap', 'screening', 's')
          read(line,*) keyword, str
          select case(trim(io_lower(str)))
          case('true', 't', '1')
             Smat = .true.
          case default
             Smat = .false.
          end select
       case('onsite', 'o')
          read(line,*) keyword, str
          select case(trim(io_lower(str)))
          case('true', 't', '1')
             Omat = .true.
          case default
             Omat = .false.
          end select
       case('potential', 'pot', 'pairpot')
          read(line,*) keyword, str
          select case(trim(io_lower(str)))
          case('true', 't', '1')
             pot = .true.
          case default
             pot = .false.
          end select
       case('files')
          read(line,*) keyword, n
          do i = 1, n
             read(u_info,*)
          end do
       case default
          write(0,*) "Warning: ignoring unknown keyword: ", trim(keyword)
       end select

    end do readloop

    close(u_info)

  end subroutine tbp_parse_info

  !--------------------------------------------------------------------!

  subroutine tbp_readline(u_in, line, iostat)

    implicit none

    integer,            intent(in)  :: u_in
    character(len=*),   intent(out) :: line
    integer,            intent(out) :: iostat

    do
       read(u_in, '(A)', iostat=iostat) line
       if (iostat /= 0) exit
       line = trim(adjustl(line))
       if (len_trim(line) == 0) cycle
       if (scan(line(1:1), '#!;/') /= 0) cycle
       exit
    end do

  end subroutine tbp_readline

  !--------------------------------------------------------------------!

  subroutine tbp_parse_matparams(line, NPMAX, l1, l2, m, r0, r1, param, np)

    implicit none

    character(len=*),                   intent(in)  :: line
    integer,                            intent(in)  :: NPMAX
    integer,                            intent(out) :: l1, l2, m
    double precision,                   intent(out) :: r0, r1
    double precision, dimension(NPMAX), intent(out) :: param
    integer,                            intent(out) :: np

    integer :: pos
    
    pos = 1
    call io_readnext(line, pos, l1)
    call io_readnext(line, pos, l2)
    call io_readnext(line, pos, m)

    call tbp_parse_params(line(pos:len(line)), NPMAX, r0, r1, param, np)

  end subroutine tbp_parse_matparams

  !--------------------------------------------------------------------!

  subroutine tbp_parse_params(line, NPMAX, r0, r1, param, np)

    implicit none

    character(len=*),                   intent(in)  :: line
    integer,                            intent(in)  :: NPMAX
    double precision,                   intent(out) :: r0, r1
    double precision, dimension(NPMAX), intent(out) :: param
    integer,                            intent(out) :: np

    integer :: pos, ip
    
    pos = 1
    call io_readnext(line, pos, r0)
    call io_readnext(line, pos, r1)

    np = NPMAX
    readnext : do ip = 1, NPMAX
       call io_readnext(line, pos, param(ip))
       if (pos == 0) then
          np = ip - 1
          exit readnext
       end if
    end do readnext

  end subroutine tbp_parse_params

  !--------------------------------------------------------------------!
  !                      unit conversion routines                      !
  !--------------------------------------------------------------------!

  function convert_energy(unit) result(fac)

    implicit none

    character(len=*), intent(in) :: unit
    double precision             :: fac

    fac = 1.0d0

    ! the default energy unit is Hartree (a.u.)
    select case(trim(io_lower(unit)))
    case('rydberg', 'ry')
       fac = 0.5d0
    case('evolt','ev')
       fac = 0.0367493088677d0
    end select

  end function convert_energy

  !--------------------------------------------------------------------!

  function convert_length(unit) result(fac)

    implicit none

    character(len=*), intent(in) :: unit
    double precision             :: fac

    fac = 1.0d0

    ! the default length unit is Bohr
    select case(trim(io_lower(unit)))
    case('a', 'angstroem', 'angstrom')
       fac = 1.88972598501d0
    end select

  end function convert_length

  !--------------------------------------------------------------------!

  function convert_mass(unit) result(fac)

    implicit none

    character(len=*), intent(in) :: unit
    double precision             :: fac

    fac = 1.0d0

    ! the default mass unit is electron mass `m_e'
    select case(trim(io_lower(unit)))
    case('amu', 'u')
       fac = 1822.82232834d0
    case('kg')
    end select

  end function convert_mass


end module tbparam
