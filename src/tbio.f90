module tbio

  !--------------------------------------------------------------------!
  !           common I/O routines of the tb package programs           !
  !--------------------------------------------------------------------!
  ! 2010-11-17 Alexander Urban (AU)                                    !
  ! 2011-02-20 AU --- Introducing arglib as alternative interface to   !
  !                   command line arguments.  The old interface is    !
  !                   depricated and shall not be used any longer.     !
  ! 2011-03-28 AU --- All TB package tools have been moved to the new  !
  !                   interface and the old one was removed.           !
  !--------------------------------------------------------------------!

  use arglib,    only: args_init,     &
                       args_final,    &
                       args_env,      &
                       args_switch,   &
                       args_get,      &
                       nargs

  use constants, only: CURDIR,        &
                       DIRSEP

  use io,        only: io_adjustl,    &
                       io_replace,    &
                       io_trim

  implicit none

  public  :: tbio_init, tbio_final,  &
             tbio_default,           &
             tbio_filecheck

  private :: tbio_default_i,         &
             tbio_default_l,         &
             tbio_default_c,         &
             tbio_default_d,         &
             set_default,            &
             set_default_i,          &
             set_default_l,          &
             set_default_c,          &
             set_default_d,          &
             def_index

  integer,     parameter, private :: lenname = 64
  integer,     parameter, private :: lendes  = 64
  integer,     parameter, private :: lenswi  = 256
  integer,     parameter, private :: lenval  = 256

  !--------------------------------------------------------------------!

  character(len=lenname), private :: tbio_name
  character(len=lenval),  private :: tbio_descr

  integer,                public  :: wdescr
  integer,                public  :: wswitch

  !------------------------- input data type --------------------------!

  type, private :: tbio_data
     character(len=lenname) :: name     = ''
     character(len=lendes)  :: descr    = ''
     character              :: type     = ''
     character(len=lenswi)  :: switches = ''
     character(len=lenval)  :: env      = ''
     character(len=lenval)  :: help     = ''
     integer                :: value_i  = 0
     logical                :: value_l  = .false.
     character(len=lenval)  :: value_c  = ''
     double precision       :: value_d  = 0.0d0
  end type tbio_data
  type(tbio_data),        dimension(:), allocatable, private :: defaults
  integer,                                           private :: nDefaults = 0
  integer,                                           private :: idef      = 0
  character(len=lenname), dimension(:), allocatable, private :: option
  integer,                                           private :: nOptions  = 0

  !--------------------------------------------------------------------!
  !                       tbio_args() interface                        !
  !                                                                    !
  !    look up specified switches, return and print out the result     !
  !--------------------------------------------------------------------!

  interface tbio_args
     module procedure tbio_args_i, tbio_args_l, tbio_args_c, tbio_args_d
  end interface

  !--------------------------------------------------------------------!
  !                      tbio_default() interface                      !
  !                                                                    !
  ! get the default value for a given name via:                        !
  !                                                                    !
  ! call tbio_default(name, value)                                     !
  ! --> implementations available: integer, logical, character, double !
  !--------------------------------------------------------------------!

  interface tbio_default
     module procedure tbio_default_i, tbio_default_l, tbio_default_c, &
                      tbio_default_d
  end interface

  !--------------------------------------------------------------------!
  !                      set_default() interface                       !
  !                                                                    !
  ! add entry to the `defaults' list via: call set_default(name, value)!
  ! --> implementations available: integer, logical, character, double !
  !--------------------------------------------------------------------!

  interface set_default
     module procedure set_default_i, set_default_l, set_default_c, &
                      set_default_d
  end interface




contains !=============================================================!




  subroutine tbio_init(name, descr, args, extra)

    implicit none

    !------------------------------------------------------------------!
    ! name       : program name (e.g. exe.x)                           !
    ! descr      : short program description                           !
    ! args(:)    : list of allowed command line options                !
    !------------------------------------------------------------------!

    character(len=*),               intent(in) :: name
    character(len=*),               intent(in) :: descr
    character(len=*), dimension(:), intent(in) :: args
    character(len=*), optional,     intent(in) :: extra

    character(len=512) :: switches
    logical            :: help
    integer            :: i, iopt, stat

    ! save program information:
    tbio_name  = trim(name)
    tbio_descr = trim(descr)
    nOptions   = size(args)
    allocate(option(nOptions))
    option(:)  = args(:)

    !------------------------------------------------------------------!

    nDefaults = 22
    idef      = 0
    allocate(defaults(nDefaults))

    !---------------------------- general -----------------------------!
    ! format     : type of the input file (mbpp, ?)                    !
    ! in         : input file name                                     !
    ! out        : output file name                                    !
    ! dir        : directory path with TB parameters                   !
    ! set        : name of the TB parameter set                        !
    ! time       : turn on timing, if available                        !
    ! niter      : number of iterations                                !
    ! tol        : tolerance, treshhold for convergence                !
    ! files      : .true., if data files are listed                    !
    ! plot       : name of output file with data for plotting          !
    ! eig        : name of file with eigenvalues                       !
    ! restart    : name of file with restart information               !
    !------------------------------------------------------------------!

    call set_default(name     = 'dir',                                   &
                     descr    = 'Parameters directory',                  &
                     value    = CURDIR // DIRSEP // 'tbparam' // DIRSEP, &
                     switches = '--dir:-d',                              &
                     help     = 'set the parameters directory path',     &
                     env      = 'TB_PARAMDIR')
    call set_default(name     = 'format',                                &
                     descr    = 'Input file format',                     &
                     value    = 'mbpp',                                  &
                     switches = '--format:-f',                           &
                     help     = 'set the input file format')
    call set_default(name     = 'in',                                    &
                     descr    = 'Input file name',                       &
                     value    = 'INP',                                   &
                     switches = '--in:-i',                               &
                     help     = 'set input file')
    call set_default(name     = 'out',                                   & 
                     descr    = 'Output file name',                      &
                     value    = 'output.dat',                            &
                     switches = '--out:-o',                              &
                     help     = 'set output file name')
    call set_default(name     = 'param',                                 &
                     descr    = 'TB parameter set',                      &
                     value    = 'proj01',                                &
                     switches = '--param:-p',                            &
                     help     = 'choose a TB parameter set')
    call set_default(name     = 'time',                                  &
                     descr    = "Write timing info to `timing.log'.",    &
                     value    = .false.,                                 &
                     switches = '--time',                                &
                     help     = 'enable timing (output to timing.log)')
    call set_default(name     = 'niter',                                 &
                     descr    = 'Max. number of iterations',             &
                     value    = 10000,                                   &
                     switches = '--niter',                               &
                     help     = 'set number of iterations')
    call set_default(name     = 'tol',                                   &
                     descr    = 'Optimization threshold',                &
                     value    = 1.0d-6,                                  &
                     switches = '--tol',                                 &
                     help     = 'set optimization threshold')
    call set_default(name     = 'files',                                 &
                     descr    = 'Data files given.',                     &
                     value    = .false.,                                 &
                     switches = '--files:--train',                       &
                     help     = 'list of data files')
    call set_default(name     = 'plot',                                  &
                     descr    = 'Output file for plot',                  &
                     value    = 'output.plt',                            &
                     switches = '--plot',                                &
                     help     = 'specify file for plot data')
    call set_default(name     = 'eig',                                   &
                     descr    = 'File with eigenvalues',                 &
                     value    = 'output.dat.eig',                        &
                     switches = '--eig',                                 &
                     help     = 'specify file with eigenvalues')
    call set_default(name     = 'restart',                               &
                     descr    = 'File with restart information',         &
                     value    = '',                                      &
                     switches = '--restart',                             &
                     help     = 'specify file with restart information')

    !---------------------------- notb.f90 ----------------------------!
    ! occup      : file with occupation numbers                        !
    !------------------------------------------------------------------!

    call set_default(name     = 'occup',                              &
                     descr    = 'File with occupation numbers',       &
                     value    = '',                                   &
                     switches = '--occup',                            &
                     help     = 'specify file with occupation numbers')

    !--------------------------- bands.f90 ----------------------------!
    ! gnuplot    : .true., if a gnuplot script shall be created        !
    ! plotbz     : .true., if Brillouin zone shall be plotted          !
    ! nountangle : .true., if bands shall not be untangled             !
    ! ezero      : energy zero point                                   !
    !------------------------------------------------------------------!

    call set_default(name     = 'gnuplot',                         &
                     descr    = 'Gnuplot file will be created.',   &
                     value    = .false.,                           &
                     switches = '--gnuplot',                       &
                     help     = 'enable output of a gnuplot script')
    call set_default(name     = 'plotbz',                          &
                     descr    = 'Brillouin zone plot.',            &
                     value    = .false.,                           &
                     switches = '--plotbz',                        &
                     help     = 'enable Brillouin zone plot')
    call set_default(name     = 'nountangle',                      &
                     descr    = 'Untangling of bands disabled.',   &
                     value    = .false.,                           &
                     switches = '--nountangle',                    &
                     help     = 'disable untangling of bands')
    call set_default(name     = 'ezero',                           &
                     descr    = 'Energy zero [eV]',                &
                     value    = 0.0d0,                             &
                     switches = '--ezero:-z',                      &
                     help     = 'set sero point of energy [eV]')

    !---------------------------- fit.f90 -----------------------------!
    ! funct      : function type to be used                            !
    ! opt        : optimization method                                 !
    ! data       : data file name                                      !
    ! deriv      : file with derivatives                               !
    !------------------------------------------------------------------!

    call set_default(name     = 'funct',                      &
                     descr    = 'Function type',              &
                     value    = 'polyexp',                    &
                     switches = '--funct',                    &
                     help     = 'set type of fitting function')
    call set_default(name     = 'opt',                        &
                     descr    = 'Optimization method',        &
                     value    = 'simplex',                    &
                     switches = '--opt',                      &
                     help     = 'choose optimization method')
    call set_default(name     = 'data',                       &
                     descr    = 'Data file',                  &
                     value    = 'DATA',                       &
                     switches = '--data',                     &
                     help     = 'specify data file')
    call set_default(name     = 'deriv',                      &
                     descr    = 'File with derivatives',      &
                     value    = '',                           &
                     switches = '--deriv',                    &
                     help     = 'specify file with derivatives')

    !--------------------------- potfit.f90 ---------------------------!
    ! nnodes     : number of abscissas for the spline interpolation    !
    !------------------------------------------------------------------!

    call set_default(name     = 'nnodes',                      &
                     descr    = 'Number of nodes/abscissas',   &
                     value    = 5,                             &
                     switches = '--nnodes',                    &
                     help     = 'specify number of nodes/abscissas')

    if (idef /= nDefaults) then
       write(0,*) "Error: inconsistent number of defaults in `tbio'."
       call tbio_final()
       stop
    end if

    !------------------------------------------------------------------!

    switches = '--help:-h'
    wswitch  = len_trim(switches)
    wdescr   = 0
    do iopt = 1, nOptions
       i = def_index(option(iopt))
       switches = trim(switches) // ':' // trim(defaults(i)%switches)
       wswitch  = max(wswitch,  len_trim(defaults(i)%switches))
       if (defaults(i)%type /= 'l') then
          wdescr = max(wdescr, len_trim(defaults(i)%descr))
       end if
    end do
    switches = trim(switches) // ':' // trim(extra)

    call args_init(check=trim(switches), status=stat)
    call args_switch('--help:-h', value=help)
    if (stat /= 0 .or. help) then
       call tbio_print_usage()
       call tbio_final()
       stop
    end if

    !------------------------------------------------------------------!

    write(*,*)

  end subroutine tbio_init

  !--------------------------------------------------------------------!

  subroutine tbio_final()

    implicit none

    write(*,*)
    call args_final()
    if (allocated(defaults)) deallocate(defaults)
    if (allocated(option)) deallocate(option)

  end subroutine tbio_final

  !--------------------------------------------------------------------!
  !                    print customized help screen                    !
  !--------------------------------------------------------------------!

  subroutine tbio_print_usage()

    implicit none

    integer               :: iopt, i
    character(len=lenval) :: def

    write(*,*)
    write(*,*) 'Usage: ', trim(tbio_name), ' [options]'
    write(*,'(1x,70("-"))')
    write(*,*) trim(tbio_descr)
    write(*,'(1x,70("-"))')
    write(*,*)
    write(*,*) 'Options:'
    write(*,*)
    write(*,*) io_trim('--help -h',wswitch), &
               '  show this usage information screen'
    do iopt = 1, nOptions
       i = def_index(option(iopt))
       select case(defaults(i)%type)
       case('i')
          def = trim(io_adjustl(defaults(i)%value_i))
       case('l')
          def = 'off'
       case('c')
          def = trim(adjustl(defaults(i)%value_c))
       case('d')
          def = trim(io_adjustl(defaults(i)%value_d,8))
       end select
       write(*,*) io_trim(io_replace(defaults(i)%switches,':',' '),wswitch), &
                  '  ', trim(defaults(i)%help), ' (default: ', trim(def), ')'
       if (len_trim(defaults(i)%env) > 0) then
          write(*,*) io_trim(' ',wswitch), '  ',                     &
                     'Alternatively set the environment variable: ', &
                     trim(defaults(i)%env)
       end if
    end do
    write(*,'(1x,70("-"))')
    write(*,*)

  end subroutine tbio_print_usage

  !--------------------------------------------------------------------!
  !                       check if a file exists                       !
  !                                                                    !
  ! The subroutine stops the program if the queried file does not exist!
  !--------------------------------------------------------------------!

  subroutine tbio_filecheck(fname)

    implicit none

    character(len=*), intent(in) :: fname
    logical                      :: fexists

    inquire(file=trim(fname), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: File not found: ", trim(fname)
       call tbio_final()
       stop
    end if

  end subroutine tbio_filecheck





  !============================= private ==============================!



  !--------------------------------------------------------------------!
  !            search the defaults list for a certain name             !
  !--------------------------------------------------------------------!

  function def_index(name) result(idx)

    implicit none

    character(len=*), intent(in) :: name

    integer :: idx
    integer :: i

    idx = 0
    do i = 1, nDefaults
       if (trim(defaults(i)%name) == trim(adjustl(name))) then
          idx = i
          exit
       end if
    end do

    if (idx == 0) then
       write(0,*) "Error: unknown data name in `tbio': ", trim(name)
       call tbio_final()
       stop
    end if

  end function def_index

  !--------------------------------------------------------------------!
  !            Implementation of the tbio_args() interface             !
  !--------------------------------------------------------------------!

  subroutine tbio_args_i(name, value, default, silent)

    implicit none

    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: value
    integer, optional, intent(in)  :: default
    logical, optional, intent(in)  :: silent

    integer :: envval
    integer :: i, stat

    i = def_index(name)

    if (defaults(i)%type /= 'i') then
       write(0,*) "Error: incompatible types in `tbio'."
       call tbio_final()
       stop
    end if

    stat = 1
    if (len_trim(defaults(i)%env) > 0) then
       call args_env(defaults(i)%env, value=envval, stat=stat)
    end if

    if (stat == 0) then
       value = envval
    else if (present(default)) then
       value = default
    else
       value = defaults(i)%value_i
    end if

    call args_switch(trim(defaults(i)%switches), value=value)

    if ((.not. present(silent)) .or. (.not. silent)) then
       write(*,*) io_trim(defaults(i)%descr,wdescr), ' : ', trim(io_adjustl(value))
    end if

  end subroutine tbio_args_i

  !--------------------------------------------------------------------!

  subroutine tbio_args_l(name, value, default, silent)

    implicit none

    character(len=*),  intent(in)  :: name
    logical,           intent(out) :: value
    logical, optional, intent(in)  :: default
    logical, optional, intent(in)  :: silent
    
    logical :: envval
    integer :: i

    i = def_index(name)

    if (defaults(i)%type /= 'l') then
       write(0,*) "Error: incompatible types in `tbio'."
       call tbio_final()
       stop
    end if

    envval = .false.
    if (len_trim(defaults(i)%env) > 0) then
       call args_env(defaults(i)%env, value=envval)
    end if

    if (envval) then
       value = .true.
    else if (present(default)) then
       value = default
    else
       value = defaults(i)%value_l
    end if

    call args_switch(trim(defaults(i)%switches), value=value)

    if ((.not. present(silent)) .or. (.not. silent)) then
       if (value) then
          write(*,*) trim(defaults(i)%descr)
       end if
    end if

  end subroutine tbio_args_l

  !--------------------------------------------------------------------!

  subroutine tbio_args_c(name, value, default, silent)

    implicit none

    character(len=*),           intent(in)  :: name
    character(len=*),           intent(out) :: value
    character(len=*), optional, intent(in)  :: default
    logical,          optional, intent(in)  :: silent
    
    character(len=lenval) :: envval
    integer :: i, stat

    i = def_index(name)

    if (defaults(i)%type /= 'c') then
       write(0,*) "Error: incompatible types in `tbio'."
       call tbio_final()
       stop
    end if

    stat = 1
    if (len_trim(defaults(i)%env) > 0) then
       call args_env(defaults(i)%env, value=envval, stat=stat)
    end if

    if (stat == 0) then
       value = envval
    else if (present(default)) then
       value = default
    else
       value = trim(defaults(i)%value_c)
    end if

    call args_switch(trim(defaults(i)%switches), value=value)

    if ((len_trim(value) > 0) .and. (.not. present(silent))) then
       write(*,*) io_trim(defaults(i)%descr,wdescr), ' : ', &
                  trim(adjustl(value))
    end if

  end subroutine tbio_args_c

  !--------------------------------------------------------------------!

  subroutine tbio_args_d(name, value, default, silent, digits)

    implicit none

    character(len=*),           intent(in)  :: name
    double precision,           intent(out) :: value
    double precision, optional, intent(in)  :: default
    logical,          optional, intent(in)  :: silent
    integer,          optional, intent(in)  :: digits

    character(len=50) :: frmt
    double precision  :: envval
    integer :: i, j, stat

    i = def_index(name)

    if (defaults(i)%type /= 'd') then
       write(0,*) "Error: incompatible types in `tbio'."
       call tbio_final()
       stop
    end if

    stat = 1
    if (len_trim(defaults(i)%env) > 0) then
       call args_env(defaults(i)%env, value=envval, stat=stat)
    end if

    if (stat == 0) then
       value = envval
    else if (present(default)) then
       value = default
    else
       value = defaults(i)%value_d
    end if

    call args_switch(trim(defaults(i)%switches), value=value)

    if (present(digits)) then
       j = digits
    else
       j = 6
    end if

    if ((.not. present(silent)) .or. (.not. silent)) then
       frmt = '(1x,A' // trim(io_adjustl(wdescr)) // '," : ",ES' &
           // trim(io_adjustl(j+7)) // '.' // trim(io_adjustl(j)) // ')'
       write(*,frmt) io_trim(defaults(i)%descr,wdescr), value
    end if

  end subroutine tbio_args_d

  !--------------------------------------------------------------------!
  !           implementation of the set_default() interface            !
  !--------------------------------------------------------------------!

  subroutine set_default_i(name, descr, value, switches, help, env)

    implicit none

    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: descr
    integer,                    intent(in) :: value
    character(len=*),           intent(in) :: switches
    character(len=*),           intent(in) :: help
    character(len=*), optional, intent(in) :: env

    if (idef >= nDefaults) then
       write(0,*) "Error: number of defaults exceeds allocated memory."
       call tbio_final()
       stop
    end if

    idef = idef + 1
    defaults(idef)%name     = trim(adjustl(name))
    defaults(idef)%type     = 'i'
    defaults(idef)%value_i  = value
    defaults(idef)%descr    = trim(adjustl(descr))
    defaults(idef)%switches = trim(adjustl(switches))
    defaults(idef)%help     = trim(adjustl(help))

    if (present(env)) then
       defaults(idef)%env   = trim(adjustl(env))
    end if

  end subroutine set_default_i

  !--------------------------------------------------------------------!

  subroutine set_default_l(name, descr, value, switches, help, env)

    implicit none

    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: descr
    logical,                    intent(in) :: value
    character(len=*),           intent(in) :: switches
    character(len=*),           intent(in) :: help
    character(len=*), optional, intent(in) :: env

    if (idef >= nDefaults) then
       write(0,*) "Error: number of defaults exceeds allocated memory."
       call tbio_final()
       stop
    end if

    idef = idef + 1
    defaults(idef)%name     = trim(adjustl(name))
    defaults(idef)%type     = 'l'
    defaults(idef)%value_l  = value
    defaults(idef)%descr    = trim(adjustl(descr))
    defaults(idef)%switches = trim(adjustl(switches))
    defaults(idef)%help     = trim(adjustl(help))

    if (present(env)) then
       defaults(idef)%env   = trim(adjustl(env))
    end if

  end subroutine set_default_l

  !--------------------------------------------------------------------!

  subroutine set_default_c(name, descr, value, switches, help, env)

    implicit none

    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: descr
    character(len=*),           intent(in) :: value
    character(len=*),           intent(in) :: switches
    character(len=*),           intent(in) :: help
    character(len=*), optional, intent(in) :: env

    if (idef >= nDefaults) then
       write(0,*) "Error: number of defaults exceeds allocated memory."
       call tbio_final()
       stop
    end if

    idef = idef + 1
    defaults(idef)%name     = trim(adjustl(name))
    defaults(idef)%type     = 'c'
    defaults(idef)%value_c  = trim(value)
    defaults(idef)%descr    = trim(adjustl(descr))
    defaults(idef)%switches = trim(adjustl(switches))
    defaults(idef)%help     = trim(adjustl(help))

    if (present(env)) then
       defaults(idef)%env   = trim(adjustl(env))
    end if

  end subroutine set_default_c

  !--------------------------------------------------------------------!

  subroutine set_default_d(name, descr, value, switches, help, env)

    implicit none

    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: descr
    double precision,           intent(in) :: value
    character(len=*),           intent(in) :: switches
    character(len=*),           intent(in) :: help
    character(len=*), optional, intent(in) :: env

    if (idef >= nDefaults) then
       write(0,*) "Error: number of defaults exceeds allocated memory."
       call tbio_final()
       stop
    end if

    idef = idef + 1
    defaults(idef)%name     = trim(adjustl(name))
    defaults(idef)%type     = 'd'
    defaults(idef)%value_d  = value
    defaults(idef)%descr    = trim(adjustl(descr))
    defaults(idef)%switches = trim(adjustl(switches))
    defaults(idef)%help     = trim(adjustl(help))

    if (present(env)) then
       defaults(idef)%env   = trim(adjustl(env))
    end if

  end subroutine set_default_d

  !--------------------------------------------------------------------!
  !           implementation of the tbio_default() interface           !
  !--------------------------------------------------------------------!

  subroutine tbio_default_i(name, value)

    implicit none

    character(len=*), intent(in)  :: name
    integer,          intent(out) :: value

    integer :: i

    i = def_index(name)
    if (defaults(i)%type /= 'i') then
       write(0,*) "Error: Incompatible data types in `tbio' for: ", trim(name)
       call tbio_final()
       stop
    end if

    value = defaults(i)%value_i

  end subroutine tbio_default_i

  !--------------------------------------------------------------------!

  subroutine tbio_default_l(name, value)

    implicit none

    character(len=*), intent(in)  :: name
    logical,          intent(out) :: value

    integer :: i

    i = def_index(name)
    if (defaults(i)%type /= 'l') then
       write(0,*) "Error: Incompatible data types in `tbio' for: ", trim(name)
       call tbio_final()
       stop
    end if

    value = defaults(i)%value_l

  end subroutine tbio_default_l

  !--------------------------------------------------------------------!

  subroutine tbio_default_c(name, value)

    implicit none

    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: value

    integer :: i

    i = def_index(name)
    if (defaults(i)%type /= 'c') then
       write(0,*) "Error: Incompatible data types in `tbio' for: ", trim(name)
       call tbio_final()
       stop
    end if

    value = trim(defaults(i)%value_c)

  end subroutine tbio_default_c

  !--------------------------------------------------------------------!

  subroutine tbio_default_d(name, value)

    implicit none

    character(len=*), intent(in)  :: name
    double precision, intent(out) :: value

    integer :: i

    i = def_index(name)
    if (defaults(i)%type /= 'd') then
       write(0,*) "Error: Incompatible data types in `tbio' for: ", trim(name)
       call tbio_final()
       stop
    end if

    value = defaults(i)%value_d

  end subroutine tbio_default_d


end module tbio
