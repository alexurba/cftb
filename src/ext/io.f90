module io

  !--------------------------------------------------------------------!
  !   I/O procedures for general use                                   !
  !--------------------------------------------------------------------!
  !                                                                    !
  ! character test functions:                                          !
  !                                                                    !
  ! io_isfigure(chr)   : returns .true. if chr is a figure             !
  ! io_isletter(chr)   : returns .true. if chr is a letter             !
  !                                                                    !
  ! String modification:                                               !
  !                                                                    !
  ! io_lower(str)        : returns string str in lower case only       !
  ! io_upper(str)        : returns string str in upper case only       !
  ! io_replace(s1,s2,s3) : substitute s2 by s3 in string s1            !
  ! io_trim(str,n)       : trim string str to length n                 !
  ! io_center(str,n)     : center a string str on a new string with    !
  !                        length n                                    !
  !                                                                    !
  ! General I/O:                                                       !
  !                                                                    !
  ! io_adjustl(x)      : adjustl() function for integer and double     !
  ! io_readval(s,n,v)  : searches string s for a key with name n and   !
  !                      returns its value in v                        !
  ! io_readnext(s,p,v) : read next values from a string starting with  !
  !                      position p                                    !
  !                                                                    !
  ! File I/O:                                                          !
  !                                                                    !
  ! io_unit(u_try)     : returns an unconnected unit number >= u_try   !
  !                      If u_try is not provided, a default value is  !
  !                      used.                                         !
  !                                                                    ! 
  !--------------------------------------------------------------------!
  ! 2010-07-05 Alexander Urban (AU)                                    !
  ! 2011-02-21 AU --- added `io_trim()' and `io_replace'               !
  ! 2011-11-18 AU --- added `io_center()'                              !
  ! 2012-01-19 AU --- added `io_unit()' to request an unconnected unit !
  !--------------------------------------------------------------------!

  implicit none

  public :: io_adjustl,  &
            io_center,   &
            io_isfigure, &
            io_isletter, &
            io_readnext, &
            io_readval,  &
            io_replace,  &
            io_trim,     &
            io_lower,    &
            io_upper      


  !--------------------------------------------------------------------!
  !                       io_adjustl (interface)                       !
  !                                                                    !
  ! Use io_adjustl to adjust not only strings but also different types !
  ! for convenient output.                                             !
  !                                                                    !
  ! Available implementations:                                         !
  !   adjustl_i(int)         -- integer                                !
  !   adjustl_d(dp, digits)  -- double precision (digits is optional)  !
  !--------------------------------------------------------------------!

  interface io_adjustl
     module procedure adjustl_i, adjustl_d
  end interface


  !--------------------------------------------------------------------!
  !                      io_readnext (interface)                       !
  !                                                                    !
  ! Use the readnext interface to read the next entry/value from a     !
  ! string, given a starting position pos:                             !
  !                                                                    !
  ! input:   string = '  1  2  3  '; pos = 1                           !
  !          call io_readnext(string, pos, val)                        !
  !                                                                    !
  ! This will result in pos == 4 (right after the value) and val == 1  !
  !                                                                    !
  ! Available implementations:                                         !
  !  readnext_i1(string, pos, value)    -- single integer value        !
  !  readnext_in(string, pos, value, n) -- n integer values            !
  !  readnext_d1(string, pos, value)    -- single double prec. value   !
  !  readnext_dn(string, pos, value, n) -- n double precision values   !
  !  readnext_c1(string, pos, value)    -- single character string     !
  !--------------------------------------------------------------------!

  interface io_readnext
     module procedure readnext_i1, readnext_in, readnext_d1, &
                      readnext_dn, readnext_c1
  end interface


  !--------------------------------------------------------------------!
  !                      io_readval (interface)                        !
  !                                                                    !
  ! The readval interface can be used to search a character string for !
  ! key/value pairs, separated by an equal (=) sign.                   !
  !   Example:   string='  natoms=3   name = O  '                      !
  !              call io_readval(string, 'name', atom_type)            !
  ! The interface expects the search string `string', the name of the  !
  ! key 'name', and a suitable variable 'value'.                       !
  !                                                                    !
  ! Available implementations:                                         !
  !  readval_i1(string, name, value)    -- single integer value        !
  !  readval_in(string, name, value, n) -- n integer values            !
  !  readval_d1(string, name, value)    -- single double prec. value   !
  !  readval_dn(string, name, value, n) -- n double precision values   !
  !  readval_c1(string, name, value)    -- single character string     !
  !--------------------------------------------------------------------!

  interface io_readval
     module procedure readval_i1, readval_in, readval_d1, readval_dn, &
                      readval_c1, readval_l
  end interface

contains


  !--------------------------------------------------------------------!
  !                      adjustl (implementation)                      !
  !                                                                    !
  !          Please use the module interface `io_adjustl()'.           !
  !--------------------------------------------------------------------!

  function adjustl_i(int) result(str)

    implicit none

    integer, intent(in) :: int
    character(len=50)   :: str

    write(str, *) int
    str = trim(adjustl(str))

  end function adjustl_i

  !--------------------------------------------------------------------!

  function adjustl_d(dp, digits) result(str)

    implicit none

    double precision,            intent(in) :: dp
    integer,           optional, intent(in) :: digits
    
    character(len=50)            :: frmt, str

    if (present(digits)) then
       write(frmt, *) digits
       frmt = '(F50.' // trim(adjustl(frmt)) // ')'
    else
       frmt = '(F50.2)'
    end if

    write(str, frmt) dp
    str = trim(adjustl(str))

  end function adjustl_d

  !--------------------------------------------------------------------!
  !                              io_trim                               !
  !--------------------------------------------------------------------!

  function io_trim(str, n) result(out)

    implicit none

    character(len=*), intent(in) :: str
    integer,          intent(in) :: n
    character(len=n)             :: out

    integer :: i

    out = trim(str)

    do i = len_trim(out) + 1, n
       out = out // ' '
    end do

  end function io_trim

  !--------------------------------------------------------------------!
  !                             io_replace                             !
  !                                                                    !
  ! str = io_replace(str, sub, ins)                                    !
  ! --> will substitute every occurance of the string `sub' in the     !
  !     longer string `str' by the string `ins'                        !
  !--------------------------------------------------------------------!

  function io_replace(str, sub, ins) result(out)

    implicit none

    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: sub
    character(len=*), intent(in) :: ins

    character(len=len(str))      :: out

    integer :: i1, i2, j1, i, j
    integer :: lsub, lstr, lins
    
    lsub = len(sub)
    lstr = len(str)
    lins = len(ins)

    out = ''

    i1 = 1
    j1 = 1
    do
       i2 = index(str(i1:lstr), sub)
       if (i2 == 0) exit

       i = min(lstr, i1 + i2 - 2)
       j = min(lstr, j1 + i2 - 2)

       out(j1:j) = str(i1:i)
       out(j+1:j+lins) = ins

       i1 = i + lsub + 1
       j1 = j + lins + 1
       if ((i1 > lstr) .or. (i2 > lstr)) exit
    end do

    out(j1:lstr) = str(i1:lstr)   

  end function io_replace

  !--------------------------------------------------------------------!
  !                     readnext (implementation)                      !
  !                                                                    !
  !          Please use the module interface `io_readnext()'.          !
  !--------------------------------------------------------------------!

  subroutine readnext_i1(string, pos, value)

    implicit none

    character(len=*), intent(in)    :: string
    integer,          intent(inout) :: pos
    integer,          intent(out)   :: value

    integer :: i1, i2

    i2 = len_trim(string)
    
    if (pos > i2) then
       pos   = 0
       value = 0
       return
    end if

    i1 = pos
    do
       if (string(i1:i1) /= ' ') exit
       i1 = i1 + 1
    end do

    read(string(i1:i2), *) value
    pos = scan(string(i1:i2), ' ')
    if (pos == 0) then
       pos = i2 + 1
    else
       pos = pos + i1 - 1
    end if

  end subroutine readnext_i1

  !--------------------------------------------------------------------!

  subroutine readnext_in(string, pos, value, n)

    implicit none

    integer,               intent(in)    :: n
    character(len=*),      intent(in)    :: string
    integer,               intent(inout) :: pos
    integer, dimension(n), intent(out)   :: value

    integer :: i

    value(1:n) = 0

    if (pos > len_trim(string)) then
       pos        = 0
       return
    end if

    do i = 1, n
       call readnext_i1(string, pos, value(i))
    end do

  end subroutine readnext_in

  !--------------------------------------------------------------------!

  subroutine readnext_d1(string, pos, value)

    implicit none

    character(len=*), intent(in)    :: string
    integer,          intent(inout) :: pos
    double precision, intent(out)   :: value

    integer :: i1, i2

    i2 = len_trim(string)
    
    if (pos > i2) then
       pos   = 0
       value = 0.0d0
       return
    end if

    i1 = pos
    do
       if (string(i1:i1) /= ' ') exit
       i1 = i1 + 1
    end do

    read(string(i1:i2), *) value
    pos = scan(string(i1:i2), ' ')
    if (pos == 0) then
       pos = i2 + 1
    else
       pos = pos + i1 - 1
    end if

  end subroutine readnext_d1

  !--------------------------------------------------------------------!

  subroutine readnext_dn(string, pos, value, n)

    implicit none

    integer,                        intent(in)    :: n
    character(len=*),               intent(in)    :: string
    integer,                        intent(inout) :: pos
    double precision, dimension(n), intent(out)   :: value

    integer :: i

    value(1:n) = 0.0d0

    if (pos > len_trim(string)) then
       pos = 0
       return
    end if

    do i = 1, n
       call readnext_d1(string, pos, value(i))
    end do

  end subroutine readnext_dn

  !--------------------------------------------------------------------!

  subroutine readnext_c1(string, pos, value)

    implicit none

    character(len=*), intent(in)    :: string
    integer,          intent(inout) :: pos
    character(len=*), intent(out)   :: value

    integer :: i1, i2

    i2 = len_trim(string)
    
    if (pos > i2) then
       pos   = 0
       value = ' '
       return
    end if

    i1 = pos
    do
       if (string(i1:i1) /= ' ') exit
       i1 = i1 + 1
    end do

    read(string(i1:i2), *) value

    pos = scan(string(i1:i2), ' ')
    if (pos == 0) then
       pos = i2 + 1
    else
       pos = pos + i1 - 1
    end if


  end subroutine readnext_c1

  !--------------------------------------------------------------------!
  !                      readval (implementation)                      !
  !                                                                    !
  !          Please use the module interface `io_readval()'.           !
  !--------------------------------------------------------------------!

  subroutine readval_i1(string, name, value)

    implicit none

    character(len=*), intent(in)    :: string, name
    integer,          intent(inout) :: value

    integer :: slen
    integer :: i, j

    slen = len_trim(string)

    i = index(string, trim(name))
    if (i == 0) then
       return
    end if
    i = i + len_trim(name)
    j = index(string(i:slen), '=')
    if (j == 0) then
       write(0,*) 'Error: invalid string in readval.'
       write(0,*) '       ', string
       stop
    end if
    i = i + j
    read(string(i:slen), *) value

  end subroutine readval_i1

  !--------------------------------------------------------------------!

  subroutine readval_in(string, name, value, n)

    implicit none

    integer,               intent(in)    :: n
    character(len=*),      intent(in)    :: string, name
    integer, dimension(n), intent(inout) :: value

    integer :: slen
    integer :: i, j

    slen = len_trim(string)

    i = index(string, trim(name))
    if (i == 0) then
       return
    end if
    i = i + len_trim(name)
    j = index(string(i:slen), '=')
    if (j == 0) then
       write(0,*) 'Error: invalid string in readval.'
       write(0,*) '       ', string
       stop
    end if
    i = i + j
    read(string(i:slen), *) value(1:n)

  end subroutine readval_in

  !--------------------------------------------------------------------!

  subroutine readval_d1(string, name, value)

    implicit none

    character(len=*), intent(in)    :: string, name
    double precision, intent(inout) :: value

    integer :: slen
    integer :: i, j

    slen = len_trim(string)

    i = index(string, trim(name))
    if (i == 0) then
       return
    end if
    i = i + len_trim(name)
    j = index(string(i:slen), '=')
    if (j == 0) then
       write(0,*) 'Error: invalid string in readval.'
       write(0,*) '       ', string
       stop
    end if
    i = i + j
    read(string(i:slen), *) value

  end subroutine readval_d1

  !--------------------------------------------------------------------!

  subroutine readval_dn(string, name, value, n)

    implicit none

    integer,                        intent(in)    :: n
    character(len=*),               intent(in)    :: string, name
    double precision, dimension(n), intent(inout) :: value

    integer :: slen
    integer :: i, j

    slen = len_trim(string)

    i = index(string, trim(name))
    if (i == 0) then
       return
    end if
    i = i + len_trim(name)
    j = index(string(i:slen), '=')
    if (j == 0) then
       write(0,*) 'Error: invalid string in readval.'
       write(0,*) '       ', string
       stop
    end if
    i = i + j
    read(string(i:slen), *) value(1:n)

  end subroutine readval_dn

  !--------------------------------------------------------------------!

  subroutine readval_c1(string, name, value)

    implicit none

    character(len=*), intent(in)    :: string, name
    character(len=*), intent(inout) :: value

    integer :: slen
    integer :: i, j

    slen = len_trim(string)

    i = index(string, trim(name))
    if (i == 0) then
       return
    end if
    i = i + len_trim(name)
    j = index(string(i:slen), '=')
    if (j == 0) then
       write(0,*) 'Error: invalid string in readval.'
       write(0,*) '       ', string
       stop
    end if
    i = i + j
    read(string(i:slen), *) value

  end subroutine readval_c1

  !--------------------------------------------------------------------!

  subroutine readval_l(string, name, value)

    implicit none

    character(len=*), intent(in)  :: string, name
    logical,          intent(out) :: value

    integer :: i

    i = index(string, trim(name))
    if (i == 0) then
       value = .false.
    else
       value = .true.
    end if

  end subroutine readval_l

  !--------------------------------------------------------------------!
  !                      Character test functions                      !
  !--------------------------------------------------------------------!

  function io_isletter(chr) result(check)

    character, intent(in) :: chr
    logical               :: check
    integer               :: ichk

    character(len=*), parameter :: letters='abcdefghijklmnopqrstuvwxyz' &
                                        // 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    ichk = scan(letters,chr)
    if (ichk == 0) then
       check = .false.
    else
       check = .true.
    end if

  end function io_isletter

  !--------------------------------------------------------------------!

  function io_isfigure(chr) result(check)

    character, intent(in) :: chr
    logical               :: check
    integer               :: ichk

    character(len=*), parameter :: figures='1234567890'

    ichk = scan(figures,chr)
    if (ichk == 0) then
       check = .false.
    else
       check = .true.
    end if

  end function io_isfigure

  !--------------------------------------------------------------------!
  !                  Lower and upper case conversion                   !
  !--------------------------------------------------------------------!

  function io_lower(str_in) result(str_out)

    implicit none

    character(len=*), intent(in) :: str_in
    character(len=len(str_in))   :: str_out

    integer, parameter :: ilowerA = ichar('a')
    integer, parameter :: iupperA = ichar('A')
    integer, parameter :: iupperZ = ichar('Z')

    integer :: i, ichr, nchr, iconv

    iconv = ilowerA - iupperA

    nchr = len(str_in)
    do i = 1, nchr
       ichr = ichar(str_in(i:i))
       if ((ichr >= iupperA) .and. (ichr <= iupperZ)) then
          str_out(i:i) = char(ichr + iconv)
       else
          str_out(i:i) = str_in(i:i)
       end if
    end do

  end function io_lower

  !--------------------------------------------------------------------!

  function io_upper(str_in) result(str_out)

    implicit none

    character(len=*), intent(in) :: str_in
    character(len=len(str_in))   :: str_out

    integer, parameter :: ilowerA = ichar('a')
    integer, parameter :: ilowerZ = ichar('z')
    integer, parameter :: iupperA = ichar('A')

    integer :: i, ichr, nchr, iconv

    iconv = iupperA - ilowerA

    nchr = len(str_in)
    do i = 1, nchr
       ichr = ichar(str_in(i:i))
       if ((ichr >= ilowerA) .and. (ichr <= ilowerZ)) then
          str_out(i:i) = char(ichr + iconv)
       else
          str_out(i:i) = str_in(i:i)
       end if
    end do

  end function io_upper

  !--------------------------------------------------------------------!
  !                    io_center - center a string                     !
  !--------------------------------------------------------------------!

  function io_center(str, n) result(str2)

    implicit none
    
    character(len=*), intent(in) :: str
    integer,          intent(in) :: n
    character(len=n)             :: str2

    integer :: l1

    str2 = ""
    l1 = len_trim(adjustl(str))
    if (l1 > n) return

    l1 = (n - l1)/2

    str2 = repeat(' ',l1) // trim(adjustl(str))

  end function io_center

  !--------------------------------------------------------------------!
  !             io_unit - request unconnected unit number              !
  !--------------------------------------------------------------------!

  function io_unit(u_try) result(u)

    implicit none

    !------------------------------------------------------------------!
    ! If the optional argument 'u_try' is provided, its value will be  !
    ! the first unit number to check.  Otherwise a default value (2)   !
    ! is used.  The unit number is successively increases by 1 until   !
    ! an unconnected unit is found.                                    !
    !------------------------------------------------------------------!

    integer, optional, intent(in) :: u_try
    integer                       :: u
    integer,            parameter :: u_ini = 20
    logical                       :: uexists, uopened

    if (present(u_try)) then
       u = u_try
    else
       u = u_ini
    end if

    search : do
       inquire(unit=u, exist=uexists)
       if (uexists) then
          inquire(unit=u, opened=uopened)
          if (.not. uopened) exit search
       end if
       u = u + 1
       if (u >= 100) then
          write(0,*) "Error: unable to find unused unit u < 100. (io_unit)"
          stop
       end if
    end do search

  end function io_unit


end module io
