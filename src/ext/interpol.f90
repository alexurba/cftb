!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Library with interpolation procedures.                             !!
!! ------------------------------------------------------------------ !!
!! Interpolation with cubic splines:                                  !!
!!  - deriv2()     compute 2nd derivative                             !!
!!  - cs()         compute interpolated value (needs 2nd derivative)  !!
!! Convenience interface:                                             !!
!!  - cs_init()    initialize cubic spline interpolator               !!
!!  - cs_final()   finalize cubic splines interpolator                !!
!!  - cs_interp()  the interpolation function (user interface)        !!
!! ------------------------------------------------------------------ !!
!! 2010-03-28 Alexander Urban (AU)                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module interpol

  public :: cs_init,   &
            cs_final,  &
            cs_interp, &
            cs,        &
            cs_d1,     &
            cs_d2,     &
            deriv2


  double precision, parameter,                 public  :: yp_max = 1.0d30

  !--------------------------------------------------------------------!
  ! x_p, y_p     : pointers to the tabulated function values           !
  ! y2(ix)       : second derviatives                                  !
  ! nx           : number of values; 1 <= ix <= nx                     !
  !--------------------------------------------------------------------!

  double precision, dimension(:), pointer,     private :: x_p, y_p
  double precision, dimension(:), allocatable, private :: y2l
  integer,                                     private :: nx = 0 

  interface cs_init
     module procedure cs_init1, cs_init2
  end interface

contains

  !--------------------------------------------------------------------!
  !                 Computation of the 2nd Derivative                  !
  !--------------------------------------------------------------------!

  subroutine deriv2(x, y, n, yp1, ypn, y2)
    !! compute the second derivatives of a tabulated function
    !! adapted from Press et al., Numerical recipes in F77, "spline"

    implicit none

    !------------------------------------------------------------------!
    ! x(ix)      : ix-th x-value of the tabulated function             !
    ! y(ix)      : function value at x = x(ix)                         !
    ! n          : dimension; 1 <= ix <= n                             !
    ! yp1, ypn   : 1st derivatives of y(x) at boundaries x(1) and x(n) !
    ! y2(ix)     : second derivative at x = x(ix)                      !
    !------------------------------------------------------------------!

    integer,                         intent(in)  :: n
    double precision, dimension(n),  intent(in)  :: x, y
    double precision,                intent(in)  :: yp1, ypn
    double precision, dimension(n),  intent(out) :: y2
    
    integer                                      :: i, k
    double precision                             :: p, qn, sig, un
    double precision, dimension(:), allocatable  :: u

    allocate(u(n))

    if (yp1 >= yp_max) then
       ! lower natural boundary conditions, if the 1st derivative at
       ! x(1) is larger than 1x10^30:
       y2(1) = 0.0d0
       u(1)  = 0.0d0
    else
       y2(1) = -0.5d0
       u(1)  = (3.0d0/(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
    endif

    do i = 2, n-1
       sig   = (x(i) - x(i-1))/(x(i+1) - x(i-1))
       p     = sig*y2(i-1) + 2.0d0
       y2(i) = (sig - 1.0d0)/p
       u(i)  = (6.0d0*((y(i+1) - y(i))/(x(i+1) - x(i)) &
                     - (y(i) - y(i-1))/(x(i) - x(i-1)) &
               )/(x(i+1) - x(i-1)) - sig*u(i-1))/p
    end do
  
    if (ypn >= yp_max) then
       ! lower natural boundary conditions, if the 1st derivative at
       ! x(n) is larger than 1x10^30:
       qn = 0.0d0
       un = 0.0d0
    else
       qn = 0.5d0
       un = (3.0d0/(x(n) - x(n-1)))*(ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))
    endif

    y2(n) = ( un - qn * u(n-1) ) / ( qn * y2(n-1) + 1.0d0 )

    do k = n-1, 1, -1
       y2(k) = y2(k) * y2(k+1) + u(k)
    end do

    deallocate(u)

  end subroutine deriv2

  !--------------------------------------------------------------------!
  !                  Interpolation with Cubic Splines                  !
  !--------------------------------------------------------------------!

  subroutine cs_init1(x, y, n, yp1, ypn)
    !! initialize cubic splines interpolator: allocate memory and
    !! connect pointers

    implicit none

    !------------------------------------------------------------------!
    ! x, y     : tabulated function y(x(i)), for 1 <= n <= n           !
    ! n        : number of tabulated function values                   !
    ! yp1, ypn : 1st derivative y'(x) at the boundaries x(1) and x(n)  !
    !------------------------------------------------------------------!

    double precision, dimension(:), target, intent(in) :: x, y
    integer,                                intent(in) :: n
    double precision,                       intent(in) :: yp1, ypn
    
    nx = n
    allocate(y2l(nx))

    x_p => x(:)
    y_p => y(:)

    call deriv2(x_p, y_p, nx, yp1, ypn, y2l)

  end subroutine cs_init1

  !--------------------------------------------------------------------!

  subroutine cs_init2(x, y, n)
    !! initialize cubic splines interpolator: allocate memory and
    !! connect pointers

    implicit none

    !------------------------------------------------------------------!
    ! x, y     : tabulated function y(x(i)), for 1 <= n <= n           !
    ! n        : number of tabulated function values                   !
    !------------------------------------------------------------------!

    double precision, dimension(:), target, intent(in) :: x, y
    integer,                                intent(in) :: n

    double precision, parameter :: yp_max = 1.0d30

    nx = n
    allocate(y2l(nx))

    x_p => x(:)
    y_p => y(:)

    call deriv2(x_p, y_p, nx, yp_max, yp_max, y2l)

  end subroutine cs_init2

  !--------------------------------------------------------------------!

  subroutine cs_final()
    !! finalize the cubic splines interpolator:
    !! free memory and disconnect pointers

    implicit none

    if (allocated(y2l)) deallocate(y2l)
    nx = 0

    x_p => null()
    y_p => null()

  end subroutine cs_final

  !--------------------------------------------------------------------!
  !                                I/O                                 !
  !--------------------------------------------------------------------!

  subroutine cs_save(x, y, y1_1, y1_n, n, filename, replace)

    implicit none

    integer,                        intent(in) :: n
    double precision, dimension(n), intent(in) :: x, y
    double precision,               intent(in) :: y1_1, y1_n
    character(len=*),               intent(in) :: filename
    logical,          optional,     intent(in) :: replace

    integer, parameter :: u_out = 60
    integer            :: i
    logical            :: fexists

    if ((.not. present(replace)) .or. (.not. replace)) then
       inquire(file=trim(filename), exist=fexists)
       if (fexists) then
          write(0,*) 'Error: file already exists: ', trim(filename)
          return
       end if
    end if

    open(u_out, file=trim(filename), status='replace', action='write')

    write(u_out,'("#",3(A15,2x))') ' x             ', &
         ' y              ', " y'             "

    write(u_out,'(1x,3(ES15.8,2x))') x(1), y(1), y1_1
    do i = 2, n - 1
       write(u_out,'(1x,2(ES15.8,2x))') x(i), y(i)
    end do
    write(u_out,'(1x,3(ES15.8,2x))') x(n), y(n), y1_n

    close(u_out)

  end subroutine cs_save

  !--------------------------------------------------------------------!

  subroutine cs_load(filename, x, y, y1_1, y1_n, n)

    implicit none

    !------------------------------------------------------------------!
    ! for n<=0 --> just count the number of abscissas                  !
    !------------------------------------------------------------------!

    integer,                        intent(inout) :: n
    character(len=*),               intent(in)    :: filename
    double precision, dimension(n), intent(out)   :: x, y
    double precision,               intent(out)   :: y1_1, y1_n

    integer, parameter :: u_in = 50

    logical            :: fexists
    character(len=128) :: line
    integer            :: nnodes, inode
    integer            :: status

    inquire(file=trim(filename), exist=fexists)
    if (.not. fexists) then
       write(0,*) 'Error: file not found: ', trim(filename)
       return
    end if

    open(u_in, file=trim(filename), status='old', action='read')

    nnodes = 0
    countloop : do
       read(u_in, '(A)', iostat=status) line
       if (status /= 0) exit countloop
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle countloop
       if (len_trim(line) == 0) cycle countloop
       nnodes = nnodes + 1
    end do countloop

    if (n <= 0) then
       n = nnodes
    else if (n < nnodes) then
       write(0,*) "Error: arrays too small in `cs_load()'."
    else
       rewind(u_in)
       inode = 0
       do while (inode < nnodes)
          read(u_in, '(A)', iostat=status) line
          if (status /= 0) exit
          line = trim(adjustl(line))
          if (line(1:1) == '#') cycle 
          if (len_trim(line) == 0) cycle
          inode = inode + 1
          if (inode == 1) then
             read(line, *) x(inode), y(inode), y1_1
          else if (inode == nnodes) then
             read(line, *) x(inode), y(inode), y1_n
          else
             read(line, *) x(inode), y(inode)
          end if
       end do
    end if

    close(u_in)

  end subroutine cs_load

  !--------------------------------------------------------------------!
  !                 interpolation convenience function                 !
  !--------------------------------------------------------------------!

  function cs_interp(x0) result(y0)
    !! interpolate the tabulated function at point x0
    !! --> convenience interface; if speed is desired directly call cs()

    implicit none

    double precision,                intent(in) :: x0
    double precision                            :: y0

    if (allocated(y2l)) then
       call cs(x_p, y_p, y2l, nx, x0, y0)
    else
       write(0,*) 'Error: Interpolation module not initialized!'
       y0 = 0.0d0
    end if

  end function cs_interp
  
  !--------------------------------------------------------------------!
  !                   reduce the number of abscissas                   !
  !--------------------------------------------------------------------!

  subroutine cs_simplify(x, y, y2, n, xs, ys, y2s, ns, &
                         mad_max, rmsd_max, mad, rmsd)

    implicit none

    integer,                        intent(in)  :: n
    double precision, dimension(n), intent(in)  :: x, y, y2
    double precision, dimension(n), intent(out) :: xs, ys, y2s
    integer,                        intent(out) :: ns
    double precision,               intent(in)  :: mad_max, rmsd_max
    double precision,               intent(out) :: mad, rmsd

    double precision, dimension(n) :: xw, yw, y2w
    double precision               :: r, dx, x0, x0s
    double precision               :: min_slope, slope, diff
    double precision               :: y1_1, y1_n
    double precision               :: y1_a, y1_b, y1_c
    integer                        :: nsteps, i, xi, nw
    double precision               :: mad_prev, rmsd_prev

    xw(:)  = x(:)
    yw(:)  = y(:)
    y2w(:) = y2(:)
    nw     = n
    call cs_d1(x,y,y2,n, x(1), y1_1)
    call cs_d1(x,y,y2,n, x(n), y1_n)
    
    nsteps = 2*n
    dx = (x(n) - x(1))/dble(nsteps - 1)

    rmsd      = 0.0d0
    mad       = 0.0d0
    rmsd_prev = 0.0d0
    mad_prev  = 0.0d0
    do while ((mad <= mad_max) .and. (rmsd <= rmsd_max) .and. (nw > 4))
       
       ns = nw
       xs(1:ns)  = xw(1:ns)
       ys(1:ns)  = yw(1:ns)
       y2s(1:ns) = y2w(1:ns)

       ! search for three abscissas to be reduced to two:
       call cs_d1(xs,ys,y2s,ns, xs(2), y1_a)
       call cs_d1(xs,ys,y2s,ns, xs(3), y1_b)
       call cs_d1(xs,ys,y2s,ns, xs(4), y1_c)
       min_slope = abs(y1_b-y1_a) + abs(y1_b-y1_c)
       xi = 3
       do i = 4, ns-2
          y1_a = y1_b
          y1_b = y1_c
          call cs_d1(xs,ys,y2s,ns, xs(i+1), y1_c)
          slope = abs(y1_b-y1_a) + abs(y1_b-y1_c)
          if (slope < min_slope) then
             min_slope = slope
             xi = i
          end if
       end do

       ! remove one abscissa:
       xw(xi-1) = 0.5d0*(xs(xi-1) + xs(xi))
       xw(xi)   = 0.5d0*(xs(xi) + xs(xi+1))
! Note: the use of the reduced set for interpolation will have
!       a smoothing effect on the resulting function
!
!       call cs(x,y,y2,n, xw(xi-1), yw(xi-1))
!       call cs(x,y,y2,n, xw(xi),   yw(xi))
       call cs(xs,ys,y2s,ns, xw(xi-1), yw(xi-1))
       call cs(xs,ys,y2s,ns, xw(xi),   yw(xi))
       xw(xi+1:ns-1) = xs(xi+2:ns)
       yw(xi+1:ns-1) = ys(xi+2:ns)
       nw = ns - 1
       call deriv2(xw, yw, nw, y1_1, y1_n, y2w)

       ! calculate MAD and RMSD:
       mad_prev  = mad
       rmsd_prev = rmsd
       rmsd = 0.0d0
       mad  = 0.0d0
       r = x(1)
       do i = 1, nsteps
          call cs(x,y,y2,n, r, x0)
          call cs(xw,yw,y2w,nw, r, x0s)
          diff = x0 - x0s
          rmsd = rmsd + diff*diff
          mad  = max(mad, abs(diff))
          r = r + dx
       end do
       rmsd = sqrt(rmsd/dble(nsteps))

    end do

    if (ns > 4) then
       mad  = mad_prev
       rmsd = rmsd_prev
    end if
    
  end subroutine cs_simplify

  !--------------------------------------------------------------------!
  !                     main cubic spline routine                      !
  !--------------------------------------------------------------------!

  subroutine cs(x, y, y2, n, x0, y0)
    !! interpolation with cubic splines at (x0,y0)
    !! adapted from Press et al., Numerical recipes in F77, "splint"


    implicit none

    !------------------------------------------------------------------!
    ! x(ix)   : ix-th x-value of the tabulated function                !
    ! y(ix)   : function value for x = x(ix)                           !
    ! y2(ix)  : 2nd derivative for x = x(ix)                           !
    ! n       : dimension; size of function table                      !
    ! x0      : x-value for the interpolation                          !
    ! y0      : interpolated function value at x0                      !
    !------------------------------------------------------------------!

    integer,                        intent(in)  :: n
    double precision, dimension(n), intent(in)  :: x, y, y2
    double precision,               intent(in)  :: x0
    double precision,               intent(out) :: y0

    integer                                     :: k, khi, klo
    double precision                            :: a, b, h 

    klo = 1
    khi = n

    do while (khi - klo > 1)
       k = (khi+klo)/2
       if ( x(k) > x0 ) then
          khi = k
       else
          klo = k
       endif
    end do

    h = x(khi) - x(klo)

    if(h == 0.0d0) then
       write(0,*) 'Error: redundant x value in table (cs).'
       call cs_final()
       stop
    end if

    a = (x(khi) - x0) / h
    b = (x0 - x(klo)) / h

    y0 = a*y(klo) + b*y(khi) &
       + ((a*a*a - a) * y2(klo) + (b*b*b - b)*y2(khi)) &
       *(h*h)/6.0d0

  end subroutine cs

  !--------------------------------------------------------------------!

  subroutine cs_d1(x, y, y2, n, x0, y0)
    !! analytic 1st derivative of the splines


    implicit none

    !------------------------------------------------------------------!
    ! x(ix)   : ix-th x-value of the tabulated function                !
    ! y(ix)   : function value for x = x(ix)                           !
    ! y2(ix)  : 2nd derivative for x = x(ix)                           !
    ! n       : dimension; size of function table                      !
    ! x0      : x-value for the interpolation                          !
    ! y0      : interpolated function value at x0                      !
    !------------------------------------------------------------------!

    integer,                        intent(in)  :: n
    double precision, dimension(n), intent(in)  :: x, y, y2
    double precision,               intent(in)  :: x0
    double precision,               intent(out) :: y0

    integer                                     :: k, khi, klo
    double precision                            :: a, b, h 

    klo = 1
    khi = n

    do while (khi - klo > 1)
       k = (khi+klo)/2
       if ( x(k) > x0 ) then
          khi = k
       else
          klo = k
       endif
    end do

    h = x(khi) - x(klo)

    if(h == 0.0d0) then
       write(0,*) 'Error: redundant x value in table (cs).'
       call cs_final()
       stop
    end if

    a = (x(khi) - x0) / h
    b = (x0 - x(klo)) / h

    y0 = ( y(khi) - y(klo) )/h &
       + ( (-3.0d0*a*a + 1.0d0)*y2(klo) + (3.0d0*b*b - 1.0d0)*y2(khi) ) &
       *h/6.0d0

  end subroutine cs_d1

  !--------------------------------------------------------------------!

  subroutine cs_d2(x, y2, n, x0, y0)
    !! analytic 2nd derivative of the splines
    !!
    !! attention: normally it makes not much sense to use the second
    !! derivative of the cubic spline functions since it is just a
    !! linear interpolation of the true second derivative of the
    !! interpolated function !!!


    implicit none

    !------------------------------------------------------------------!
    ! x(ix)   : ix-th x-value of the tabulated function                !
    ! y2(ix)  : 2nd derivative for x = x(ix)                           !
    ! n       : dimension; size of function table                      !
    ! x0      : x-value for the interpolation                          !
    ! y0      : interpolated function value at x0                      !
    !------------------------------------------------------------------!

    integer,                        intent(in)  :: n
    double precision, dimension(n), intent(in)  :: x, y2
    double precision,               intent(in)  :: x0
    double precision,               intent(out) :: y0

    integer                                     :: k, khi, klo
    double precision                            :: b, h 

    klo = 1
    khi = n

    do while (khi - klo > 1)
       k = (khi+klo)/2
       if ( x(k) > x0 ) then
          khi = k
       else
          klo = k
       endif
    end do

    h = x(khi) - x(klo)

    if(h == 0.0d0) then
       write(0,*) 'Error: redundant x value in table (cs).'
       call cs_final()
       stop
    end if

    ! a = (x(khi) - x0) / h
    b = (x0 - x(klo)) / h

    y0 = y2(klo) + b*( y2(khi) - y2(klo) )

  end subroutine cs_d2

end module interpol
