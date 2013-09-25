program fit

  !--------------------------------------------------------------------!
  ! Fit datapoints with an analytical function.                        !
  !--------------------------------------------------------------------!
  ! 2010-11-30 Alexander Urban                                         !
  !--------------------------------------------------------------------!


  use bondint,   only: polyexp, polyexp_cut

  use tbio,      only: tbio_init,            &
                       tbio_final,           &
                       tbio_args

  use simplex,   only: simplex_init,         &
                       simplex_final,        &
                       simplex_step,         &
                       simplex_reset

  !--------------------------------------------------------------------!

  implicit none

  integer,          parameter :: u_data  = 20
  integer,          parameter :: u_deriv = 21
  integer,          parameter :: u_in    = 22
  integer,          parameter :: u_out   = 23
  integer,          parameter :: u_log   = 24
  character(len=*), parameter :: logfile = 'rmsd.log'

  integer                                     :: np, nfix, ifix
  double precision, dimension(:), allocatable :: param, pini, popt
  logical,          dimension(:), allocatable :: fix
  integer                                     :: ndata
  double precision, dimension(:), allocatable :: xdata, ydata
  integer                                     :: nderiv
  double precision, dimension(:), allocatable :: xderiv, yderiv

  character(len=100) :: infile, outfile, ftype, optmeth
  character(len=100) :: datfile, derivfile, plotfile
  character(len=100) :: line
  logical            :: usederiv

  integer            :: i, ip, iter, niter
  double precision   :: rmsd_target, rmsd0, rmsd1
  logical            :: conv
  integer            :: stat

  integer            :: irep
  integer, parameter :: nrep = 3

  !----------------------------- options ------------------------------!

  call tbio_init(name  = 'fit.x',                                      &
                 descr = 'fit data points with a non-linear function', & 
                 args  = (/ 'in   ', 'out  ', 'funct', 'opt  ',        &
                            'data ', 'deriv', 'niter', 'tol  ' /),     &
                 extra = '')

  call tbio_args('in',     infile)
  call tbio_args('out',    outfile, default='fit.out')
  call tbio_args('funct',  ftype)
  call tbio_args('opt',    optmeth)
  call tbio_args('data',   datfile)
  call tbio_args('deriv',  derivfile)
  call tbio_args('niter',  niter)
  call tbio_args('tol',    rmsd_target)

  call tbio_final()

  if (len_trim(derivfile) == 0) then
     usederiv = .false.
  else
     usederiv = .true.
  end if

  plotfile = trim(outfile) // '.plt'

  !-------------------------- initialization --------------------------!

  ! read data points:
  open(u_data, file=trim(datfile), status='old', action='read')
  ndata = 0
  do
     read(u_data, '(A)', iostat=stat) line
     if (stat /= 0) exit
     if (len_trim(line) == 0) cycle
     ndata = ndata + 1
  end do
  rewind(u_data)
  allocate(xdata(ndata), ydata(ndata))
  do i = 1, ndata
     read(u_data, *) xdata(i), ydata(i)
  end do
  close(u_data)

  ! if used: read derivative points:
  nderiv = 0
  if (usederiv) then
     open(u_deriv, file=trim(derivfile), status='old', action='read')
     do
        read(u_deriv, '(A)', iostat=stat) line
        if (stat /= 0) exit
        if (len_trim(line) == 0) cycle
        nderiv = nderiv + 1
     end do
     rewind(u_deriv)
     allocate(xderiv(nderiv), yderiv(nderiv))
     do i = 1, nderiv
        read(u_deriv, *) xderiv(i), yderiv(i)
     end do
     close(u_deriv)
  end if

  ! read fit.x input file:
  open(u_in, file=trim(infile), status='old', action='read')
  read(u_in, *) np
  allocate(pini(np), fix(np))
  nfix = 0
  do i = 1, np
     read(u_in,*) pini(i), ifix
     if (ifix == 1) then
        nfix = nfix + 1
        fix(i) = .true.
     else
        fix(i) = .false.
     end if
  end do
  close(u_in)

  allocate(param(np-nfix), popt(np-nfix))
  ip = 0
  do i = 1, np
     if (.not. fix(i)) then
        ip = ip + 1
        param(ip) = pini(i)
     end if
  end do

  open(u_log, file=trim(logfile), status='replace', action='write')
  write(u_log, '("#",A12,2x,A12,2x,A12)') 'RMSD (sum)  ', 'RMSD (data) ', &
                                          'RMSD (deriv)'

  !---------------------------- regression ----------------------------!

  call simplex_init(np-nfix)

  iter = 0
  repeat : do irep = 1, nrep
     call simplex_reset()
     conv = .false.
     rmsd0 = 0.0d0
     do i = 1, (niter - iter)
        iter = iter + 1
        rmsd1 = get_rmsd(ndata, xdata, ydata, usederiv, nderiv,  &
                         xderiv, yderiv, ftype, np, nfix, param, &
                         pini, fix, u_log)
        if ((abs(rmsd1 - rmsd0) < rmsd_target) .and. (iter > 50)) then
           conv = .true.
           exit
        end if
        call simplex_step(param, rmsd1, np-nfix, rmsd_target, popt, rmsd0, conv)
        if (conv) exit
        rmsd0 = rmsd1
     end do
  end do repeat

  call simplex_reset()
  call simplex_final()

  !------------------------------ output ------------------------------!

  open(u_out, file=trim(outfile), status='replace', action='write')
  write(u_out, *) np
  ip = 0
  do i = 1, np
     if (fix(i)) then
        write(u_out, *) pini(i), 1
     else
        ip = ip + 1
        write(u_out, *) param(ip), 0
     end if
  end do
  write(u_out, *)
  write(u_out, '(1x,70("-"))')
  write(u_out, *)
  write(*,*)
  if (conv) then
     write(u_out, '(1x,"Converged after ",I8," iterations.")') iter
     write(*,     '(1x,"Converged after ",I8," iterations.")') iter
  else
     write(u_out, '(1x,"NOT converged after ",I8," iterations.")') iter
     write(*,     '(1x,"NOT converged after ",I8," iterations.")') iter
  end if
  write(u_out, '(1x,"Final RMSD: ", F15.6)') rmsd1
  write(*,     '(1x,"Final RMSD (check ",A," for details): ", F15.6)') &
        trim(logfile), rmsd1
  write(*,*)
  write(u_out, *)
  ip = 0
  do i = 1, 2
     if (fix(i)) then
        write(u_out, '(1x,F6.3)', advance='no') pini(i)
     else
        ip = ip + 1
        write(u_out, '(1x,F6.3)', advance='no') param(ip)
     end if
  end do
  do i = 3, np
     if (fix(i)) then
        write(u_out, '(1x,ES14.6)', advance='no') pini(i)
     else
        ip = ip + 1
        write(u_out, '(1x,ES14.6)', advance='no') param(ip)
     end if
  end do
  write(u_out, *)
  write(u_out, *)
  close(u_out)

  call plot_fit(ftype, np, nfix, param, pini, fix)

  !--------------------------- finalization ---------------------------!

  write(u_log, '("#",A12,2x,A12,2x,A12)') 'RMSD (sum)  ', 'RMSD (data) ', &
                                          'RMSD (deriv)'
  close(u_log)
  deallocate(pini, fix, param, popt)
  deallocate(xdata, ydata)
  if (usederiv) deallocate(xderiv, yderiv)

contains

  function get_rmsd(ndata, xdata, ydata, usederiv, nderiv, xderiv, yderiv, &
                    ftype, np, nf, pnow, pini, fix, u_log) result(rmsd)

    implicit none

    integer,                             intent(in) :: ndata
    double precision, dimension(ndata),  intent(in) :: xdata, ydata
    logical,                             intent(in) :: usederiv
    integer,                             intent(in) :: nderiv
    double precision, dimension(nderiv), intent(in) :: xderiv, yderiv
    character(len=*),                    intent(in) :: ftype
    integer,                             intent(in) :: np, nf
    double precision, dimension(np-nf),  intent(in) :: pnow
    double precision, dimension(np),     intent(in) :: pini
    logical,          dimension(np),     intent(in) :: fix
    integer,                             intent(in) :: u_log
    double precision                                :: rmsd

    double precision, dimension(np) :: param
    double precision                :: rmsd_data, rmsd_deriv
    double precision                :: x, y
    double precision                :: f, df, d2f, d3f
    integer                         :: i, ip

!    double precision                :: k, w

    ip = 0
    do i = 1, np
       if (fix(i)) then
          param(i) = pini(i)
       else
          ip = ip + 1
          param(i) = pnow(ip)
       end if
    end do

!    k = log(25.0d0)/(param(2) - param(1))

    rmsd_data = 0.0d0
    do i = 1, ndata
       x = xdata(i)
       y = ydata(i)
       if (x < param(1)) cycle
       select case(trim(ftype))
       case('polyexp')
          call eval_polyexp(x, np, param, f, df, d2f, d3f)
       case default
          stop 999
       end select
!       w = exp(-k*x)
!       rmsd_data = rmsd_data + w*(y - f)*(y - f)
       rmsd_data = rmsd_data + (y - f)*(y - f)
    end do
    rmsd_data = sqrt(rmsd_data)

    rmsd_deriv = 0.0d0
    if (usederiv) then
       do i = 1, nderiv
          x = xderiv(i)
          y = yderiv(i)
          if (x < param(1)) cycle
          select case(trim(ftype))
          case('polyexp')
             call eval_polyexp(x, np, param, f, df, d2f, d3f)
          case default
             write(0,*) 'Error: unknown function type: ', trim(ftype)
             stop
          end select
!          w = exp(-k*x)
!          rmsd_deriv = rmsd_deriv + w*(y - df)*(y - df)
          rmsd_deriv = rmsd_deriv + (y - df)*(y - df)
       end do
       rmsd_deriv = sqrt(rmsd_deriv)
    end if

!!! FIXME: looks like the rmsd is never divided by the number of points.

    rmsd = rmsd_data + rmsd_deriv

    write(u_log, '(1x,3(ES12.6,2x))') rmsd, rmsd_data, rmsd_deriv

  end function get_rmsd

!!$  !--------------------------------------------------------------------!
!!$  !                        "minimum margin" fit                        !
!!$  !--------------------------------------------------------------------!
!!$
!!$  function get_margin(ndata, xdata, ydata, usederiv, nderiv, xderiv, yderiv, &
!!$                    ftype, np, nf, pnow, pini, fix, u_log) result(M)
!!$
!!$    implicit none
!!$
!!$    integer,                             intent(in) :: ndata
!!$    double precision, dimension(ndata),  intent(in) :: xdata, ydata
!!$    logical,                             intent(in) :: usederiv
!!$    integer,                             intent(in) :: nderiv
!!$    double precision, dimension(nderiv), intent(in) :: xderiv, yderiv
!!$    character(len=*),                    intent(in) :: ftype
!!$    integer,                             intent(in) :: np, nf
!!$    double precision, dimension(np-nf),  intent(in) :: pnow
!!$    double precision, dimension(np),     intent(in) :: pini
!!$    logical,          dimension(np),     intent(in) :: fix
!!$    integer,                             intent(in) :: u_log
!!$    double precision                                :: M
!!$
!!$    double precision, dimension(np) :: param
!!$    double precision                :: x, y
!!$    double precision                :: f, df, d2f, d3f
!!$    integer                         :: i, ip
!!$
!!$    double precision                :: k, w
!!$
!!$    ip = 0
!!$    do i = 1, np
!!$       if (fix(i)) then
!!$          param(i) = pini(i)
!!$       else
!!$          ip = ip + 1
!!$          param(i) = pnow(ip)
!!$       end if
!!$    end do
!!$
!!$
!!$  end function get_margin

  !--------------------------------------------------------------------!

  subroutine plot_fit(ftype, np, nf, pnow, pini, fix)

    implicit none

    character(len=*),                   intent(in) :: ftype
    integer,                            intent(in) :: np, nf
    double precision, dimension(np-nf), intent(in) :: pnow
    double precision, dimension(np),    intent(in) :: pini
    logical,          dimension(np),    intent(in) :: fix

    integer, parameter :: u_plt = 44
    double precision, dimension(np) :: param
    double precision                :: x, dx
    double precision                :: f, df, d2f, d3f
    integer                         :: i, n, ip

    open(u_plt, file=trim(plotfile), status='replace', action='write')

    ip = 0
    do i = 1, np
       if (fix(i)) then
          param(i) = pini(i)
       else
          ip = ip + 1
          param(i) = pnow(ip)
       end if
    end do

    dx = 0.1d0
    n = nint((param(2)-param(1))/dx) + 1
    x = param(1)
    
    do i = 1, n      
       select case(trim(ftype))
       case('polyexp')
          call eval_polyexp(x, np, param, f, df, d2f, d3f)
       case default
          stop 999
       end select
       write(u_plt, '(1x,F12.6,2x,4(ES15.8,1x))') x, f, df, d2f, d3f
       x = x + dx
    end do

    close(u_plt)

  end subroutine plot_fit

  !--------------------------------------------------------------------!

  subroutine eval_polyexp(x, np, param, y, dy, d2y, d3y)
    
    implicit none

    double precision,                intent(in)  :: x
    integer,                         intent(in)  :: np
    double precision, dimension(np), intent(in)  :: param
    double precision,                intent(out) :: y, dy, d2y, d3y

    double precision                :: rcut

    double precision                :: fx, dfx, d2fx, d3fx
    double precision                :: tx, dtx, d2tx, d3tx
    double precision                :: fc, dfc, d2fc, d3fc
   
    rcut = param(2)

    if (x >= rcut) then
       y  = 0.0d0
       dy = 0.0d0
    else

       ! polynomial x exponential decay --> at the cut-off radius:
       call polyexp(rcut, param(3), np-4, param(4:np), fc, dfc, d2fc, d3fc)

       ! polynomial x exponential decay --> at x:
       call polyexp(x, param(3), np-4, param(4:np), fx, dfx, d2fx, d3fx)
       ! value at x for the 2nd degree Taylor expansion around rcut:
       call polyexp_cut(x, rcut, fc, dfc, d2fc, tx, dtx, d2tx, d3tx)

       ! subtract Taylor expansion from function:
       y   = fx   - tx
       dy  = dfx  - dtx
       d2y = d2fx - d2tx
       d3y = d3fx - d3tx

    end if
    
  end subroutine eval_polyexp


end program fit
