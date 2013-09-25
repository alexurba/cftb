module simplex

  !--------------------------------------------------------------------!
  ! Multidimensional optimization using the downhill simplex method    !
  ! (Nelder and Mead, Computer Journal 7 (1965) 308).                  !
  ! Don't use this algorithm if you have gradients available!!!        !
  !                                                                    !
  ! Two implementations are provided:                                  !
  !                                                                    !
  ! (1) A fast optimization procedure `simplex_min', which expects a   !
  !     set of starting vectors along with the target function.  The   !
  !     convenience interface `simplex_optimize' can be used with just !
  !     a single initial guess.                                        !
  !                                                                    !
  ! (2) A step-wise optimization is possible with the `simplex_step'   !
  !     procedure.  The routine expects a single vector and the        !
  !     corresponding function value.  The next trial vector is        !
  !     returned for evaluation. This implementation saves the current !
  !     state of the algorithm either to the memory, or to a file,     !
  !     which might come handy if the function evaluation takes much   !
  !     more time then one  optimization step.                         !
  !                                                                    !
  !     In order to allocate memory to save the current status, the    !
  !     subroutine simplex_init() has to be called.  Deallocation is   !
  !     done in simplex_final().  None of those is necessary if a save !
  !     file is desired.                                               !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! 2010-08-04 Alexander Urban (AU)                                    !
  ! 2011-03-07 AU --- memory based step-wise optimization              !
  !--------------------------------------------------------------------!

  implicit none
  save

  public  :: simplex_optimize,  &
             simplex_reset,     &
             simplex_step,      &
             simplex_min,       &
             simplex_init,      &
             simplex_final

  private :: simplex_hilo,      &
             simplex_start,     &
             simplex_try

  integer,          parameter, private :: u_save   = 20
  character(len=*), parameter, private :: savefile = 'SIMPLEX.save'

  logical,                                       private :: memopt = .false.

  integer,                                       private :: sav_ndim
  integer,                                       private :: sav_last
  integer,                                       private :: sav_ninit
  integer,                                       private :: sav_ihi, sav_inhi, sav_ilo
  
  double precision, dimension(:,:), allocatable, private :: sav_p
  double precision, dimension(:),   allocatable, private :: sav_psum
  double precision, dimension(:),   allocatable, private :: sav_y

contains


  !--------------------------------------------------------------------!
  !    user interface for optimizations with the SIMPLEX algorithm     !
  !--------------------------------------------------------------------!

  subroutine simplex_optimize(f, x, yopt, ndim, iter, ftol)

    implicit none
    
    !------------------------------------------------------------------!
    ! ndim : number of dimensions of the optimization problem          !
    ! x(i) : on entry: i-th parameter of the initial guess point       !
    !        on exit : optimized solution vector                       !
    ! yopt : function value for the optimized x vector                 !
    ! iter : on entry: max. number of iterations                       !
    !        on exit : actual number of iterations needed              !
    ! ftol : convergence criterion (max. change of the function value) !
    !------------------------------------------------------------------!

    integer,                           intent(in)    :: ndim
    double precision, dimension(ndim), intent(inout) :: x
    double precision,                  intent(out)   :: yopt
    integer,                           intent(inout) :: iter
    double precision,                  intent(in)    :: ftol

    interface
       function f(x,n) result(y)
         double precision, dimension(n), intent(in) :: x
         integer,                        intent(in) :: n
         double precision                           :: y
       end function f
    end interface

    integer                                       :: itmax
    double precision, dimension(:,:), allocatable :: p
    double precision, dimension(:),   allocatable :: y

    itmax = iter

    allocate(p(ndim+1,ndim), y(ndim+1))

    call simplex_start(x, ndim, p, y, f, lambda=1.0d-1)
    call simplex_min(p, y, ndim, ftol, f, iter)

    ! always restart SIMPLEX one time:
    if (iter < itmax) then
       itmax = iter
       x(1:ndim) = p(1, 1:ndim)
       call simplex_start(x, ndim, p, y, f, lambda=1.0d-3)
       call simplex_min(p, y, ndim, ftol, f, iter)   
    end if
    iter = iter + itmax

    ! save best solution found:
    x(1:ndim) = p(1, 1:ndim)
    yopt      = y(1)

    deallocate(p, y)

  end subroutine simplex_optimize

  !--------------------------------------------------------------------!
  !                   initialization of the simplex                    !
  !                                                                    !
  ! (ndim + 1) trial vectors are constructed in proximity to the       !
  ! initial guess vector x0.                                           !
  ! The function values for the trial vectors are stored in y.         !
  !--------------------------------------------------------------------!

  subroutine simplex_start(x0, ndim, p, y, f, lambda)

    implicit none
    
    integer,                                  intent(in)  :: ndim
    double precision, dimension(ndim),        intent(in)  :: x0
    double precision, dimension(ndim+1,ndim), intent(out) :: p
    double precision, dimension(ndim+1),      intent(out) :: y
    double precision,                         intent(in)  :: lambda

    interface
       function f(x,n) result(y)
         double precision, dimension(n), intent(in) :: x
         integer,                        intent(in) :: n
         double precision                           :: y
       end function f
    end interface

    integer                     :: itrial

    p(1,1:ndim) = x0(1:ndim)
    y(1)        = f(x0,ndim)
    do itrial = 2, ndim + 1
       p(itrial,1:ndim)   = x0(1:ndim)
       p(itrial,itrial-1) = p(itrial,itrial-1)+lambda*x0(itrial-1)
       y(itrial)          = f(p(itrial,1:ndim),ndim)
    end do

  end subroutine simplex_start

  !--------------------------------------------------------------------!
  !                 Simplex algorithm (implementation)                 !
  !                                                                    !
  !      (cp. (C) Copr. 1986-92 Numerical Recipes Software ,4-#)       !
  !--------------------------------------------------------------------!

  subroutine simplex_min(p, y, ndim, ftol, func, iter)

    implicit none

    !------------------------------------------------------------------!
    ! ndim   : number of dimensions of the optimization problem        !
    ! p(m,n) : on entry: n-th component of the m-th trial vector       !
    !          on exit : same, but new (ndim+1) new points that are    !
    !                    within the convergence criterion              !
    ! y(m)   : y(n) = func(p(m,:))                                     !
    ! ftol   : convergence criterion                                   !
    ! iter   : on entry: max. number of iterations                     !
    !          on exit : actual number of iterations performed         !
    !------------------------------------------------------------------!

    integer,                                  intent(in)    :: ndim
    double precision, dimension(ndim+1,ndim), intent(inout) :: p
    double precision, dimension(ndim+1),      intent(inout) :: y
    double precision,                         intent(in)    :: ftol
    integer,                                  intent(inout) :: iter

    interface
       function func(x,n) result(y)
         double precision, dimension(n), intent(in) :: x
         integer,                        intent(in) :: n
         double precision                           :: y
       end function func
    end interface

    integer                           :: itmax
    integer                           :: i, ihi, ilo, inhi, j, n
    double precision                  :: rtol, swap, ysave, ytry
    double precision, dimension(ndim) :: psum

    itmax = iter
    iter  = 0

    ! sum up same components of the trials
    ! --> redo this after contraction of the simplex
    do n = 1, ndim
       psum(n) = sum(p(:,n))
    end do

    opt : do

       ! determine best (lowest), worst (highest)
       ! and 2nd worst (next highest) trials:
       ilo = 1
       if (y(1) > y(2)) then
          ihi  = 1
          inhi = 2
       else
          ihi  = 2
          inhi = 1
       end if
       do i = 1, ndim + 1
          if(y(i) <= y(ilo)) ilo = i
          if(y(i) >  y(ihi)) then
             inhi = ihi
             ihi  = i
          else if(y(i) > y(inhi)) then
             if(i /= ihi) inhi = i
          end if
       end do

       ! range from highest to lowest trial:
       rtol = 2.0d0*abs(y(ihi) - y(ilo))/( abs(y(ihi)) + abs(y(ilo)) )

       ! check, whether convergence criterion has been met, or the
       ! max. number of iterations has been exceeded:
       if ( (rtol < ftol) .or. (iter >= itmax)) then
          ! on exit: p(1) and y(1) contain the best trial:
          swap   = y(1)
          y(1)   = y(ilo)
          y(ilo) = swap
          do n = 1, ndim
             swap     = p(1,n)
             p(1,n)   = p(ilo,n)
             p(ilo,n) = swap
          end do
          exit opt
       end if

       iter = iter + 2

       ! reflect simplex from the high point:
       ! (= extrapolation by factor -1 through the face of the simplex)
       ytry = simplex_try(p, y, psum, ndim, func, ihi, -1.0d0)
    
       if (ytry <= y(ilo)) then
          ! result is better then former low point: 
          ! --> try additional extrapolation by factor 2.0 (expansion)
          ytry = simplex_try(p, y, psum, ndim, func, ihi, 2.0d0)
       else if (ytry >= y(inhi)) then
          ! result is worse then the former second highest point:
          ! --> try something else
          ysave = y(ihi)

          ! 1-D contraction by factor 0.5:
          ytry  = simplex_try(p, y, psum, ndim, func, ihi, 0.5d0)

          if (ytry >= ysave) then
             ! still not better:
             ! --> contract around the lowest (best) point:
             do i = 1, ndim + 1
                if (i /= ilo) then
                   do j = 1, ndim
                      psum(j) = 0.5*(p(i,j) + p(ilo,j))
                      p(i,j)  = psum(j)
                   end do
                   y(i) = func(psum,ndim)
                end if
             end do
             iter = iter + ndim

             ! update psum:
             do n = 1, ndim
                psum(n) = sum(p(:,n))
             end do
          end if

       else
          iter = iter - 1
       end if

    end do opt

  end subroutine simplex_min

  !--------------------------------------------------------------------!
  !                  trial simplex optimization step                   !
  !--------------------------------------------------------------------!

  function simplex_try(p, y, psum, ndim, func, ihi, fac) result(ytry)

    implicit none

    !------------------------------------------------------------------!
    ! ndim   : number of dimensions of the optimization problem        !
    ! p(m,n) : on entry: n-th component of the m-th trial vector       !
    !          on exit : same, but updated if the optimization trial   !
    !                    step was successful                           !
    ! psum(n) = sum(p(:,n))                                            !
    ! y(m)   : y(n) = func(p(m,:))                                     !
    ! ihi    : number of the highest (worst) trial point               !
    ! fac    : factor for contraction or expansion, e.g.               !
    !          fac = -1.0 --> reflection (w/o contraction/expansion)   !
    !                 0.5 --> contraction by 1/2                       !
    !                 2.0 --> expansion by 2                           !
    !------------------------------------------------------------------!

    integer,                                  intent(in)    :: ndim, ihi
    double precision, dimension(ndim+1,ndim), intent(inout) :: p
    double precision, dimension(ndim),        intent(inout) :: psum
    double precision, dimension(ndim+1),      intent(inout) :: y
    double precision,                         intent(in)    :: fac
    double precision                                        :: ytry

    interface
       function func(x,n) result(y)
         double precision, dimension(n), intent(in) :: x
         integer,                        intent(in) :: n
         double precision                           :: y
       end function func
    end interface

    integer                           :: j
    double precision                  :: fac1, fac2
    double precision, dimension(ndim) :: ptry

    fac1 = (1.0d0 - fac)/ndim
    fac2 = fac1 - fac
    do j = 1, ndim
       ptry(j) = psum(j)*fac1 - p(ihi,j)*fac2
    end do
    ytry = func(ptry,ndim)
    if (ytry < y(ihi)) then
       y(ihi) = ytry
       do j = 1, ndim
          psum(j)  = psum(j) - p(ihi,j) + ptry(j)
          p(ihi,j) = ptry(j)
       end do
    end if

  end function simplex_try


  !====================================================================!



  !--------------------------------------------------------------------!
  !         allocate memory for optimization without save file         !
  !--------------------------------------------------------------------!

  subroutine simplex_init(ndim)

    implicit none

    integer, intent(in) :: ndim

    allocate(sav_p(ndim+1, ndim), sav_psum(ndim), sav_y(ndim+1))

    sav_ndim    = ndim
    sav_last    = 0
    sav_ninit   = 0
    sav_p(:,:)  = 0.0d0
    sav_psum(:) = 0.0d0
    sav_y(:)    = 0.0d0
    sav_ihi     = 1
    sav_inhi    = 1
    sav_ilo     = 1
    
    memopt = .true.

  end subroutine simplex_init

  !--------------------------------------------------------------------!

  subroutine simplex_final()

    implicit none

    if (memopt) then
       deallocate(sav_p, sav_psum, sav_y)
       memopt = .false.
    end if

  end subroutine simplex_final

  !--------------------------------------------------------------------!

  subroutine simplex_save(ndim, ninit, p, psum, y, ihi, inhi, ilo, last)

    implicit none


    integer,                                  intent(in) :: ndim
    integer,                                  intent(in) :: ninit
    double precision, dimension(ndim+1,ndim), intent(in) :: p
    double precision, dimension(ndim),        intent(in) :: psum
    double precision, dimension(ndim+1),      intent(in) :: y
    integer,                                  intent(in) :: ihi, ilo, inhi
    integer,                                  intent(in) :: last

    if (memopt) then
       if (ndim /= sav_ndim) then
          write(0,*) "Error: something wrong in `simplex_restore()'"
          call simplex_final()
          stop
       end if
       sav_ninit   = ninit  
       sav_p(:,:)  = p(:,:) 
       sav_psum(:) = psum(:)
       sav_y(:)    = y(:)   
       sav_ihi     = ihi    
       sav_inhi    = inhi   
       sav_ilo     = ilo    
       sav_last    = last   
    else
       open(u_save, file=savefile, status='replace', form='unformatted', action='write')
       write(u_save) ndim
       write(u_save) ninit
       write(u_save) p(1:ndim+1,1:ndim), psum(1:ndim), y(1:ndim+1)
       write(u_save) ihi, inhi, ilo, last
       close(u_save)
    end if

  end subroutine simplex_save

  !--------------------------------------------------------------------!

  subroutine simplex_restore(ndim, ninit, p, psum, y, ihi, inhi, ilo, last)

    implicit none

    integer,                                  intent(in)  :: ndim
    integer,                                  intent(out) :: ninit
    double precision, dimension(ndim+1,ndim), intent(out) :: p
    double precision, dimension(ndim),        intent(out) :: psum
    double precision, dimension(ndim+1),      intent(out) :: y
    integer,                                  intent(out) :: ihi, ilo, inhi
    integer,                                  intent(out) :: last

    integer :: ndim2
    logical :: fexists

    if (memopt) then
       if (ndim /= sav_ndim) then
          write(0,*) "Error: something wrong in `simplex_restore()'"
          call simplex_final()
          stop
       end if
       ninit   = sav_ninit
       p(:,:)  = sav_p(:,:)
       psum(:) = sav_psum(:)
       y(:)    = sav_y(:)
       ihi     = sav_ihi
       inhi    = sav_inhi
       ilo     = sav_ilo
       last    = sav_last
    else
       inquire(file=savefile, exist=fexists)
       if (fexists) then
          open(u_save, file=savefile, status='old', form='unformatted', action='read')
          read(u_save) ndim2
          if (ndim2 /= ndim) then
             write(0,*) 'Error: incompatible save file (SIMPLEX).'
             close(u_save)
             call simplex_final()
             stop
          end if
          read(u_save) ninit
          read(u_save) p(1:ndim+1,1:ndim), psum(1:ndim), y(1:ndim+1)
          read(u_save) ihi, inhi, ilo, last
          close(u_save)
       else
          ninit   = 0
          p(:,:)  = 0.0d0
          psum(:) = 0.0d0
          y(:)    = 0.0d0
          ihi     = 1
          inhi    = 1
          ilo     = 1
          last    = 0
       end if
    end if

  end subroutine simplex_restore


  !====================================================================!


  !--------------------------------------------------------------------!
  !            user interface for single optimization steps            !
  !--------------------------------------------------------------------!

  subroutine simplex_step(x1, y1, ndim, ftol, xlow, ylow, conv)
  
    implicit none

    !------------------------------------------------------------------!
    ! ndim : number of dimensions of the optimization problem          !
    ! x(i) : on entry: i-th parameter of the initial guess point       !
    !        on exit : optimized solution vector                       !
    ! yopt : function value for the optimized x vector                 !
    ! iter : on entry: max. number of iterations                       !
    !        on exit : actual number of iterations needed              !
    ! ftol : convergence criterion (max. change of the function value) !
    !------------------------------------------------------------------!

    integer,                           intent(in)    :: ndim
    double precision, dimension(ndim), intent(inout) :: x1
    double precision,                  intent(in)    :: y1
    double precision,                  intent(in)    :: ftol
    double precision, dimension(ndim), intent(out)   :: xlow
    double precision,                  intent(out)   :: ylow
    logical,                           intent(out)   :: conv
    
    integer   :: itrial, idim, i, j


    double precision, dimension(ndim)             :: xnew
    double precision, dimension(:,:), allocatable :: p
    double precision, dimension(:),   allocatable :: psum
    double precision, dimension(:),   allocatable :: y
    integer                                       :: ninit
    integer                                       :: ihi, inhi, ilo
    ! laststep = 0 --> initialization
    !            1 --> reflection
    !            2 --> expansion
    !            3 --> contraction
    integer                                       :: laststep, currentstep

    double precision, parameter :: lambda = 1.0d-1
    double precision            :: fac, fac1, fac2

    conv = .false.

    allocate(p(ndim+1,ndim), psum(ndim), y(ndim+1))
    
    call simplex_restore(ndim, ninit, p, psum, y, ihi, inhi, ilo, laststep)

    if (ninit > 0) then

       ! the following action depends on the last step:
       select case(laststep)

       case(0) ! initialization
          ninit    = ninit + 1
          y(ninit) = y1
          if (ninit < (ndim+1)) then
             ! continue to initialize trials:
             currentstep = 0
             xnew = p(ninit+1, 1:ndim)
          else if (ninit == (ndim+1)) then
             ! all trials initialized
             do idim = 1, ndim
                psum(idim) = sum(p(:,idim))
             end do
             ! 
             call simplex_hilo(ndim, y, ftol, ihi, inhi, ilo, conv)
             ! try reflection:
             currentstep = 1
             fac  = -1.0d0
             fac1 = (1.0d0 - fac)/dble(ndim)
             fac2 = fac1 - fac
             xnew(1:ndim) = psum(1:ndim)*fac1 - p(ihi,1:ndim)*fac2
          end if

       case(1) ! reflection
          if (y1 <= y(ilo)) then
             ! accept last step:
             psum(1:ndim)   = psum(1:ndim) - p(ihi,1:ndim)
             p(ihi, 1:ndim) = x1(1:ndim)
             psum(1:ndim)   = psum(1:ndim) + p(ihi,1:ndim)
             y(ihi)         = y1
             ! try expansion:
             currentstep = 2
             fac  = 2.0d0
             fac1 = (1.0d0 - fac)/dble(ndim)
             fac2 = fac1 - fac
             xnew(1:ndim) = psum(1:ndim)*fac1 - p(ihi,1:ndim)*fac2
          else if (y1 >= y(inhi)) then
             if (y1 < y(ihi)) then
                ! accept last step:
                psum(1:ndim)   = psum(1:ndim) - p(ihi,1:ndim)
                p(ihi, 1:ndim) = x1(1:ndim)
                psum(1:ndim)   = psum(1:ndim) + p(ihi,1:ndim)
                y(ihi)         = y1
             end if
             ! try contraction:
             currentstep = 3
             fac  = 0.5d0
             fac1 = (1.0d0 - fac)/dble(ndim)
             fac2 = fac1 - fac
             xnew(1:ndim) = psum(1:ndim)*fac1 - p(ihi,1:ndim)*fac2
          else
             ! accept last step:
             psum(1:ndim)   = psum(1:ndim) - p(ihi,1:ndim)
             p(ihi, 1:ndim) = x1(1:ndim)
             psum(1:ndim)   = psum(1:ndim) + p(ihi,1:ndim)
             y(ihi)         = y1
             call simplex_hilo(ndim, y, ftol, ihi, inhi, ilo, conv)
             ! try reflection:
             currentstep = 1
             fac  = -1.0d0
             fac1 = (1.0d0 - fac)/dble(ndim)
             fac2 = fac1 - fac
             xnew(1:ndim) = psum(1:ndim)*fac1 - p(ihi,1:ndim)*fac2
          end if

       case(2) ! expansion
          if (y1 < y(ihi)) then
             ! accept last step:
             psum(1:ndim)   = psum(1:ndim) - p(ihi,1:ndim)
             p(ihi, 1:ndim) = x1(1:ndim)
             psum(1:ndim)   = psum(1:ndim) + p(ihi,1:ndim)
             y(ihi)         = y1
             call simplex_hilo(ndim, y, ftol, ihi, inhi, ilo, conv)
          end if
          ! try reflection:
          currentstep = 1
          fac  = -1.0d0
          fac1 = (1.0d0 - fac)/dble(ndim)
          fac2 = fac1 - fac
          xnew(1:ndim) = psum(1:ndim)*fac1 - p(ihi,1:ndim)*fac2

       case(3) ! contraction
          if (y1 < y(ihi)) then
             ! accept last step:
             psum(1:ndim)   = psum(1:ndim) - p(ihi,1:ndim)
             p(ihi, 1:ndim) = x1(1:ndim)
             psum(1:ndim)   = psum(1:ndim) + p(ihi,1:ndim)
             y(ihi)         = y1
             call simplex_hilo(ndim, y, ftol, ihi, inhi, ilo, conv)
             ! try reflection:
             currentstep = 1
             fac  = -1.0d0
             fac1 = (1.0d0 - fac)/dble(ndim)
             fac2 = fac1 - fac
             xnew(1:ndim) = psum(1:ndim)*fac1 - p(ihi,1:ndim)*fac2
          else
             ! contract around the lowest (best) point:
             do i = 1, ndim + 1
                if (i /= ilo) then
                   do j = 1, ndim
                      psum(j) = 0.5*(p(i,j) + p(ilo,j))
                      p(i,j)  = psum(j)
                   end do
                end if
             end do
             ! re-initialize:
             xnew(1:ndim)  = p(ilo,1:ndim)
             p(ilo,1:ndim) = p(1,1:ndim)
             p(1,1:ndim)   = xnew(1:ndim)
             y(1)          = y(ilo)
             ninit         = 1
             currentstep   = 0
             xnew(1:ndim)  = p(2,1:ndim)
          end if

       case default
          write(0,*) 'Error: you should never come here (SIMPLEX)!'
          deallocate(p,psum,y)
          stop

       end select

    else 

       ! no restore information found --> start initialization
       p(1,1:ndim) = x1(1:ndim)
       y(1)        = y1
       do itrial = 2, ndim + 1
          p(itrial,1:ndim)   = x1(1:ndim)
          p(itrial,itrial-1) = p(itrial,itrial-1) + lambda*x1(itrial-1)
          y(itrial)          = 0.0d0
       end do
       ninit    = 1 ! one corner of the simplex
       ihi      = 0
       inhi     = 0
       ilo      = 0
       laststep = 0 ! 0 = initialization
       ! return next trial structure for evaluation:
       xnew(1:ndim) = p(2,1:ndim)
       currentstep  = 0

    end if

    ! save current status:
    call simplex_save(ndim, ninit, p, psum, y, ihi, inhi, ilo, currentstep)

    ! after the initialization always do the following:
    if (currentstep /= 0) then
       xlow(1:ndim) = p(ilo,1:ndim)
       ylow         = y(ilo)
    end if

    ! override input x vector with new x vector:
    x1(:) = xnew(:)

    deallocate(p, psum, y)

  end subroutine simplex_step

  !--------------------------------------------------------------------!
  !                          reset save file                           !
  !--------------------------------------------------------------------!
  
  subroutine simplex_reset()

    implicit none

    if (memopt) then
       sav_last    = 0
       sav_ninit   = 0
       sav_p(:,:)  = 0.0d0
       sav_psum(:) = 0.0d0
       sav_y(:)    = 0.0d0
       sav_ihi     = 1
       sav_inhi    = 1
       sav_ilo     = 1
    else
       open(u_save, file=trim(savefile), status='replace')
       close(u_save, status='delete')
    end if

  end subroutine simplex_reset

  !--------------------------------------------------------------------!
  !          determine highest, 2nd highest and lowest point           !
  !--------------------------------------------------------------------!

  subroutine simplex_hilo(ndim, y, ftol, ihi, inhi, ilo, conv)

    implicit none

    integer,                             intent(in)  :: ndim
    double precision, dimension(ndim+1), intent(in)  :: y
    double precision,                    intent(in)  :: ftol
    integer,                             intent(out) :: ihi, inhi, ilo
    logical,                             intent(out) :: conv

    double precision                                 :: rtol
    integer          :: i
   
    ilo = 1
    if (y(1) > y(2)) then
       ihi  = 1
       inhi = 2
    else
       ihi  = 2
       inhi = 1
    end if
    do i = 1, ndim + 1
       if(y(i) <= y(ilo)) ilo = i
       if(y(i) >  y(ihi)) then
          inhi = ihi
          ihi  = i
       else if(y(i) > y(inhi)) then
          if(i /= ihi) inhi = i
       end if
    end do
    
    ! check for convergence:
!    rtol = 2.0d0*abs(y(ihi) - y(ilo))/( abs(y(ihi)) + abs(y(ilo)) )
    rtol = abs(y(ihi) - y(ilo))
    if (rtol <= ftol) then
       conv = .true.
    else
       conv = .false.
    end if

  end subroutine simplex_hilo

end module simplex

