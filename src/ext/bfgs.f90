module bfgs

  !--------------------------------------------------------------------!
  ! Function minimization using the BFGS quasi-Newton method.          !
  !--------------------------------------------------------------------!
  ! The code is losely based on Bernd Meyer's `relax' program.         !
  !--------------------------------------------------------------------!
  ! 2011-04-08 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  implicit none

  integer,          parameter, private :: u_save   = 40
  character(len=*), parameter, private :: savefile = 'BFGS.save'

  !--------------------------------------------------------------------!

  double precision, dimension(3), public  :: bfgs_a0   = 1.0d-1
  double precision,               public  :: bfgs_up   = 5.0d-5
  logical,                        public  :: bfgs_conv = .false.

  !--------------------------------------------------------------------!

  double precision, dimension(:),   allocatable, private :: bfgs_X, bfgs_X_old
  double precision, dimension(:),   allocatable, private :: bfgs_F, bfgs_F_old
  double precision, dimension(:,:), allocatable, private :: bfgs_A
  integer,                                       private :: bfgs_DoF
  logical,                                       private :: bfgs_restart

  logical, private :: isInit = .false.
  
contains

  subroutine bfgs_init(DoF, alpha, update, file)
    
    implicit none

    integer,                    intent(in) :: DoF
    double precision, optional, intent(in) :: alpha
    double precision, optional, intent(in) :: update
    logical,          optional, intent(in) :: file

    if (isInit) then
       write(0,*) "Warning: redundant initialization (module `bfgs')."
       return
    end if

    if (present(alpha)) then
       bfgs_a0(:) = alpha
    end if

    if (present(update)) then
       bfgs_up = alpha
    end if

    if ((.not. present(file)) .or. (.not. file)) then

       bfgs_DoF     = DoF
       bfgs_restart = .false.

       allocate(bfgs_X(DoF), bfgs_X_old(DoF), &
                bfgs_F(DoF), bfgs_F_old(DoF), &
                bfgs_A(DoF,DoF))
    
       isInit = .true.

    end if

  end subroutine bfgs_init

  !--------------------------------------------------------------------!
  
  subroutine bfgs_final()

    implicit none
    
    if (isInit) then
       deallocate(bfgs_X, bfgs_X_old, bfgs_F, bfgs_F_old, bfgs_A)
       isInit = .false.
    end if

  end subroutine bfgs_final

  !--------------------------------------------------------------------!
  !                    restart information handling                    !
  !--------------------------------------------------------------------!

  subroutine bfgs_save(n, X, F, A)

    implicit none

    integer,                          intent(in) :: n
    double precision, dimension(n),   intent(in) :: X
    double precision, dimension(n),   intent(in) :: F
    double precision, dimension(n,n), intent(in) :: A

    if (isInit)  then
       bfgs_DoF     = n
       bfgs_X(:)    = X(:)
       bfgs_F(:)    = F(:)
       bfgs_A(:,:)  = A(:,:)
       bfgs_restart = .true.
    else
       open(u_save, file=savefile, status='replace', action='write', &
                    form='unformatted')
       write(u_save) n
       write(u_save) X(:)
       write(u_save) F(:)
       write(u_save) A(:,:)
       close(u_save)
    end if

  end subroutine bfgs_save

  !--------------------------------------------------------------------!

  subroutine bfgs_restore(n, X, F, A, restart)

    implicit none

    integer,                          intent(in)  :: n
    double precision, dimension(n),   intent(out) :: X
    double precision, dimension(n),   intent(out) :: F
    double precision, dimension(n,n), intent(out) :: A
    logical,                          intent(out) :: restart

    logical :: fexists
    integer :: n1, i

    locale : if (isInit) then
       ! load from memory:
       if (n /= bfgs_DoF) then
          write(0,*) "Error: incompatible restart information" &
                  // " in `bfgs_restore()'."
          stop
       else
          if (bfgs_restart) then
             X(:)    = bfgs_X(:)
             F(:)    = bfgs_F(:)
             A(:,:)  = bfgs_A(:,:)
             restart = .true.
          else
             X(:) = 0.0d0
             F(:) = 0.0d0
             ! Initialize inverse Jacobian matrix:
             ! FIXME: generalize this to arbitrary dimensions.
             A(1:n,1:n) = 0.D0
             do i = 1, n
                A(i,i) = bfgs_a0(mod(i-1,3)+1)
             end do
             restart = .false.
          end if
       end if
    else 
       inquire(file=savefile, exist=fexists)
       if (.not. fexists) then
          X(:) = 0.0d0
          F(:) = 0.0d0
          ! Initialize inverse Jacobian matrix:
          ! FIXME: generalize this to arbitrary dimensions.
          A(1:n,1:n) = 0.D0
          do i = 1, n
             A(i,i) = bfgs_a0(mod(i-1,3)+1)
          end do
          restart = .false.
       else
          ! load from file:
          open(u_save, file=savefile, status='old', action='read', &
                       form='unformatted')
          read(u_save) n1
          if (n /= n1) then
             write(0,*) "Error: incompatible restart information" &
                     // " in `bfgs_restore()'."
             stop
          else
             read(u_save) X(:)   
             read(u_save) F(:)   
             read(u_save) A(:,:) 
          end if
          close(u_save)
          restart = .true.
       end if
    end if locale

  end subroutine bfgs_restore

  !--------------------------------------------------------------------!
  !                implementation of the BFGS algorithm                !
  !--------------------------------------------------------------------!

  subroutine bfgs_optimize(n, X, F, tol, F_max, F_rms, conv)

    implicit none

    integer,                        intent(in)    :: n
    double precision, dimension(n), intent(inout) :: X
    double precision, dimension(n), intent(in)    :: F
    double precision,               intent(in)    :: tol
    double precision,               intent(out)   :: F_max
    double precision,               intent(out)   :: F_rms
    logical,                        intent(out)   :: conv

    double precision, dimension(n)   :: X_old, F_old
    double precision, dimension(n,n) :: A
    logical                          :: restart

    call bfgs_fstatus(n, F, F_max, F_rms)
    if (abs(F_max) < tol) then
       conv = .true.
       return
    else
       conv = .false.
    end if

    call bfgs_restore(n, X_old, F_old, A, restart)

    if (restart .and. (F_rms > bfgs_up)) then
       call bfgs_jacobian(n, X, X_old, F, F_old, A)
    end if

    call bfgs_save(n, X, F, A)

    call bfgs_parameters(n, X, F, A)

  end subroutine bfgs_optimize

  !--------------------------------------------------------------------!
  !                   update inverse Jacobian matrix                   !
  !--------------------------------------------------------------------!

  subroutine bfgs_jacobian(n, X, X_old, F, F_old, A)

    implicit none
    
    integer,                          intent(in)    :: n
    double precision, dimension(n),   intent(in)    :: X, X_old
    double precision, dimension(n),   intent(in)    :: F, F_old
    double precision, dimension(n,n), intent(inout) :: A

    double precision              :: fac, F_norm, F_scal
    double precision, dimension(n):: av, va
    double precision, dimension(n):: dX
    double precision, dimension(n):: dF
    integer                       :: i, j
    
    do i = 1, n
       dX(i) = X(i) - X_old(i)
       dF(i) = F_old(i) - F(i)
    end do

    ! av(:) = matmul(A,dF) ?
    av(1:n) = 0.0d0
    do j = 1, n
       do i = 1, n
          av(i) = av(i) + A(i,j)*dF(j)
       end do
    end do

    F_scal = 0.0d0
    F_norm = 0.0d0
    do i = 1, n
       F_scal = F_scal + dX(i)*dF(i)
       F_norm = F_norm + dF(i)*av(i)
    end do
    F_scal = 1.0d0/F_scal
    fac    = 1.0d0/F_norm

    do i = 1, n
       va(i) = F_scal*dX(i) - fac*av(i)
    end do

    do j = 1, n
       do i = 1, n
          A(i,j) = A(i,j) &
               + F_scal*dX(i)*dX(j) - fac*av(i)*av(j) + F_norm*va(i)*va(j)
       end do
    end do

  end subroutine bfgs_jacobian

  !--------------------------------------------------------------------!
  !               update parameters / optimization step                !
  !--------------------------------------------------------------------!

  subroutine bfgs_parameters(n, X, F, A)

    implicit none

    integer,                          intent(in)    :: n
    double precision, dimension(n),   intent(inout) :: X
    double precision, dimension(n),   intent(in)    :: F
    double precision, dimension(n,n), intent(in)    :: A

    integer          :: i, j
    double precision :: shft

    do i = 1, n
       shft = 0.0d0
       do j = 1, n
          shft = shft + A(i,j)*F(j)
       end do
       X(i) = X(i) + shft
    end do

  end subroutine bfgs_parameters


  !====================================================================!
  !                                                                    !
  !                        auxiliary routines                          !
  !                                                                    !
  !====================================================================!

  

  !--------------------------------------------------------------------!
  !        maximum and root mean squared values of derivatives F       !
  !--------------------------------------------------------------------!

  subroutine bfgs_fstatus(n, F, F_max, F_rms)
    
    implicit none

    integer,                        intent(in)  :: n
    double precision, dimension(n), intent(in)  :: F
    double precision,               intent(out) :: F_max
    double precision,               intent(out) :: F_rms

    integer          :: i

    F_max = 0.0d0
    F_rms = 0.0d0
    do i = 1, n
       if (abs(F(i)) > abs(F_max)) F_max = F(i)
       F_rms = F_rms + F(i)*F(i)
    end do
    F_rms = sqrt(F_rms)

  end subroutine bfgs_fstatus


end module bfgs
