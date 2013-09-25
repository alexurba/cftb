program potfit

  !--------------------------------------------------------------------!
  ! Fit a pairpotential to reproduce cohesive energies.                !
  !--------------------------------------------------------------------!
  ! 2011-02-18 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use bondint,   only: polyexp, polyexp_cut

  use constants, only: Ang2Bohr

  use interpol,  only: deriv2,              &
                       cs,                  &
                       cs_d1,               &
                       cs_d2,               &
                       cs_load,             &
                       cs_simplify

  use io,        only: io_adjustl,          &
                       io_replace,          &
                       io_trim

  use simplex,   only: simplex_step,        &
                       simplex_reset,       &
                       simplex_init,        &
                       simplex_final

  use tbio,      only: tbio_init,           &
                       tbio_final,          &
                       tbio_args,           &
                       tbio_print_usage,    &
                       args_switch,         &
                       args_get,            &
                       nargs, wdescr

  implicit none


  integer,          parameter :: u_out = 40
  integer,          parameter :: u_plt = 41
  integer,          parameter :: u_trn = 42
  integer,          parameter :: u_tst = 43

  !----------------------- boundary conditions ------------------------!
  ! rmin, rmax     : interval for the potential function               !
  !                  function and derivative will go to zero at r=rmax !
  ! y1_1           : first derivative at r=rmin (degree of freedom)    !
  ! y1_n           : first derivative at r=rmax --> y1_n = 0.0d0       !
  !                                                                    !
  !---------------------------- data files ----------------------------!
  ! nfiles         : total number of data files                        !
  ! ntrain, ntest  : number of training/testing data files             !
  ! files(i)       : file name of the i-th file                        !
  ! type1, type2   : atom types for the pair potential                 !
  ! iE, nE         : counter for all structures in memory              !
  ! nE_train/_test : number of structures in training/testing set      !
  !                                                                    !
  !------------------------------ bonds -------------------------------!
  ! disthist       : histogram of bond lengths                         !
  ! nhist          : number of bins for the histogram                  !
  !                                                                    !
  !------------------------ energy evaluation -------------------------!
  ! From all data files the relevant distances between rmin and rmax   !
  ! and the selected atom types are stored in memory.                  !
  ! ndist_max      : max. number of different distances per data file  !
  ! ndist_inc      : increment of ndist_max, if more memory has to be  !
  !                  allocated                                         !
  ! ndist(i)       : number of different distances of file i           !
  ! p_dist(j,i)    : value of the j-th distance of file i              !
  ! p_count(j,i)   : count of the j-th distance of file i              !
  ! p_dist and p_count are dynamically increased if the allocated      !
  ! memory is insufficient.                                            !
  !                                                                    !
  !------------------------- fitting function -------------------------!
  ! ftype          : function type                                     !
  ! nDOF           : number of degrees of freedom                      !
  ! bias           : constant to be added to the pair potential sum    !
  !                                                                    !
  ! (1) cubic splines:                                                 !
  ! ------------------                                                 !
  ! nnodes         : number of abscissas for the spline function       !
  ! x(i), y(i)     : x/y values at the abscissa i                      !
  ! y2(i)          : d^2 y/dx^2 at the abscissas (cubic splines)       !
  !                                                                    !
  ! (2) polyexp functions:                                             !
  ! ----------------------                                             !
  ! np           : number of parameters                              !
  ! nfixed        : number of fixed parameters; nfix < n              !
  ! pini(i)    : i-th parameter                                    !
  ! fixed(i)      " .true., if the i-th parameter is fixed            !
  !                                                                    !
  !------------------------------- fit --------------------------------!
  ! E_target(i)    : target energy for the pair potential for file i   !
  ! E_now(i)       : during fit: energy calculated with current guess  !
  ! E_tb(i)        : TB-Band cohesive energy (without pot.) for file i !
  ! E_cohes(i)     : DFT cohesive energy for file i                    !
  ! rmsd_E         : rmsd(energy)                                      !
  ! rmsd_dE        : rmsd(d energy / d structure)                      !
  ! dmean_E        : mean deviation from the target energies           !
  ! dmax_E         : max. deviation from the target energy             !
  ! rmsd_err       : rmsd(mean deviation)                              !
  ! rmsd0, rmsd1   : previous and current RMSD for the energy          !
  ! param(i)       : i-th optimization parameter                       !
  !                  parameters are: y(1:nnodes-1), y1_1               !
  ! popt(i)        : after fit: optimal value for the i-th parameter   !
  ! conv           : .true. if the fit converged                       !
  !                                                                    !
  !----------------------------- options ------------------------------!
  ! infile         : name of the input file with initial guess         !
  ! outfile        : output file name                                  !
  ! niter          : max. number of iterations                         !
  ! tol            : optimization threshold                            !
  ! opt_target     : set the target function for the optimization      !
  !--------------------------------------------------------------------!

  double precision                                         :: rmin
  double precision                                         :: rmax
  double precision                                         :: y1_1
  double precision                                         :: y1_n

  integer                                                  :: nfiles
  integer                                                  :: itrain, ntrain
  integer                                                  :: itest, ntest
  integer                                                  :: iE, nE
  integer                                                  :: nE_test, nE_train
  character(len=100), dimension(:),   allocatable          :: files
  character(len=2)                                         :: type1, type2

  double precision,   dimension(:),   allocatable          :: disthist
  integer                                                  :: nhist

  integer                                                  :: ndist_max, ndist_inc
  integer,            dimension(:),   allocatable          :: ndist

  double precision,   dimension(:,:),              pointer :: p_dist
  integer,            dimension(:,:),              pointer :: p_count
  double precision,   dimension(:,:), allocatable, target  :: t1_dist, t2_dist
  integer,            dimension(:,:), allocatable, target  :: t1_count, t2_count

  character(len=50)                                        :: ftype
  integer                                                  :: nDOF
  double precision                                         :: bias

  integer                                                  :: nnodes
  double precision,   dimension(:),   allocatable          :: x, y, y2
  double precision                                         :: a

  integer                                                  :: np, nfixed
  double precision,   dimension(:),   allocatable          :: pini
  logical,            dimension(:),   allocatable          :: fixed

  double precision,   dimension(:),   allocatable          :: E_target, E_now, E_tb, E_cohes
  double precision                                         :: rmsd0, rmsd1
  double precision                                         :: rmsd_E, rmsd_dE
  double precision                                         :: dmax_E, dmean_E
  double precision                                         :: rmsd_err
  double precision,   dimension(:),   allocatable          :: param, popt
  logical                                                  :: conv

  character(len=100)                                       :: infile, outfile
  character(len=100)                                       :: plotfile, restfile
  character(len=20)                                        :: nodepos
  integer                                                  :: iter, niter, irep, nrep
  double precision                                         :: tol

  character(len=50)                                        :: opt_target
  logical                                                  :: optderiv

  character(len=100) :: str
  integer            :: i, ipos

  !----------------------------- options ------------------------------!

  call parse_command_arguments()


  !-------------------------- initialization --------------------------!  

  nrep = 5
  y1_n = 0.0d0

  ! retrieve fitting information from data files:
  ndist_max = 30
  ndist_inc = 30
  allocate(t1_dist(ndist_max, nfiles),  &
           t1_count(ndist_max, nfiles), &
           ndist(nfiles),               &
           E_target(nfiles),            &
           E_now(nfiles),               &
           E_tb(nfiles),                &
           E_cohes(nfiles))
  t1_count(:,:) = 0
  p_dist  => t1_dist
  p_count => t1_count
  ndist(:)    = 0
  E_target(:) = 0.0d0
  E_now(:)    = 0.0d0
  E_tb(:)     = 0.0d0
  E_cohes(:)  = 0.0d0

  select case(ftype)
  case('csplines')
     call init_csplines_fit()
  case('polyexp', 'sumexp', 'gsp')
     call init_fit()
  case default
     write(0,*) 'Error: unknown function type: ', trim(ftype)
     call finalize()
     stop
  end select

  call parse_data_files()

  nhist = 30
  allocate(disthist(nhist))
  call analyze_distances(nE_train, ndist_max, ndist, p_dist, p_count, &
                         rmin, rmax, nhist, disthist)


  !----------------------------- fitting ------------------------------!

  call simplex_init(nDOF)

  bias = 0.0d0
  iter = 0
  repeat : do irep = 1, nrep
     call simplex_reset()
     conv = .false.
     rmsd0 = 0.0d0
     do i = 1, (niter - iter)
        iter  = iter + 1

        ! evaluate the RMSD:
        select case(ftype)
        case('csplines')
           call deriv2(x, y, nnodes, y1_1, y1_n, y2)
           call eval_RMSD_csplines(x, y, y2, nnodes, nE_train,         &
                ndist_max, ndist, p_dist, p_count, E_tb, E_cohes,      &
                E_target, E_now, rmsd_E, rmsd_dE, dmax_E, dmean_E,     &
                rmsd_err)

           param(1:nnodes-1) = y(1:nnodes-1)
           if (optderiv) param(nnodes) = y1_1
           if (trim(nodepos) == 'opt') then
              param(nnodes+1:2*nnodes-2) = x(2:nnodes-1)
           end if

        case default
           call eval_RMSD(nDOF, param, nE_train, ndist_max,    &
                ndist, p_dist, p_count, E_tb, E_cohes, E_target,       &
                E_now, rmsd_E, rmsd_dE, dmax_E, dmean_E, rmsd_err)
        end select

        ! select target function for optimization:
        select case(opt_target)
        case('absolute','abs')
           rmsd1 = rmsd_E + rmsd_dE
        case('relative','rel')  
           rmsd1 = rmsd_err + rmsd_dE
           bias  = dmean_E
        case('differences','diff')  
           rmsd1 = rmsd_dE
           bias  = dmean_E
        end select
        if ((abs(rmsd1 - rmsd0) < tol) .and. (iter > 50)) then
           conv = .true.
           exit
        end if

        ! SIMPLEX optimization step:
        call simplex_step(param, rmsd1, nDOF, tol, popt, rmsd0, conv)

        if (trim(ftype) == 'csplines') then

           y(1:nnodes-1) = param(1:nnodes-1)
           if (optderiv) y1_1 = param(nnodes)
           if (trim(nodepos) == 'opt') then
              x(2:nnodes-1) = param(nnodes+1:2*nnodes-2)
           end if

           ! auto-derivative for the cubic splines fit:
           if (.not. optderiv) then
              ! the derivative at the first node is chosen to resemble an
              ! exponential decay:
              if ((y(1) < 0.0d0) .or. (y(2) < 0.0d0)) exit
              a    = (log(y(2)) - x(2)/x(1)*log(y(1)))/(x(1) - x(2)) - log(y(1))/x(1)
              y1_1 = -a*y(1)
              if (y1_1 > 0.0d0) then
                 y(1) = y(2)
                 y1_1 = 0.0d0
              end if
           end if

        end if

        ! check for convergence:
        if (conv) then
           select case(ftype)
           case('csplines')
              y(1:nnodes-1) = popt(1:nnodes-1)
              if (optderiv) y1_1 = popt(nnodes)
              if (trim(nodepos) == 'opt') then
                 x(2:nnodes-1) = popt(nnodes+1:2*nnodes-2)
              end if
              call deriv2(x, y, nnodes, y1_1, y1_n, y2)
              call eval_RMSD_csplines(x, y, y2, nnodes, nE_train,      &
                   ndist_max, ndist, p_dist, p_count, E_tb, E_cohes,   &
                   E_target, E_now, rmsd_E, rmsd_dE, dmax_E, dmean_E,  &
                   rmsd_err)
              exit
           case default
              call eval_RMSD(nDOF, param, nE_train, ndist_max, &
                   ndist, p_dist, p_count, E_tb, E_cohes, E_target,    &
                   E_now, rmsd_E, rmsd_dE, dmax_E, dmean_E, rmsd_err)
           end select
        end if

        rmsd0 = rmsd1
     end do
  end do repeat

  call simplex_reset()
  call simplex_final()

  ! simple evaluation of a restarted fit: no optimization
  if (iter == 0) then
     select case(ftype)
     case('csplines')
        call deriv2(x, y, nnodes, y1_1, y1_n, y2)
        call eval_RMSD_csplines(x, y, y2, nnodes, nE_train, ndist_max, &
             ndist, p_dist, p_count, E_tb, E_cohes, E_target, E_now,   &
             rmsd_E, rmsd_dE, dmax_E, dmean_E, rmsd_err)
     case default
        call eval_RMSD(nDOF, param, nE_train, ndist_max, &
             ndist, p_dist, p_count, E_tb, E_cohes, E_target,    &
             E_now, rmsd_E, rmsd_dE, dmax_E, dmean_E, rmsd_err)
     end select
  end if

  ! energy correction shift:
  select case (opt_target)
  case ('relative','rel') 
     bias  = dmean_E
  case('differences','diff')
     bias  = dmean_E
  end select


  !------------------------------ output ------------------------------!

  ! write to stdout and the output file:
  call print_results()

  ! write out datapoints for plotting:
  select case(ftype)
  case('csplines')
     call print_plot_csplines()
  case default
     call print_plot(nDOF, param)
  end select

  ! training structures energies:
  write(u_trn,'("#",8x,4(A15,2x))') 'E_tb+E_pot', 'E_train', 'E_tb', 'E_pot'
  open(u_trn, file='train.dat', status='replace', action='write')
  do iE = 1, nE_train
     if (ndist(iE) == 0) cycle
     write(u_trn,'(3x,I4,2x,4(ES15.8,2x))') iE, E_tb(iE)+E_now(iE)+bias, &
          E_cohes(iE), E_tb(iE), E_now(iE)
  end do
  close(u_trn)

  ! testing structures energies:
  if (ntest > 0) then
     i = nE_train + 1

     select case(ftype)
     case('csplines')
        call eval_RMSD_csplines(x, y, y2, nnodes, nE_test, ndist_max,   &
             ndist(i:nE), p_dist(:,i:nE), p_count(:,i:nE), E_tb(i:nE),  &
             E_cohes(i:nE), E_target(i:nE), E_now(i:nE), rmsd_E,        &
             rmsd_dE, dmax_E, dmean_E, rmsd_err)
     case default
        call eval_RMSD(nDOF, param, nE_test, ndist_max,                 &
             ndist(i:nE), p_dist(:,i:nE), p_count(:,i:nE), E_tb(i:nE),  &
             E_cohes(i:nE), E_target(i:nE), E_now(i:nE), rmsd_E,        &
             rmsd_dE, dmax_E, dmean_E, rmsd_err)
     end select

     write(u_tst,'("#",8x,4(A15,2x))') 'E_tb+E_pot', 'E_test', 'E_tb', 'E_pot'
     open(u_tst, file='test.dat', status='replace', action='write')
     do iE = nE_train + 1, nE
        if (ndist(iE) == 0) cycle
        write(u_tst,'(3x,I4,2x,4(ES15.8,2x))') iE, E_tb(iE)+E_now(iE)+bias, &
             E_cohes(iE), E_tb(iE), E_now(iE)
     end do
     close(u_tst)
     write(*,*)
  end if

  !--------------------------- finalization ---------------------------!

  call finalize()


contains !=============================================================!


  !--------------------------------------------------------------------!
  !              finalize everything / deallocate memory               !
  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none
    
    if (allocated(ndist)) then
       deallocate(ndist, E_target, E_now, E_tb, E_cohes)
    end if

    if (allocated(disthist)) deallocate(disthist)
    if (allocated(files))    deallocate(files)
    if (allocated(x))        deallocate(x, y, y2)
    if (allocated(param))    deallocate(param, popt)
    if (allocated(pini)) deallocate(pini, fixed)

    if (allocated(t1_dist)) then
       deallocate(t1_dist, t1_count)
    else
       deallocate(t2_dist, t2_count)
    end if
    
  end subroutine finalize

  !--------------------------------------------------------------------!
  !                      command line parameters                       !
  !--------------------------------------------------------------------!

  subroutine parse_command_arguments()

    implicit none

    call tbio_init(                                                     &
         name='potfit.x',                                               &
         descr='fit the repulsive pair potential',                      &
         args=(/ 'in     ', 'out    ', 'niter  ', 'tol    ', 'files  ', &
                 'funct  ', 'plot   ', 'nnodes ', 'restart'/),          &
         extra='--bond:--range:--nodepos:--test:--opttarget:--optderiv')
  
    call tbio_args('in',      infile)
    call tbio_args('restart', restfile)
    call tbio_args('out',     outfile,  default='potfit.out')
    call tbio_args('plot',    plotfile, default='potfit.dat')
    call tbio_args('niter',   niter)
    call tbio_args('tol',     tol,      default=1.0d-10)
    call tbio_args('nnodes',  nnodes)
    call tbio_args('funct',   ftype,    default='csplines')
  
    ! potential range / interval:
    call args_switch('--range', pos=ipos)
    if (ipos /= 0) then
       call args_get(ipos+1, str)
       read(str,*) rmin
       call args_get(ipos+2, str)
       read(str,*) rmax
       write(*,*) io_trim('Potential range [Bohr]', wdescr), ' : [', &
            trim(io_adjustl(rmin,2)), ', ', trim(io_adjustl(rmax)) // ']'
    else if (len_trim(restfile) == 0) then
       write(*,*)
       write(*,*) "Please specify min./max. distance for the potential."
       write(*,*) "(you can also use the `--range' switch for this purpose)"
       write(*,*)
       write(*,'(1x,"min. r : ")', advance='no')
       read(*,*) rmin
       write(*,'(1x,"max. r : ")', advance='no')
       read(*,*) rmax
       write(*,*)
    end if
  
    ! atom types / bond:
    call args_switch('--bond', pos=ipos)
    if (ipos /= 0) then
       call args_get(ipos+1, type1)
       call args_get(ipos+2, type2)
       write(*,*) io_trim('Atom types', wdescr), ' : ', &
            trim(adjustl(type1)), ' ', trim(adjustl(type2))
    else
       write(*,*)
       write(*,*) "Please specify atom types."
       write(*,*) "(you can also use the `--bond' switch for this purpose)"
       write(*,*)
       write(*,'(1x,"type 1: ")', advance='no')
       read(*,*) type1
       write(*,'(1x,"type 2: ")', advance='no')
       read(*,*) type2
       write(*,*)
    end if
  
    ! training set:
    call args_switch('--train', pos=itrain)
    if (itrain == 0) then
       write(0,*) 'Error: No training data files.'
       call tbio_print_usage()
       call tbio_final()
       stop
    end if
    ntrain = 0
    do i = itrain + 1, nargs
       call args_get(i, str)
       if (str(1:1) == '-') exit
       ntrain = ntrain + 1
    end do
  
    ! testing set:
    call args_switch('--test', pos=itest)
    ntest = 0
    if (itest > 0) then
       do i = itest + 1, nargs
          call args_get(i, str)
          if (str(1:1) == '-') exit
          ntest = ntest + 1
       end do
    end if
  
    nfiles = ntrain + ntest
    allocate(files(nfiles))
    write(*,*) io_trim('Number of data files', wdescr), ' : ', &
               trim(io_adjustl(nfiles))
    write(*,*)
    write(*,*) 'Training set:'
    write(*,*)
    do i = 1, ntrain
       call args_get(i+itrain, files(i))
       write(*,'(3x,I4,2x)', advance='no') i
       write(*,*) trim(files(i))
    end do
    if (ntest > 0) then
       write(*,*)
       write(*,*) 'Testing set (only for comparison):'
       write(*,*)
       do i = 1, ntest
          call args_get(i+itest, files(ntrain+i))
          write(*,'(3x,I4,2x)', advance='no') i
          write(*,*) trim(files(ntrain+i))
       end do
    end if
  
    !--------------------------------------------!
    ! additional (for now undocumented) switches !
    !--------------------------------------------!
   
    call args_switch('--nodepos',   value=nodepos,    default='uniform')
    call args_switch('--opttarget', value=opt_target, default='relative')
    call args_switch('--optderiv',  value=optderiv)
  
    call tbio_final()

  end subroutine parse_command_arguments

  !--------------------------------------------------------------------!
  !            print results to stdout and the output file             !
  !--------------------------------------------------------------------!

  subroutine print_results()
    
    implicit none

    if (conv) then
       write(*,*) 'Converged after ', trim(io_adjustl(iter)), ' iterations.'
    else
       write(*,*) 'NOT Converged after ', trim(io_adjustl(iter)), ' iterations.'
    end if
  
    write(*,*) 
    write(*,*) 'Fitting results:'
    write(*,*) 
    write(*,'(1x,"target function      : ",A15)') opt_target
    write(*,'(1x,"fitting function     : ",A15)') ftype
    write(*,*)
    write(*,'(1x,"rmsd(E)              : ",ES15.8)') rmsd_E
    write(*,'(1x,"rmsd(dE/dx)          : ",ES15.8)') rmsd_dE
    write(*,'(1x,"max. absolute error  : ",ES15.8)') dmax_E
    write(*,'(1x,"mean error           : ",ES15.8)') dmean_E
    write(*,'(1x,"rmsd(mean error)     : ",ES15.8)') rmsd_err
    write(*,'(1x,"(bias                : ",ES15.8,")")') bias
    write(*,*)
  
    select case(ftype)
    case('csplines')
       call print_results_csplines()
    case default
       call print_final_parameters(nDOF, param, conv, iter)
    end select

  end subroutine print_results

  !--------------------------------------------------------------------!
  !                        parse the data files                        !
  !--------------------------------------------------------------------!

  subroutine parse_data_files()

    use pbc,       only: pbc_d_max_origin,    &
                         pbc_number_of_tvecs, &
                         pbc_compute_tvecs

    use xsflib,    only: xsf_init,  &
                         xsf_final, &
                         avec,      &
                         coorat,    &
                         forlat,    &
                         nameat,    &
                         natoms,    &
                         ntype,     &
                         E_coh, E_tbband

    implicit none

    integer                                       :: nallatoms
    double precision, dimension(:,:), allocatable :: coord
    integer,          dimension(:),   allocatable :: atomtype

    double precision                              :: dmax
    integer,          dimension(:,:), allocatable :: tvec
    integer                                       :: ntvecs
    double precision, dimension(:),   allocatable :: tveclen

    double precision, dimension(3) :: vec
    integer                        :: ifile, itvec, itype1, itype2, itype
    integer                        :: i, iatom1, iatom2
    logical                        :: twotypes, newdist
    integer                        :: i1, i2
    double precision               :: dist2, rmax2, r

    double precision, parameter    :: dtol = 1.0d-6

    rmax2 = rmax*rmax

    nE       = 0
    nE_test  = 0
    nE_train = 0
    datafiles : do ifile = 1, nfiles

       call xsf_init(files(ifile))
       itype1 = -1
       itype2 = -1
       do itype = 1, ntype
          if (trim(nameat(itype)) == trim(type1)) itype1 = itype
          if (trim(nameat(itype)) == trim(type2)) itype2 = itype
       end do

       ! only, if the relevant atom types are present:
       process : if ((itype1 > 0) .and. (itype2 > 0)) then

          nE = nE + 1
          if (ifile <= ntrain) then
             nE_train = nE_train + 1
          else
             nE_test = nE_test + 1
          end if

          ! target energy (Ry to Hartree conversion):
          E_tb(nE)     = E_tbband
          E_cohes(nE)  = 0.5d0*E_coh
          E_target(nE) = 0.5d0*E_coh - E_tbband

          ! Angstroem to Bohr:
          avec(:,:) = avec(:,:)*Ang2Bohr

          ! rearrange coordinates in memory; keep only atoms 
          ! of required types:
          nallatoms = natoms(itype1)
          if (itype1 /= itype2) then
             nallatoms = nallatoms + natoms(itype2)
             twotypes = .true.
          else
             twotypes = .false.
          end if
          allocate(coord(3,nallatoms), atomtype(nallatoms))

          i1 = 1
          i2 = natoms(itype1)
          coord(:,i1:i2)  = reshape(coorat(:,:,itype1), shape(coord(:,i1:i2)))
          atomtype(i1:i2) = itype1
          if (twotypes) then
             i1 = natoms(itype1) + 1
             i2 = nallatoms
             coord(:,i1:i2)  = reshape(coorat(:,:,itype2), shape(coord(:,i1:i2)))
             atomtype(i1:i2) = itype2
          end if
          coord(:,:) = coord(:,:)*Ang2Bohr

          ! translation vectors:
          dmax   = pbc_d_max_origin(nallatoms, coord, avec)
          ntvecs = pbc_number_of_tvecs(rmax, dmax, avec)
          allocate(tvec(3,ntvecs), tveclen(ntvecs))
          call pbc_compute_tvecs(rmax, dmax, avec, ntvecs, tvec, tveclen)

          atom1 : do iatom1 = 1, nallatoms
             itype1 = atomType(iatom1)

             ! atoms in the unit cell T = (0,0,0)
             do iatom2 = iatom1+1, nallatoms
                itype2 = atomType(iatom2)
                if (twotypes .and. (itype1 == itype2)) cycle
                vec(1:3) = coord(1:3,iatom2) - coord(1:3,iatom1)
                dist2 = sum(vec*vec)
                if (dist2 > rmax2) cycle
                r = sqrt(dist2)
                newdist = .true.
                do i = 1, ndist(nE)
                   if (abs(r - p_dist(i,nE)) < dtol) then
                      p_dist(i,nE) = ( p_dist(i,nE)*dble(p_count(i,nE)) &
                                   +   r)/dble(p_count(i,nE) + 1)
                      p_count(i,nE) = p_count(i,nE) + 1
                      newdist = .false.
                   end if
                end do
                if (newdist) then
                   if (ndist(nE) >= ndist_max) call more_memory()
                   ndist(nE) = ndist(nE) + 1 
                   p_dist(ndist(nE), nE)  =  r
                   p_count(ndist(nE), nE) =  1
               end if

             end do

             ! periodic pictures of the unit cell; T /= (0,0,0):
             translations : do itvec = 1, ntvecs
                if (tveclen(itvec) - 2.0d0*dmax > rmax) exit translations
                do iatom2 = 1, nallatoms
                   itype2 = atomtype(iatom2)
                   if (twotypes .and. (itype1 == itype2)) cycle
                   vec(1:3) = dble(tvec(1,itvec))*avec(:,1) &
                            + dble(tvec(2,itvec))*avec(:,2) &
                            + dble(tvec(3,itvec))*avec(:,3)
                   vec(1:3) = vec(1:3) + (coord(1:3,iatom2) - coord(1:3,iatom1))
                   dist2 = sum(vec*vec)
                   if (dist2 > rmax2) cycle
                   r = sqrt(dist2)
                   newdist = .true.
                   do i = 1, ndist(nE)
                      if (abs(r - p_dist(i,nE)) < dtol) then
                         p_dist(i,nE) = ( p_dist(i,nE)*dble(p_count(i,nE)) &
                                      +   r)/dble(p_count(i,nE) + 1)
                         p_count(i,nE) = p_count(i,nE) + 1
                         newdist = .false.
                      end if
                   end do
                   if (newdist) then
                      if (ndist(nE) >= ndist_max) call more_memory()
                      ndist(nE) = ndist(nE) + 1
                      p_dist(ndist(nE), nE)  =  r
                      p_count(ndist(nE), nE) =  1
                   end if

                end do
             end do translations
          end do atom1

          deallocate(coord, atomtype, tvec, tveclen)
       end if process

       call xsf_final()

    end do datafiles

  end subroutine parse_data_files

  !--------------------------------------------------------------------!
  !                  enlarge allocated memory region                   !
  !--------------------------------------------------------------------!

  subroutine more_memory()

    implicit none

    if (allocated(t1_dist)) then
       allocate(t2_dist(ndist_max+ndist_inc,  nfiles),  &
                t2_count(ndist_max+ndist_inc, nfiles)   )
       t2_dist(1:ndist_max,:)  = t1_dist(1:ndist_max,:)
       t2_count(1:ndist_max,:) = t1_count(1:ndist_max,:)
       t2_count(ndist_max+1:ndist_max+ndist_inc,:) = 0
       p_dist  => t2_dist
       p_count => t2_count
       deallocate(t1_dist, t1_count)
    else
       allocate(t1_dist(ndist_max+ndist_inc,  nfiles),  &
                t1_count(ndist_max+ndist_inc, nfiles)   )
       t1_dist(1:ndist_max,:)  = t2_dist(1:ndist_max,:)
       t1_count(1:ndist_max,:) = t2_count(1:ndist_max,:)
       t1_count(ndist_max+1:ndist_max+ndist_inc,:) = 0
       p_dist  => t1_dist
       p_count => t1_count
       deallocate(t2_dist, t2_count)
    end if

    ndist_max = ndist_max + ndist_inc

  end subroutine more_memory

  !--------------------------------------------------------------------!
  !                       distances distribution                       !
  !--------------------------------------------------------------------!

  subroutine analyze_distances(nE, nd, ndist, dist, count, r0, r1, nr, p)

    implicit none
    
    integer,                             intent(in)  :: nE, nd
    integer,          dimension(nE),     intent(in)  :: ndist
    double precision, dimension(nd, nE), intent(in)  :: dist
    integer,          dimension(nd, nE), intent(in)  :: count
    double precision,                    intent(in)  :: r0, r1
    integer,                             intent(in)  :: nr
    double precision, dimension(nr),     intent(out) :: p
    
    integer, parameter :: u_hist = 50

    double precision :: r, dr
    integer          :: ir, iE, id

    p(:) = 0

    dr = (r1-r0)/dble(nr)
    r  = r0
    do ir = 1, nr
       do iE = 1, nE
       do id = 1, ndist(iE)
          if ((dist(id,iE) > r) .and. (dist(id,iE) <= r+dr)) then
             p(ir) = p(ir) + count(id,iE)
          end if
       end do
       end do
       r = r + dr
    end do

    p(:) = p(:)/sum(p(:))/dr

    open(u_hist, file='disthist.dat', status='replace', action='write')
    r  = r0 + 0.5d0*dr
    do ir = 1, nr
       write(u_hist,'(1x,2(ES15.8,2x))') r, p(ir)
       r = r + dr
    end do
    close(u_hist)

  end subroutine analyze_distances




  !====================================================================!
  !                                                                    !
  !                       cubic splines fitting                        !
  !                                                                    !
  !====================================================================!




  !--------------------------------------------------------------------!
  !                initialize/restart cubic splines fit                !
  !--------------------------------------------------------------------!

  subroutine init_csplines_fit()
    
    implicit none

    ! restarting previous fit?
    if (len_trim(restfile) > 0) then
       nnodes = 0
       call cs_load(trim(restfile), x, y, y1_1, y1_n, nnodes)
       allocate(x(nnodes), y(nnodes), y2(nnodes))
       call cs_load(trim(restfile), x, y, y1_1, y1_n, nnodes)
       rmin = x(1)
       rmax = x(nnodes)
       if (nodepos == 'opt') then
          nDOF = 2*nnodes - 2
       else
          nDOF = nnodes
       end if
       allocate(param(nDOF), popt(nDOF))
    else
       ! nodes/abscissas placing:
       nDOF = 0
       allocate(x(nnodes), y(nnodes), y2(nnodes))
       select case(nodepos)
       case('auto')
          nDOF = nnodes
          call uniform_nodes(rmin, rmax, nnodes, x)
       case('uniform')
          nDOF = nnodes
          call uniform_nodes(rmin, rmax, nnodes, x)
       case('opt')
          nDOF = 2*nnodes - 2
          call uniform_nodes(rmin, rmax, nnodes, x)
       case default
          write(0,*) 'Error: unknown node position command: ', trim(nodepos)
          stop
       end select
       allocate(param(nDOF), popt(nDOF))
       call read_initial_guess(trim(infile), rmin, rmax, nnodes, x, y)
       ! the derivative at the first node is chosen to resemble an
       ! exponential decay:
       a    = (log(y(2)) - x(2)/x(1)*log(y(1)))/(x(1) - x(2)) - log(y(1))/x(1)
       y1_1 = -a*y(1)
    end if

    ! determine the derivative at the first abscissa:
    if (.not. optderiv) then
       nDOF = nDOF - 1
    end if

  end subroutine init_csplines_fit

  !--------------------------------------------------------------------!
  ! Initialize vector x with uniformly distributed values in the given !
  ! interval [rmin, rmax].                                             !
  !--------------------------------------------------------------------!

  subroutine uniform_nodes(rmin, rmax, n, x)

    implicit none

    double precision,               intent(in)  :: rmin, rmax
    integer,                        intent(in)  :: n
    double precision, dimension(n), intent(out) :: x

    double precision :: dx, x0
    integer          :: i

    dx   = (rmax - rmin)/dble(n-1)
    x0   = rmin
    x(1) = rmin
    do i = 2, n-1
       x0   = x0 + dx
       x(i) = x0
    end do
    x(n) = rmax

  end subroutine uniform_nodes

  !--------------------------------------------------------------------!
  ! Read those data points from an input file, that shall be used to   !
  ! construct the initial guess.                                       !
  ! The number of data points does not nedd to be equal to the number  !
  ! of fitting nodes, but the fitting interval must be covered.        !
  !--------------------------------------------------------------------!

  subroutine read_initial_guess(fname, rmin, rmax, n, x, y)

    implicit none

    character(len=*),               intent(in)  :: fname
    double precision,               intent(in)  :: rmin, rmax
    integer,                        intent(in)  :: n
    double precision, dimension(n), intent(in)  :: x
    double precision, dimension(n), intent(out) :: y

    double precision, dimension(:), allocatable :: x_tmp, y_tmp, y2_tmp
    
    integer, parameter :: u_in = 20

    character(len=256) :: line
    logical :: fexists
    integer :: nlines, ival, nvalues, stat, i

    inquire(file=trim(fname), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: File not found: ", trim(fname)
       stop
    end if

    !------------------------ read data points ------------------------!

    open(u_in, file=trim(fname), status='old', action='read')

    nlines  = 0
    nvalues = 0
    count : do
       read(u_in, '(A)', iostat=stat) line
       if (stat /= 0) exit count
       nlines = nlines + 1

       if (len_trim(line) == 0) cycle count
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle count
       nvalues = nvalues + 1
    end do count

    allocate(x_tmp(nvalues), y_tmp(nvalues), y2_tmp(nvalues))
    rewind(u_in)

    ival = 0
    process : do i = 1, nlines
       read(u_in, '(A)') line
       line = trim(adjustl(line))
       if (len_trim(line) == 0) cycle process
       if (line(1:1) == '#') cycle process

       ival = ival + 1
       read(line, *) x_tmp(ival), y_tmp(ival)
    end do process
    
    close(u_in)

    !---------------- get values at desired abscissas -----------------!

    ! check, if the desired interval is covered:
    if ((x_tmp(1) > rmin) .or. (x_tmp(nvalues) < rmax)) then
       write(0,*) "Error: Initial guess is outside of fitting interval."
       deallocate(x_tmp, y_tmp, y2_tmp)
       stop
    end if

    call deriv2(x_tmp, y_tmp, nvalues, huge(1.0d0), 0.0d0, y2_tmp)

    do i = 1, n-1
       call cs(x_tmp, y_tmp, y2_tmp, nvalues, x(i), y(i))
    end do
    y(n) = 0.0d0

    !------------------- deallocate scratch memory --------------------!

    deallocate(x_tmp, y_tmp, y2_tmp)

  end subroutine read_initial_guess

  !--------------------------------------------------------------------!
  !           calculate the RMSD of the energy (and forces)            !
  !--------------------------------------------------------------------!

  subroutine eval_RMSD_csplines(x, y, y2, n, nE, ndist_max, ndist,    &
       dist, dcount, E_tb, E_cohes, E_target, E_now, rmsd_E, rmsd_dE, &
       dmax_E, dmean_E, rmsd_err)

    implicit none

    integer,                                       intent(in)  :: n
    double precision, dimension(n),                intent(in)  :: x, y, y2
    integer,                                       intent(in)  :: nE
    integer,                                       intent(in)  :: ndist_max
    integer,          dimension(nfiles),           intent(in)  :: ndist
    double precision, dimension(ndist_max,nfiles), intent(in)  :: dist
    integer,          dimension(ndist_max,nfiles), intent(in)  :: dcount
    double precision, dimension(nfiles),           intent(in)  :: E_tb
    double precision, dimension(nfiles),           intent(in)  :: E_cohes
    double precision, dimension(nfiles),           intent(in)  :: E_target
    double precision, dimension(nfiles),           intent(out) :: E_now
    double precision,                              intent(out) :: rmsd_E
    double precision,                              intent(out) :: rmsd_dE
    double precision,                              intent(out) :: dmax_E
    double precision,                              intent(out) :: dmean_E
    double precision,                              intent(out) :: rmsd_err

    double precision, parameter :: dE_max = 0.05d0

    integer          :: iE, idist
    double precision :: E, dE1, dE2

    rmsd_E   = 0.0d0
    rmsd_dE  = 0.0d0
    dmax_E   = 0.0d0
    dmean_E  = 0.0d0

    do iE = 1, nE

       E_now(iE) = 0.0d0
       do idist = 1, ndist(iE)
          call cs(x, y, y2, n, dist(idist,iE), E)
          E_now(iE) = E_now(iE) + dcount(idist,iE)*E
       end do

       if (iE > 1) then
          dE1 = E_cohes(iE) - E_cohes(iE-1)
          if (abs(dE1) < dE_max) then
             dE2 = (E_now(iE)+E_tb(iE)) - (E_now(iE-1)+E_tb(iE-1))
             E = dE1 - dE2

             if (sign(1.0d0,dE1) /= sign(1.0d0,dE2)) then
                E = exp(E*E)
             end if

             rmsd_dE = rmsd_dE + E*E
          end if
       end if

       E = E_target(iE) - E_now(iE)
       dmean_E = dmean_E + E
       dmax_E  = max(dmax_E, abs(E))
       rmsd_E  = rmsd_E + E*E

    end do

    rmsd_E  = sqrt(rmsd_E/dble(nE))
    rmsd_dE = sqrt(rmsd_dE/dble(nE))
    dmean_E = dmean_E/dble(nE)

    rmsd_err = 0.0d0
    do iE = 1, nE
       E = (E_target(iE) - E_now(iE)) - dmean_E
       rmsd_err = rmsd_err + E*E
    end do
    rmsd_err = sqrt(rmsd_err/dble(nE))

  end subroutine eval_RMSD_csplines

  !--------------------------------------------------------------------!
  !                      data points for plotting                      !
  !--------------------------------------------------------------------!

  subroutine print_plot_csplines()
    
    implicit none

    double precision :: x0, dx, y0, y0_1, y0_2
    integer          :: i

    ! finer grained for plotting:
    open(u_plt, file=trim(plotfile), status='replace', action='write')
    x0 = rmin
    dx = (rmax - rmin)/dble(100 - 1)
    call deriv2(x, y, nnodes, y1_1, y1_n, y2)
    do i = 1, 100
       call cs(x, y, y2, nnodes, x0, y0)
       call cs_d1(x, y, y2, nnodes, x0, y0_1)
       call cs_d2(x, y2, nnodes, x0, y0_2)
       write(u_plt,'(1x,4(ES15.8,2x))') x0, y0, y0_1, y0_2
       x0 = x0 + dx
    end do
    close(u_plt)

  end subroutine print_plot_csplines

  !--------------------------------------------------------------------!
  !                          print final fit                           !
  !--------------------------------------------------------------------!
  
  subroutine print_results_csplines()

    implicit none

    open(u_out, file=trim(outfile), status='replace', action='write')

    if (conv) then
       write(u_out,'("# Converged after")', advance='no')
       write(u_out,*) trim(io_adjustl(iter)), ' iterations.'
    else
       write(u_out,'("# NOT Converged after")', advance='no')
       write(u_out,*) trim(io_adjustl(iter)), ' iterations.'
    end if

    if (.not. optderiv) then
       write(*,'(1x,"potential continues a*exp(-b*x) with a,b = ",2(ES15.8,2x))') &
            y(1)*exp(-y1_1/y(1)*x(1)), -y1_1/y(1)
       write(*,*)
    end if

    write(u_out,'("#",3(A15,2x))') io_trim(' x',15), io_trim(' y',15), io_trim(" y'", 15)
    write(u_out,'(1x,3(ES15.8,2x))') x(1), y(1), y1_1
    do i = 2, nnodes - 1
       write(u_out,'(1x,2(ES15.8,2x))') x(i), y(i)
    end do
    write(u_out,'(1x,3(ES15.8,2x))') x(nnodes), y(nnodes), y1_n
    write(u_out,*)
    write(u_out,'("# ",2(F4.2,2x))', advance='no') rmin, rmax
    write(u_out,'(2(ES15.8,2x))', advance='no') y1_1, y1_n
    do i = 1, nnodes
       write(u_out,'(2(ES15.8,2x))', advance='no') x(i), y(i)
    end do
    write(u_out,*)

    close(u_out)

  end subroutine print_results_csplines




  !====================================================================!
  !                                                                    !
  !                      general fitting function                      !
  !                                                                    !
  !====================================================================!




  !--------------------------------------------------------------------!
  !               initialize/restart fitting parameters                !
  !--------------------------------------------------------------------!

  subroutine init_fit()

    implicit none

    integer, parameter :: u_in = 20
    integer            :: ip, ifix

    rmin = 0.0d0
    rmax = 0.0d0
    nDOF = 0

    if (len_trim(restfile) > 0) then
       open(u_in, file=trim(restfile), status='old', action='read')
    else
       open(u_in, file=trim(infile), status='old', action='read')
    end if

    read(u_in, *) np
    allocate(pini(np), fixed(np))
    nfixed = 0
    do i = 1, np
       read(u_in,*) pini(i), ifix
       if (ifix == 1) then
          nfixed = nfixed + 1
          fixed(i) = .true.
       else
          fixed(i) = .false.
       end if
    end do
    close(u_in)

    nDOF = np - nfixed
    allocate(param(nDOF), popt(nDOF))

    ip = 0
    do i = 1, np
       if (.not. fixed(i)) then
          ip = ip + 1
          param(ip) = pini(i)
       end if
    end do

    rmin = pini(1)
    rmax = pini(2)

  end subroutine init_fit

  !--------------------------------------------------------------------!
  !             evaluate RMSD for given fitting function               !
  !--------------------------------------------------------------------!

  subroutine eval_RMSD(n, pnow, nE, ndist_max, ndist, dist,    &
       dcount, E_tb, E_cohes, E_target, E_now, rmsd_E, rmsd_dE,        &
       dmax_E, dmean_E, rmsd_err) 

    implicit none

    integer,                                       intent(in)  :: n
    double precision, dimension(n),                intent(in)  :: pnow
    integer,                                       intent(in)  :: nE
    integer,                                       intent(in)  :: ndist_max
    integer,          dimension(nfiles),           intent(in)  :: ndist
    double precision, dimension(ndist_max,nfiles), intent(in)  :: dist
    integer,          dimension(ndist_max,nfiles), intent(in)  :: dcount
    double precision, dimension(nfiles),           intent(in)  :: E_tb
    double precision, dimension(nfiles),           intent(in)  :: E_cohes
    double precision, dimension(nfiles),           intent(in)  :: E_target
    double precision, dimension(nfiles),           intent(out) :: E_now
    double precision,                              intent(out) :: rmsd_E
    double precision,                              intent(out) :: rmsd_dE
    double precision,                              intent(out) :: dmax_E
    double precision,                              intent(out) :: dmean_E
    double precision,                              intent(out) :: rmsd_err

    double precision, parameter       :: dE_max = 0.05d0

    double precision, dimension(np) :: param
    integer                           :: ip, i
    integer                           :: iE, idist
    double precision                  :: E, dE1, dE2
    double precision                  :: dEdr, d2Edr2, d3Edr3

    rmsd_E   = 0.0d0
    rmsd_dE  = 0.0d0
    dmax_E   = 0.0d0
    dmean_E  = 0.0d0

    ip = 0
    do i = 1, np
       if (fixed(i)) then
          param(i) = pini(i)
       else
          ip = ip + 1
          param(i) = pnow(ip)
       end if
    end do

    do iE = 1, nE

       E_now(iE) = 0.0d0
       do idist = 1, ndist(iE)
    
          call eval_fitfunct(dist(idist,iE), ftype, np, param, E, &
                             dEdr, d2Edr2, d3Edr3)

          E_now(iE) = E_now(iE) + dcount(idist,iE)*E
       end do

       if (iE > 1) then
          dE1 = E_cohes(iE) - E_cohes(iE-1)
          if (abs(dE1) < dE_max) then
             dE2 = (E_now(iE)+E_tb(iE)) - (E_now(iE-1)+E_tb(iE-1))
             E = dE1 - dE2

             if (sign(1.0d0,dE1) /= sign(1.0d0,dE2)) then
                E = exp(E*E)
             end if

             rmsd_dE = rmsd_dE + E*E
          end if
       end if

       E = E_target(iE) - E_now(iE)
       dmean_E = dmean_E + E
       dmax_E  = max(dmax_E, abs(E))
       rmsd_E  = rmsd_E + E*E

    end do

    rmsd_E  = sqrt(rmsd_E/dble(nE))
    rmsd_dE = sqrt(rmsd_dE/dble(nE))
    dmean_E = dmean_E/dble(nE)

    rmsd_err = 0.0d0
    do iE = 1, nE
       E = (E_target(iE) - E_now(iE)) - dmean_E
       rmsd_err = rmsd_err + E*E
    end do
    rmsd_err = sqrt(rmsd_err/dble(nE))

  end subroutine eval_RMSD

  !--------------------------------------------------------------------!
  !                    write out data for plotting                     !
  !--------------------------------------------------------------------!

  subroutine print_plot(n, pnow)
    
    implicit none

    integer,                        intent(in) :: n
    double precision, dimension(n), intent(in) :: pnow

    double precision                :: x0, dx
    integer                         :: i, ip

    double precision                :: rmin, rmax
    double precision                :: E, dE, d2E, d3E
    double precision, dimension(np) :: param

    integer,              parameter :: npoints = 100

    ip = 0
    do i = 1, np
       if (fixed(i)) then
          param(i) = pini(i)
       else
          ip = ip + 1
          param(i) = pnow(ip)
       end if
    end do

    rmin = param(1)
    rmax = param(2)

    open(u_plt, file=trim(plotfile), status='replace', action='write')
    x0 = rmin
    dx = (rmax - rmin)/dble(100 - 1)
    do i = 1, 100
       call eval_fitfunct(x0, ftype, np, param, E, dE, d2E, d3E)
       write(u_plt,'(1x,4(ES15.8,2x))') x0, E, dE, d2E
       x0 = x0 + dx
    end do
    close(u_plt)
    
   end subroutine print_plot

  !--------------------------------------------------------------------!
  !                          print final fit                           !
  !--------------------------------------------------------------------!
  
  subroutine print_final_parameters(n, pnow, conv, iter)

    implicit none

    integer,                        intent(in) :: n
    double precision, dimension(n), intent(in) :: pnow
    logical,                        intent(in) :: conv
    integer,                        intent(in) :: iter

    integer :: ip, i

    open(u_out, file=trim(outfile), status='replace', action='write')

    write(u_out, *) np
    ip = 0
    do i = 1, np
       if (fixed(i)) then
          write(u_out, *) pini(i), 1
       else
          ip = ip + 1
          write(u_out, *) pnow(ip), 0
       end if
    end do

    write(u_out, *)
    write(u_out, '(1x,70("-"))')
    write(u_out, *)

    if (conv) then
       write(u_out,'("# Converged after")', advance='no')
       write(u_out,*) trim(io_adjustl(iter)), ' iterations.'
    else
       write(u_out,'("# NOT Converged after")', advance='no')
       write(u_out,*) trim(io_adjustl(iter)), ' iterations.'
    end if

    write(u_out, *)
    ip = 0
    do i = 1, 2
       if (fixed(i)) then
          write(u_out, '(1x,F6.3)', advance='no') pini(i)
       else
          ip = ip + 1
          write(u_out, '(1x,F6.3)', advance='no') pnow(ip)
       end if
    end do
    do i = 3, np
       if (fixed(i)) then
          write(u_out, '(1x,ES14.6)', advance='no') pini(i)
       else
          ip = ip + 1
          write(u_out, '(1x,ES14.6)', advance='no') pnow(ip)
       end if
    end do
    write(u_out, *)
    write(u_out, *)

    close(u_out)

  end subroutine print_final_parameters




  !====================================================================!
  !                                                                    !
  !                        evaluation routines                         !
  !                                                                    !
  !====================================================================!




  subroutine eval_fitfunct(x, ftype, np, param, y, dy, d2y, d3y)
    
    implicit none

    double precision,                intent(in)  :: x
    character(len=*),                intent(in)  :: ftype
    integer,                         intent(in)  :: np
    double precision, dimension(np), intent(in)  :: param
    double precision,                intent(out) :: y, dy, d2y, d3y
    
    select case(ftype)
    case('gsp')
       call eval_GSP(x, np, param, y, dy, d2y, d3y)
    case('polyexp')
       call eval_polyexp(x, np, param, y, dy, d2y, d3y)
    case('sumexp')
       call eval_sumexp(x, np, param, y, dy, d2y, d3y)
    case default
       stop 999
    end select

  end subroutine eval_fitfunct

  !--------------------------------------------------------------------!
  !                evaluation of the polyexp functions                 !
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
       y   = 0.0d0
       dy  = 0.0d0
       d2y = 0.0d0
       d3y = 0.0d0
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

  !--------------------------------------------------------------------!
  !             evaluation of the `sumexp' fitting function            !
  !--------------------------------------------------------------------!

  subroutine eval_sumexp(x, np, param, y, dy, d2y, d3y)

    implicit none

    double precision,                intent(in)  :: x
    integer,                         intent(in)  :: np
    double precision, dimension(np), intent(in)  :: param
    double precision,                intent(out) :: y, dy, d2y, d3y

    double precision :: rmin, rcut, r0
    double precision :: a
    integer          :: i
    double precision :: c, dc, d2c, d3c

    if (mod(np-3,2) /= 0) then
       write(0,*) "Error: Function type `sumexp' expects two parameters"
       write(0,*) "       per exponential function."
       call finalize()
       stop
    end if

    rmin = param(1)
    rcut = param(2)
    r0   = param(3)

    y   = 0.0d0
    dy  = 0.0d0
    d2y = 0.0d0
    d3y = 0.0d0

    if (x < rcut) then

       do i = 4, np, 2
          a   = param(i)*exp(param(i+1)*x)
          y   = y   + a
          a   = param(i+1)*a
          dy  = dy  + a
          a   = param(i+1)*a
          d2y = d2y + a
          a   = param(i+1)*a
          d3y = d3y + a
       end do

       if (r0 < rcut) then
          call cutoff_cos(x,r0, rcut, c, dc, d2c, d3c)
          y   =  y*c
          dy  = dy*c  + y*dc
          d2y = d2y*c + y*d2c + 2.0d0*dy*dc
          d3y = d3y*c + y*d3c + 3.0d0*(d2y*dc + dy*d2c)
       end if

    end if
    
  end subroutine eval_sumexp

  !--------------------------------------------------------------------!
  !                            GSP function                            !
  !--------------------------------------------------------------------!

  subroutine eval_GSP(x, np, param, y, dy, d2y, d3y)

    implicit none

    double precision,                intent(in)  :: x
    integer,                         intent(in)  :: np
    double precision, dimension(np), intent(in)  :: param
    double precision,                intent(out) :: y, dy, d2y, d3y

    double precision :: rmin, rcut, r0cut
    double precision :: y0, r0, rc, n, nc, r
    double precision :: a, b
    double precision :: t, dt, d2t, d3t

    if (np /= 8) then
       write(0,*) "Error: Function type `GSP' expects exactly 8 parameters."
       call finalize()
       stop
    end if

    rmin  = param(1)
    rcut  = param(2)

    y   = 0.0d0
    dy  = 0.0d0
    d2y = 0.0d0
    d3y = 0.0d0

    if (x < rcut) then

       r0cut = param(3)

       y0    = param(4)
       r0    = param(5)
       rc    = param(6)
       n     = param(7)
       nc    = param(8)

       a = y0*(r0/rc)**n * exp(n*(r0/rc)**nc)
       b = -n/rc**nc

       if (x < r0cut) then
          ! use the actual GSP function:
          r   = x
          y   = a*exp(b*r**nc)
          dy  = y*b*nc*r**(nc-1.0d0)
          d2y = dy*(dy/y + 1.0d0/r*(nc-1.0d0))
          d3y = d2y*d2y/dy + dy*(dy/y - 1.0d0/r)/r*(nc-1.0d0)
       else
          ! replace GSP function with tail function:
          r = r0cut
          y   = a*exp(b*r**nc)
          dy  = y*b*nc*r**(nc-1.0d0)
          call cutoff_poly3(x, r0cut, rcut, y, dy, t, dt, d2t, d3t)
          y   = t
          dy  = dt
          d2y = d2t
          d3y = d3t
       end if
       
    end if

  end subroutine eval_GSP

  !--------------------------------------------------------------------!
  !                      cosine cut-off function                       !
  !--------------------------------------------------------------------!

  subroutine cutoff_cos(x, r0, rc, y, dy, d2y, d3y)

    implicit none

    double precision, intent(in)  :: x, r0, rc
    double precision, intent(out) :: y, dy, d2y, d3y

    double precision,   parameter :: pi = 3.14159265358979d0

    double precision :: a, a2

    if (x<r0) then
       y   = 1.0d0
       dy  = 0.0d0
       d2y = 0.0d0
       d3y = 0.0d0
    else if (x>rc) then
       y   = 0.0d0
       dy  = 0.0d0
       d2y = 0.0d0
       d3y = 0.0d0
    else
       a   = pi/(rc-r0)
       a2  = a*a
       y   =  0.5d0*(cos((x-r0)*a) + 1.0d0)
       dy  =  0.5d0*sin((x-r0)*a)*a
       d2y = -(y - 0.5d0)*a2
       d3y = -dy*a2
    end if

  end subroutine cutoff_cos

  !--------------------------------------------------------------------!
  !                      polynomial tail cut-off                       !
  !--------------------------------------------------------------------!

  subroutine cutoff_poly3(r, r1, rc, f1, df1, t, dt, d2t, d3t)

    implicit none

    !------------------------------------------------------------------!
    ! The 3rd order polynomial tail function t(r) takes a function     !
    ! f(r) smoothly to zero.  f(r) will be replaced by the polynomial  !
    ! on the interval [r1,rc].                                         !
    !                                                                  !
    ! t(r) = A*(rc - r)**2 + B*(rc - r)**3                             !
    ! with:  A =  [ 3*f(r1) - df(r1)/dr * (r1-rc) ]/(rc-r1)**2         !
    !        B = -[ 2*f(r1) + df(r1)/dr * (rc-r1) ]/(rc-r1)**3         !
    !                                                                  !
    ! r   : evaulate the function at r                                 !
    ! rc  : cut-off radius                                             !
    ! r1  : t(r) is defined on [r1,rc]                                 !
    ! f1  : f(r1)                                                      !
    ! df1 : df(r1)/dr                                                  !
    !------------------------------------------------------------------!

    double precision, intent(in)  :: r, r1, rc, f1, df1
    double precision, intent(out) :: t, dt, d2t, d3t

    double precision :: A, B
    double precision :: k1, k2

    k1 = (rc-r1)
    k2 = k1*k1
    A =  ( 3.0d0*f1 + k1*df1 )/k2
    B = -( 2.0d0*f1 + k1*df1 )/(k1*k2)

    k1 = (rc-r)
    t   =  ( A + B*k1 )*k1*k1
    dt  = -( 2.0d0*A + 3.0d0*B*k1 )*k1
    d2t =  2.0d0*A + 6.0d0*B*k1
    d3t = -6.0d0*B

  end subroutine cutoff_poly3

end program potfit
