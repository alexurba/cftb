program tb

  !--------------------------------------------------------------------!
  ! tb.x -- tight-binding calculations                                 !
  !--------------------------------------------------------------------!
  ! 2011-03-19 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use bfgs,      only: bfgs_init,             &
                       bfgs_final,            &
                       bfgs_optimize

  use geometry,  only: geo_print_coordinates, &
                       cooLattice,            &
                       doForces,              &
                       doGeoOpt,              &
                       nAtoms,                &
                       optIter,               &
                       optTol,                &
                       optAlpha,              &
                       optAll,                &
                       optNumActive,          &
                       optActive,             &
                       optNumDoF,             &
                       optDoF

  use io,        only: io_adjustl

  use tbinter,   only: tbi_init,              &
                       tbi_final,             &
                       tbi_reset,             &
                       tbi_energy_analysis,   &
                       tbi_energies,          &
                       tbi_forces,            &
                       tbi_kpoint_sampling,   &
                       tbi_print_eigenvalues, &
                       tbi_print_energies,    &
                       tbi_print_forces,      &
                       tbi_print_mulliken,    &
                       tbi_print_occupation,  &
                       tbi_print_system,      &
                       tbi_setup_T0

  use tbio,      only: tbio_init,             &
                       tbio_final,            &
                       tbio_args,             &
                       args_switch

  implicit none

  !----------------------------- energies -----------------------------!
  ! E_fermi    : Fermi energy                                          !
  ! E_smear    : energy correction due to smearing of the Fermi edge   !
  ! E_band     : band sum = sum of eigenvalues of the occupied bands   !
  ! E_atoms    : sum of the energies of the free atoms in the system   !
  ! E_pairpot  : energy due to the pair potential                      !
  ! E_coh      : cohesive energy (TB-Band model)                       !
  !                                                                    !
  ! Energy partitioning:                                               !
  !                                                                    !
  ! E_prom     : promotion energy                                      !
  ! E_cf       : crystal field energy                                  !
  ! E_polar    : polarization energy                                   !
  ! E_cov      : covalent bond energy                                  !
  ! E_bond     : cohesive energy in the TB-Bond model                  !
  !                                                                    !
  !------------------------ files and options -------------------------!
  ! infile     : input file name                                       !
  ! outfile    : output file name                                      !
  ! occupfile  : file with band occupations (optional)                 !
  ! infiletype : file format of the input file                         !
  ! paramset   : name of the TB parameter set                          !
  ! paramdir   : directory path to look for TB parameter files         !
  ! timing     : if .true. --> collect timing information              !
  !--------------------------------------------------------------------!

  double precision   :: E_fermi, E_smear
  double precision   :: E_band, E_atoms, E_pairpot
  double precision   :: E_coh

  double precision   :: E_prom, E_cf, E_polar, E_cov, E_bond

  character(len=100) :: infile, outfile, occupfile
  character(len=100) :: infiletype
  character(len=100) :: paramset, paramdir
  logical            :: timing

  integer            :: niters, iter
  double precision   :: F_max, F_rms
  logical            :: conv

  logical            :: energy_analysis   = .true.

  logical            :: print_energies    = .true.
  logical            :: print_forces      = .true.
  logical            :: print_charges     = .true.
  logical            :: print_eigenvalues = .true.
  logical            :: print_occupation  = .true.

  character(len=512) :: flags

  double precision, dimension(:,:), pointer   :: forCart
  double precision, dimension(:), allocatable :: X, F

  !----------------------------- options ------------------------------!

  flags = '--eanalysis:--printoccup:--printev:--printcharge'
  
  call tbio_init(name  = 'tb.x',                                    &
                 descr = 'run tight-binding calculations',          & 
                 args  = (/ 'in    ', 'out   ', 'format', 'param ', &
                            'dir   ', 'time  ', 'occup ' /),        &
                 extra = trim(flags) )

  call tbio_args('in',     infile)
  call tbio_args('out',    outfile)
  call tbio_args('format', infiletype)
  call tbio_args('param',  paramset)
  call tbio_args('time',   timing)
  call tbio_args('occup',  occupfile)
  call tbio_args('dir',    paramdir)

  call args_switch('--eanalysis',   value=energy_analysis)

  call args_switch('--printoccup',  value=print_occupation)
  call args_switch('--printev',     value=print_eigenvalues)
  call args_switch('--printcharge', value=print_charges)

  call tbio_final()

  !-------------------------- initialization --------------------------!

  call tbi_init(infile, infiletype, paramset, paramdir, outfile, &
                occupfile, timing, matrices=.false.)

  call tbi_print_system()
    
  if (doGeoOpt) then
     niters     = optIter
     if (.not. optAll) allocate(X(optNumDoF), F(optNumDoF))
     call bfgs_init(optNumDoF, alpha=optAlpha)
     call print_geoopt_header(niters, optTol, optNumDoF)
  else
     niters = 1
  end if

  !--------------------------- calculation ----------------------------!

  iter = 0
  optimize : do
     iter = iter + 1

     call tbi_reset()
     call tbi_setup_T0()
     call tbi_kpoint_sampling()
     call tbi_energies(E_fermi, E_smear, E_band, E_atoms, E_pairpot, E_coh)

!DEBUG
     if (print_forces) call tbi_print_forces()
!END DEBUG

     if (doForces) call tbi_forces(F=forCart)

     if (doGeoOpt) then
        if (.not. optAll) then
           call compress_data(optNumActive, optActive, optDoF, nAtoms, &
                              cooLattice, forCart, optNumDof, X, F)
           call bfgs_optimize(optNumDoF, X, F, optTol, &
                              F_max, F_rms, conv)
           call update_coordinates(optNumActive, optActive, optDoF, &
                                   optNumDof, X, nAtoms, cooLattice)
        else
           call bfgs_optimize(optNumDoF, cooLattice, forCart, optTol, &
                              F_max, F_rms, conv)
        end if
        write(*,'(1x,I3,3x,3(ES15.8,1x))') iter, E_coh, F_max, F_rms
        if (conv) exit optimize
     end if

     if (iter >= niters) exit optimize
  end do optimize

  if (energy_analysis) then
     call tbi_energy_analysis(E_prom, E_cf, E_polar, E_cov, E_bond)
  end if

  !------------------------------ output ------------------------------!

  if (doGeoOpt) then
     write(*,*)
     if (conv) then
        write(*,*) 'Converged after ' // trim(io_adjustl(iter)) // &
                   ' iterations.'
     else
        write(*,*) 'NOT converged after ' // trim(io_adjustl(iter)) // &
                   ' iterations.'
     end if
     write(*,*)
     write(*,*)
  end if

  if (print_eigenvalues) call tbi_print_eigenvalues()
  if (print_occupation)  call tbi_print_occupation()
  if (print_charges)     call tbi_print_mulliken()
  if (print_energies)    call tbi_print_energies()
  if (print_forces)      call tbi_print_forces()
  if (doGeoOpt)          call geo_print_coordinates()
  
  !--------------------------- finalization ---------------------------!

  if (doGeoOpt) then
     call bfgs_final()
  end if
  call tbi_final()

contains

  subroutine print_geoopt_header(niters, tol, nDoF)

    implicit none

    integer,          intent(in) :: niters
    double precision, intent(in) :: tol
    integer,          intent(in) :: nDoF

    write(*,*) 'Geometry Optimization'
    write(*,*) '====================='
    write(*,*) 
    write(*,'(1x,"Max. optimization steps : ",A)') trim(io_adjustl(niters))
    write(*,'(1x,"Optimization threshold  : ",ES8.2)') tol
    write(*,'(1x,"Degrees of freedom      : ",A)') trim(io_adjustl(nDoF))
    write(*,*) 
    write(*,'(7x,3(A15,1x))') 'coh. energy  ', 'max. force  ', &
                              'RMS force   '
    write(*,'(7x,3(A15,1x))') '[Ha]     ', '[Ha/Bohr]   ', &
                              '[Ha/Bohr]   '

  end subroutine print_geoopt_header

  !--------------------------------------------------------------------!
  !        geometry optimization for certain degrees of freedom        !
  !--------------------------------------------------------------------!

  subroutine compress_data(nActive, active, DoF, nAtoms, X_in, F_in, n, X, F)

    implicit none

    integer,                                intent(in)  :: nActive
    integer,          dimension(nActive),   intent(in)  :: active
    integer,          dimension(3,nActive), intent(in)  :: DoF
    integer,                                intent(in)  :: nAtoms
    double precision, dimension(3,nAtoms),  intent(in)  :: X_in, F_in
    integer,                                intent(in)  :: n
    double precision, dimension(n),         intent(out) :: X, F

    integer :: iact, iat, idim, ix

    ix = 0
    do iact = 1 , nActive
       iat = active(iact)
       do idim = 1, 3
          if (DoF(idim,iact) == 1) then
             ix = ix + 1
             X(ix) = X_in(idim,iat)
             F(ix) = F_in(idim,iat)
          end if
       end do
    end do

  end subroutine compress_data

  !--------------------------------------------------------------------!

  subroutine update_coordinates(nActive, active, DoF, n, X, nAtoms, X_out)

    implicit none

    integer,                                intent(in)  :: nActive
    integer,          dimension(nActive),   intent(in)  :: active
    integer,          dimension(3,nActive), intent(in)  :: DoF
    integer,                                intent(in)  :: n
    double precision, dimension(n),         intent(in)  :: X
    integer,                                intent(in)  :: nAtoms
    double precision, dimension(3,nAtoms),  intent(out) :: X_out

    integer :: iact, iat, idim, ix

    ix = 0
    do iact = 1 , nActive
       iat = active(iact)
       do idim = 1, 3
          if (DoF(idim,iact) == 1) then
             ix = ix + 1
             X_out(idim,iat) = X(ix)
          end if
       end do
    end do

  end subroutine update_coordinates

end program tb
