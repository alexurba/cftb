module energy

  !--------------------------------------------------------------------!
  !                     energy related procedures                      !
  !--------------------------------------------------------------------!
  ! 2010-11-03 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use constants, only: SQRT_PI, SQRT_PI_INV
  use funclib,   only: erfc
  use io,        only: io_adjustl
  use sortlib,   only: argsort
  
  implicit none

contains

  !--------------------------------------------------------------------!
  !                     determine the Fermi level                      !
  !--------------------------------------------------------------------!

  subroutine nrg_fermi_level(nelecs, nkpts, wkpt, nevs, eval, metallic, &
                             smearwidth, Ef, corr, iter)

    implicit none

    !------------------------------------------------------------------!
    ! nelecs          : number of electrons                            !
    ! nkpts           : number of k-points                             !
    ! wkpt(ikpt)      : weight of ikpt-th k-point                      !
    ! nevs            : number of eigenvalues = number of bands        !
    ! eval(iev,ikpt)  : iev-th eigenvalue at k-point ikpt              !
    ! metallic        : .true., if system has no bandgap               !
    ! smearwidth      : energy width for smearing of the Fermi edge    !
    ! Ef              : Fermi energy                                   !
    ! corr            : broadening correction                          !
    ! iter            : iterations to convergence (= -1 if failed)     !
    !------------------------------------------------------------------!

    integer,                                  intent(in)  :: nelecs
    integer,                                  intent(in)  :: nkpts
    double precision, dimension(nkpts),       intent(in)  :: wkpt
    integer,                                  intent(in)  :: nevs
    double precision, dimension(nevs, nkpts), intent(in)  :: eval
    logical,                                  intent(in)  :: metallic
    double precision,                         intent(in)  :: smearwidth
    double precision,                         intent(out) :: Ef
    double precision,                         intent(out) :: corr
    integer,                                  intent(out) :: iter

    integer,          parameter :: NITMAX = 5000
    double precision, parameter :: CONV   = 1.0d-5
    double precision, parameter :: DE0    = 0.001d0

    double precision :: ndiff, ndiff0
    double precision :: E, Ef_try, dE, nstates
    double precision :: smearwidth_inv

    integer          :: ikpt, iev

    if ((nkpts < 1) .or. (nevs < 1)) then
       write(0,*) "Error: no eigenvalues available in `nrg_fermi_level()'."
       stop
    end if

    ! start with a guess:
    Ef_try = eval(nelecs/2, 1)
    Ef     = Ef_try
    dE     = DE0
    corr   = 0.0d0

    iter = 0
    if ((.not. metallic) .and. (nkpts==1)) return

    ! search fo Ef in one direction:
    iter   = 1
    ndiff0 = 0.0d0
    search : do
       nstates = 0.0d0
       do ikpt = 1, nkpts
       do iev  = 1, nevs
          E = eval(iev, ikpt)
          nstates = nstates + wkpt(ikpt)*nrg_occupancy(E, Ef_try, metallic, smearwidth)
       end do
       end do
       ndiff = dble(nelecs) - 2.0d0*nstates
       if (abs(ndiff) < CONV) then
          Ef = Ef_try
          return
       else if (ndiff0 == 0.0d0) then
          ndiff0 = ndiff
       else if (ndiff0*ndiff < 0.0d0) then
          exit search
       else if (ndiff < 0.0d0) then
          dE     = 2.0d0*dE
          Ef_try = Ef_try - dE
       else
          dE     = 2.0d0*dE
          Ef_try = Ef_try + dE
       end if
       iter = iter + 1
       if (iter > NITMAX) then
          write(0,*) "Error: unable to determine Fermi level " &
                     // "in `nrg_fermi_level()' (1)."
          stop
       end if
    end do search

    ! refine with nested intervals:
    refine : do

       iter = iter + 1

       dE = 0.5d0*dE
       if (ndiff < 0.0d0) then
          Ef_try = Ef_try - dE
       else
          Ef_try = Ef_try + dE
       end if

       nstates = 0.0d0
       do ikpt = 1, nkpts
       do iev  = 1, nevs
          E = eval(iev, ikpt)
          nstates = nstates + wkpt(ikpt)*nrg_occupancy(E, Ef_try, metallic, smearwidth)
       end do
       end do
       ndiff = dble(nelecs) - 2.0d0*nstates

       if (abs(ndiff) < CONV) then
          exit refine
       end if

       if (iter > NITMAX) then
          write(0,*) "Error: unable to determine Fermi level " &
                     // "in `nrg_fermi_level()' (2)."
          stop
       end if
    end do refine
 
    Ef = Ef_try

    ! broadening correction:
    corr = 0.0d0
    if (metallic) then
       smearwidth_inv = 1.0d0/smearwidth
       do ikpt = 1, nkpts
          do iev  = 1, nevs
             E = (Ef - eval(iev, ikpt))*smearwidth_inv
             corr = corr + wkpt(ikpt)*exp(-E*E)
          end do
       end do
       corr = -smearwidth*SQRT_PI_INV*corr
    end if

  end subroutine nrg_fermi_level

  !--------------------------------------------------------------------!
  !                   ocuupation number distribution                   !
  !--------------------------------------------------------------------!

  function nrg_occupancy(E, flevel, metallic, smearwidth) result(f)

    implicit none

    double precision, intent(in) :: E
    double precision, intent(in) :: flevel
    logical,          intent(in) :: metallic
    double precision, intent(in) :: smearwidth
    double precision             :: f

    double precision, parameter  :: EPS_here = 1.0d-6

    if (metallic) then
       f = 0.5d0*(erfc((E - flevel)/smearwidth))
    else if ((E - flevel) < EPS_here) then
       f = 1.0d0
    else
       f = 0.0d0
    end if

  end function nrg_occupancy

  !--------------------------------------------------------------------!
  !      gaussian broadening of the delta function delta(ev - E)       !
  !--------------------------------------------------------------------!

  function nrg_gaussian(E, ev, smearwidth) result(g)

    implicit none

    double precision, intent(in) :: E
    double precision, intent(in) :: ev
    double precision, intent(in) :: smearwidth
    double precision             :: g

    double precision :: arg, sw
    
    sw  = 1.0d0/smearwidth
    arg = (E - ev)*sw
    g   = SQRT_PI_INV*sw*exp(-arg*arg)

  end function nrg_gaussian

  !--------------------------------------------------------------------!
  !                      sum of atomic energies                        !
  !--------------------------------------------------------------------!

  subroutine nrg_atomic(nTypes, nAtomsOfType, nElec, nlmax, nl, level, &
                        l_of_il, Eatoms, Qatoms)

    implicit none

    integer,                                   intent(in)  :: nTypes
    integer,          dimension(nTypes),       intent(in)  :: nAtomsOfType
    integer,          dimension(nTypes),       intent(in)  :: nElec
    integer,                                   intent(in)  :: nlmax
    integer,          dimension(nTypes),       intent(in)  :: nl
    double precision, dimension(nlmax,nTypes), intent(in)  :: level
    integer,          dimension(nlmax,nTypes), intent(in)  :: l_of_il
    double precision,                          intent(out) :: Eatoms
    double precision, dimension(nlmax,nTypes), intent(out) :: Qatoms

    integer, dimension(nlmax) :: idx
    integer                   :: itype, l, il, il2, nel, norb
    double precision          :: f_occup
    double precision          :: Etype

    Eatoms      = 0.0d0
    Qatoms(:,:) = 0.0d0

    do itype = 1, nTypes
       Etype = 0.0d0
       nel   = nElec(itype)
       call argsort(level(1:nl(itype),itype),idx(1:nl(itype)))
       do il2 = 1, nl(itype)

          il      = idx(il2)
          l       = l_of_il(il, itype)
          norb    = 2*l + 1
          f_occup = dble(min(nel, 2*norb))

          Etype            = Etype + f_occup*level(il, itype)
          Qatoms(il,itype) = f_occup

          nel = nel - f_occup
          if (nel == 0) exit

       end do
       Eatoms = Eatoms + dble(nAtomsOfType(itype))*Etype
    end do
    
  end subroutine nrg_atomic
  
end module energy


!!$  function nrg_atomic(nTypes, nAtomsOfType, nElec, nlmax, nl, level) result(E)
!!$
!!$    implicit none
!!$
!!$    integer,                                   intent(in) :: nTypes
!!$    integer,          dimension(nTypes),       intent(in) :: nAtomsOfType
!!$    integer,          dimension(nTypes),       intent(in) :: nElec
!!$    integer,                                   intent(in) :: nlmax
!!$    integer,          dimension(nTypes),       intent(in) :: nl
!!$    double precision, dimension(nlmax,nTypes), intent(in) :: level
!!$    double precision                                      :: E
!!$
!!$    integer, dimension(nlmax) :: idx
!!$    integer                   :: itype, il, nel
!!$    double precision          :: Etype
!!$
!!$    E = 0.0d0
!!$
!!$    do itype = 1, nTypes
!!$       Etype = 0.0d0
!!$       nel   = nElec(itype)
!!$       call argsort(level(1:nl(itype),itype),idx(1:nl(itype)))
!!$       do il = 1, nl(itype)
!!$          Etype = Etype + dble(min(nel,2))*level(idx(il), itype)
!!$          nel = nel - min(nel,2)
!!$          if (nel == 0) exit
!!$       end do
!!$       E = E + dble(nAtomsOfType(itype))*Etype
!!$    end do
!!$    
!!$  end function nrg_atomic
