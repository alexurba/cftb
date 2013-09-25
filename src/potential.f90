module potential

  !--------------------------------------------------------------------!
  ! Evaluate real space interaction potentials.  Especially used for   !
  ! the (repulsive) TB pair potential.                                 !
  !--------------------------------------------------------------------!
  ! 2010-11-07 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use bondint,   only: bi_pairpot,       &
                       bi_get_isInit,    &
                       bi_get_hasV,      &
                       bi_get_r_max_pot

  implicit none

  public  :: pot_init,      &
             pot_final,     &
             pot_pairpot


  !--------------------------------------------------------------------!
  ! pot_isInit     : .true., if the module has been initialized        !
  ! pot_hasPairPot : .true., if a TB pair potential is available       !
  !                                                                    !
  ! pot_r_max      : upper radial cut-off for all pot. interactions    !
  !--------------------------------------------------------------------!

  logical,                                      private :: pot_isInit     = .false.
  logical,                                      private :: pot_hasPairPot = .false.

  double precision,                             public  :: pot_r_max
  double precision,                             public  :: pot_r_max2

contains

  !--------------------------------------------------------------------!
  !                    initialization/finalization                     !
  !--------------------------------------------------------------------!

  subroutine pot_init()

    implicit none

    if (pot_isInit) then
       write(0,*) "Warning: module `potential' has already been initialized."
       return
    end if

    if (bi_get_isInit()) then
       pot_r_max      = bi_get_r_max_pot()
       pot_r_max2     = pot_r_max*pot_r_max
       pot_hasPairPot = bi_get_hasV()
    else
       write(0,*) 'Error: bond integral tables not available in ' &
            // "`pot_init()'."
       stop
    end if

    pot_isInit = .true.

  end subroutine pot_init

  !--------------------------------------------------------------------!

  subroutine pot_final()

    implicit none

    pot_isInit = .false.

  end subroutine pot_final

  !--------------------------------------------------------------------!
  !            potential energy from the TB pair potential             !
  !--------------------------------------------------------------------!

  subroutine pot_pairpot(natoms, coord, atomType, avec, nT, T, &
                         Tlen, dmax, energy, forces)

    implicit none

    integer,                               intent(in)  :: natoms
    double precision, dimension(3,natoms), intent(in)  :: coord
    integer,          dimension(natoms),   intent(in)  :: atomType
    double precision, dimension(3,3),      intent(in)  :: avec
    integer,                               intent(in)  :: nT
    integer,          dimension(3,nT),     intent(in)  :: T
    double precision, dimension(nT),       intent(in)  :: Tlen
    double precision,                      intent(in)  :: dmax

    double precision,                                intent(out) :: energy
    double precision, dimension(3,natoms), optional, intent(out) :: forces

    double precision, dimension(3) :: vec
    double precision               :: dmax2
    double precision               :: dist2, Eval, dEval, r
    integer                        :: itype1, itype2, iatom1, iatom2, iT
    logical                        :: calcF

    if (.not. pot_hasPairPot) then
       write(0,*) "Warning: no pair potential parametrization available " &
                  // "in `pot_pairpot()'."
       return
    end if

    if (present(forces)) then
       calcF = .true.
       forces(:,:) = 0.0d0
    else
       calcF = .false.
    end if

    dmax2 = 2.0d0*dmax
    dmax2 = dmax2*dmax2

    energy = 0.0d0
    atom1 : do iatom1 = 1, natoms
       itype1 = atomType(iatom1)

       ! atoms in the unit cell T = (0,0,0)
       do iatom2 = iatom1+1, natoms
          itype2 = atomType(iatom2)

          vec(1:3) = coord(1:3,iatom2) - coord(1:3,iatom1)
          vec(1:3) = vec(1)*avec(1:3,1) &
                   + vec(2)*avec(1:3,2) &
                   + vec(3)*avec(1:3,3)

          dist2 = sum(vec*vec)
          if (dist2 > pot_r_max2) cycle

          r = sqrt(dist2)
          call bi_pairpot(r, itype1, itype2, V=Eval, dV=dEval)
          energy = energy + Eval

          if (calcF) then
             vec(1:3) = vec(1:3)/r
             forces(1:3,iatom1) = forces(1:3,iatom1) + vec(1:3)*dEval
             forces(1:3,iatom2) = forces(1:3,iatom2) - vec(1:3)*dEval
!!$             forces(1:3,iatom1) = forces(1:3,iatom1) - vec(1:3)*dEval
!!$             forces(1:3,iatom2) = forces(1:3,iatom2) + vec(1:3)*dEval
          end if
          
       end do

       ! periodic pictures of the unit cell; T /= (0,0,0):
       translations : do iT = 1, nT

          if (Tlen(iT) - 2.0d0*dmax > pot_r_max) exit translations

          do iatom2 = 1, natoms
             itype2 = atomType(iatom2)

             vec(1:3) = coord(1:3,iatom2) - coord(1:3,iatom1)
             vec(1:3) = vec(1:3) + dble(T(1:3,iT))
             vec(1:3) = vec(1)*avec(1:3,1) &
                      + vec(2)*avec(1:3,2) &
                      + vec(3)*avec(1:3,3)

             dist2 = sum(vec*vec)
             if (dist2 > pot_r_max2) cycle
             
             r = sqrt(dist2)
             call bi_pairpot(r, itype1, itype2, V=Eval, dV=dEval)
             energy = energy + Eval

             if (calcF) then
                vec(1:3) = vec(1:3)/r
                forces(1:3,iatom1) = forces(1:3,iatom1) + vec(1:3)*dEval
                forces(1:3,iatom2) = forces(1:3,iatom2) - vec(1:3)*dEval
!!$                forces(1:3,iatom1) = forces(1:3,iatom1) - vec(1:3)*dEval
!!$                forces(1:3,iatom2) = forces(1:3,iatom2) + vec(1:3)*dEval
             end if

          end do

       end do translations

    end do atom1

  end subroutine pot_pairpot

  

end module potential
