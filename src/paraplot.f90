program paraplot

  !--------------------------------------------------------------------!
  ! paraplot -- plot TB parameters                                     !
  !--------------------------------------------------------------------!
  ! This tool can be used to get plottable data files containing the   !
  ! TB parameters of a selected parameter set.                         !
  !--------------------------------------------------------------------!
  ! 2011-03-11 Alexander Urban (AU)                                    !
  ! 2012-01-26 AU --- print out on-site parameters for both directions !
  !--------------------------------------------------------------------!


  use bondint,   only: bi_bondint,          &
                       bi_pairpot,          &
                       bi_get_r_min_pot,    &
                       bi_get_r_max_pot,    &
                       bi_get_r_min,        &
                       bi_get_r_max

  use tbio,      only: tbio_init,           &
                       tbio_final,          &
                       tbio_args,           &
                       args_switch,         &
                       args_get
                
  use tbparam,   only: tbp_init,            &
                       tbp_final,           &
                       tbp_H, tbp_S,        &
                       tbp_O, tbp_pot,      &
                       l_of_il, nl

  implicit none

  !---------------------------- constants -----------------------------!
  ! NLMAX     : max. number of different angular momentum channels per !
  !             atom type                                              !
  !--------------------------------------------------------------------!

  integer,   parameter :: NLMAX =   3
  integer,   parameter :: NPNTS = 100
  integer,   parameter :: u_H   =  40
  integer,   parameter :: u_S   =  41
  integer,   parameter :: u_O   =  42
  integer,   parameter :: u_V   =  43

  !---------------------------- options   -----------------------------!
  ! paramset       : name of the TB parameters set                     !
  ! paramdir       : directory with the TB parameters files            !
  !--------------------------------------------------------------------!
  
  character(len=2), dimension(2) :: atomType
  integer                        :: nTypes

  integer                        :: itype1, itype2, il1, il2
  integer                        :: l1, l2, m
  double precision               :: H, S, O, dH, dS, dO, V, dV

  double precision               :: rmax, rmin
  double precision               :: r, dr

  character(len=100) :: paramset, paramdir
  character(len=100) :: fname

  integer            :: ipos, ipt

  !----------------------------- options ------------------------------!

  call tbio_init(name  = 'paraplot.x',              &
                 descr = 'plot parameter curves',   &
                 args  = (/ 'dir  ','param' /), &
                 extra = '--bond')


  call tbio_args('dir',   paramdir)
  call tbio_args('param', paramset)

  call args_switch('--bond', pos=ipos)
  if (ipos == 0) then
     write(0,*) "Please specify two species using the `--bond' switch."
     call tbio_final()
     stop
  end if
  call args_get(ipos+1, atomType(1))
  call args_get(ipos+2, atomType(2))
  if (atomType(1) == atomType(2)) then
     nTypes = 1
  else
     nTypes = 2
  end if

  call tbio_final()

  !-------------------------- initialization --------------------------!

  call tbp_init(NLMAX, nTypes, atomType, paramset, paramdir)

  rmin = bi_get_r_min()
  rmax = bi_get_r_max()

  dr = (rmax - rmin)/dble(NPNTS - 1)


  !------------------------------ output ------------------------------!

  do itype1 = 1, nTypes
  do itype2 = itype1, nTypes

     ! pair potential:
     if (tbp_pot) then
        fname = trim(paramset)//'_'//trim(atomType(itype1))//'-'&
              //trim(atomType(itype2))//'.dat'
        open(u_V, file='V_'//trim(fname), status='replace', action='write')
        r = rmin
        do ipt = 1, NPNTS
           call bi_pairpot(r, itype1, itype2, V=V, dV=dV)
           write(u_V, '(1x,3(ES15.8,2x))') r, V, dV
           r = r + dr
        end do
        close(u_V)
     end if

     l1_loop : do il1 = 1, nl(itype1)
        l1 = l_of_il(il1, itype1)

        ! on-site parametrization:
        if (tbp_O) then
           do il2 = 1, nl(itype1)
              l2 = l_of_il(il2,itype1)
              do m = 0, min(l1,l2)
                 fname = get_filename(itype1, itype2, l1, l2, m)
                 open(u_O, file='O_'//trim(fname), status='replace', action='write')
                 r = rmin
                 do ipt = 1, NPNTS
                    call bi_bondint(r, itype1, itype2, l1, l2, m, O=O, dO=dO)
                    write(u_O, '(1x,3(ES15.8,2x))') r, O, dO
                    r = r + dr
                 end do
                 close(u_O)
              end do
           end do
        end if

        ! H / S matrix elements:
        do il2 = 1, nl(itype2) 
           l2 = l_of_il(il2, itype2)   
           do m = 0, min(l1,l2)
              fname = get_filename(itype1, itype2, l1, l2, m)
              if (tbp_H) then
                 open(u_H, file='H_'//trim(fname), status='replace', action='write')
              end if
              if (tbp_S) then
                 open(u_S, file='S_'//trim(fname), status='replace', action='write')
              end if
              r = rmin
              do ipt = 1, NPNTS
                 if (tbp_H) then
                    call bi_bondint(r, itype1, itype2, l1, l2, m, H=H, dH=dH)
                    write(u_H, '(1x,3(ES15.8,2x))') r, H, dH
                 end if
                 if (tbp_S) then
                    call bi_bondint(r, itype1, itype2, l1, l2, m, S=S, dS=dS)
                    write(u_S, '(1x,3(ES15.8,2x))') r, S, dS
                 end if
                 r = r + dr
              end do
              if (tbp_H) close(u_H)
              if (tbp_S) close(u_S)
           end do
        end do

     end do l1_loop

  end do
  end do


  !--------------------------- finalization ---------------------------!

  call tbp_final()



contains !=============================================================!

  function get_filename(itype1, itype2, l1, l2, m) result(fname)

    implicit none

    integer, intent(in) :: itype1, itype2
    integer, intent(in) :: l1, l2, m
    character(len=100)  :: fname

    fname = trim(paramset)//'_'//trim(atomType(itype1)) & 
         //'-'//trim(atomType(itype2))//'_'

    select case(l1)
    case (0)
       fname = trim(fname) // 's'
    case (1)
       fname = trim(fname) // 'p'
    case (2)
       fname = trim(fname) // 'd'
    case default
       fname = trim(fname) // 'x'
    end select

    select case(l2)
    case (0)
       fname = trim(fname) // 's'
    case (1)
       fname = trim(fname) // 'p'
    case (2)
       fname = trim(fname) // 'd'
    case default
       fname = trim(fname) // 'x'
    end select

    select case(m)
    case (0)
       fname = trim(fname) // 's'
    case (1)
       fname = trim(fname) // 'p'
    case (2)
       fname = trim(fname) // 'd'
    case default
       fname = trim(fname) // 'x'
    end select

    fname = trim(fname) // '.dat'

  end function get_filename


end program paraplot

