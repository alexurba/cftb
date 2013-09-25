program paramod

  !--------------------------------------------------------------------!
  ! paramod --- modify TB parametrization                              !
  !--------------------------------------------------------------------!
  ! Add a distance dependent potential shift to the TB parameters.     !
  !--------------------------------------------------------------------!
  ! 2012-04-03 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!


  use bondint,   only: bi_bondint,          &
                       bi_pairpot,          &
                       bi_get_r_min_pot,    &
                       bi_get_r_max_pot,    &
                       bi_get_r_min,        &
                       bi_get_r_max

  use interpol,  only: cs_init,             &
                       cs_final,            &
                       cs_interp

  use io,        only: io_unit,             &
                       io_adjustl

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
  
  character(len=100)                                  :: paramset, paramdir
  character(len=100)                                  :: fname

  integer                                             :: u_shft
  character(len=50)                                   :: shftFile
  integer                                             :: nshfts, ishft
  double precision, dimension(:), allocatable, target :: shiftX, shiftY
  double precision                                    :: offset

  character(len=2), dimension(2)                      :: atomType
  integer                                             :: nTypes

  integer                                             :: itype1, itype2, il1, il2
  integer                                             :: l1, l2, m
  double precision                                    :: H, S, O, dH, dS, dO, V, dV

  double precision                                    :: rmax, rmin
  double precision                                    :: r, dr

  integer                                             :: ipos, ipt
  logical                                             :: fexists

  !----------------------------- options ------------------------------!

  call tbio_init(name  = 'paramod.x',              &
                 descr = 'modify TB parameters',   &
                 args  = (/ 'dir  ','param' /),    &
                 extra = '--bond:--shift')


  call tbio_args('dir',   paramdir)
  call tbio_args('param', paramset)

  ! specify the bond for which the params shall be altered:
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

  ! file with potential shifts:
  call args_switch('--shift', pos=ipos)
  if (ipos == 0) then
     write(0,*) "Please specify a file with potential shifts"
     write(0,*) "using the `--shift' switch."
     call tbio_final()
     stop
  end if
  call args_get(ipos+1, shftFile)
  inquire(file=trim(shftFile), exist=fexists)
  if (.not. fexists) then
     write(0,*) "Error: File not found: ", trim(shftFile)
     call tbio_final()
     stop
  end if

  call tbio_final()


  !----------- read potential shifts and init interpolator ------------!

  call read_potential_shifts()
  call cs_init(shiftX, shiftY, nshfts)


  !--------------------- initialize TB parameters ---------------------!

  call tbp_init(NLMAX, nTypes, atomType, paramset, paramdir)

  rmin = max(bi_get_r_min(), shiftX(1))
  rmax = max(bi_get_r_max(), shiftX(nshfts))
  offset = cs_interp(rmax)

  dr = (rmax - rmin)/dble(NPNTS - 1)


  !------------------------------ output ------------------------------!

  ! the interpolated shift data:
  u_shft = io_unit()
  open(u_shft, file=trim(shftFile)//'-interp', status='replace', action='write')
  r = rmin
  do ipt = 1, NPNTS
     write(u_shft, '(1x,2(ES15.8,2x))') r, cs_interp(r)
     r = r + dr
  end do
  close(u_shft)

  do itype1 = 1, nTypes
  do itype2 = itype1, nTypes

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
                    if ((l1 == l2) .and. (r <= shiftX(nshfts))) then
                       write(u_O, '(1x,2(ES15.8,2x))') r, O + cs_interp(r)-offset
                    else
                       write(u_O, '(1x,2(ES15.8,2x))') r, O
                    end if
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
              r = rmin
              do ipt = 1, NPNTS
                 if (tbp_H) then
                    call bi_bondint(r, itype1, itype2, l1, l2, m, H=H, dH=dH)
                    call bi_bondint(r, itype1, itype2, l1, l2, m, S=S, dS=dS)
                    if (r <= shiftX(nshfts)) then
                       write(u_H, '(1x,2(ES15.8,2x))') r, H + cs_interp(r)*S
                    else
                       write(u_H, '(1x,2(ES15.8,2x))') r, H
                    end if
                 end if
                 r = r + dr
              end do
              if (tbp_H) close(u_H)
           end do
        end do

     end do l1_loop

  end do
  end do

  !--------------------- further output to stdout ---------------------!

  write(*,*) "Offset of the atomic eigenvalues (Ha): ", &
              trim(io_adjustl(offset,8))
  write(*,*)

  !--------------------------- finalization ---------------------------!

  call tbp_final()
  call finalize()


contains !=============================================================!

  subroutine finalize()

    implicit none

    call cs_final

    if (allocated(shiftX)) then
       deallocate(shiftX, shiftY)
    end if

  end subroutine finalize

  !--------------------------------------------------------------------!
  !                  read potential shifts from file                   !
  !--------------------------------------------------------------------!

  subroutine read_potential_shifts()

    implicit none

    integer            :: status
    character(len=100) :: line

    u_shft = io_unit()
    open(u_shft, file=trim(shftFile), status='old', action='read')
    nshfts = 0
    countloop : do
       read(u_shft,'(A)',iostat=status) line
       if (status /= 0) exit countloop
       if (len_trim(line) == 0) cycle countloop
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle countloop
       nshfts = nshfts + 1
    end do countloop
    allocate(shiftX(nshfts),shiftY(nshfts))
    rewind(u_shft)
    ishft = 0
    readloop : do
       if (ishft >= nshfts) exit readloop
       read(u_shft,'(A)') line
       if (len_trim(line) == 0) cycle readloop
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle readloop
       ishft = ishft + 1
       read(line,*) shiftX(ishft), shiftY(ishft)
    end do readloop
    close(u_shft)

  end subroutine read_potential_shifts

  !--------------------------------------------------------------------!

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

    fname = trim(fname) // '.dat-shft'

  end function get_filename


end program paramod

