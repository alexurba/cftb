module timing

  implicit none
  save

  public :: tng_init,    &
            tng_final,   &
            tng_timing,  &
            tng_timing2

  double precision,            public  :: tng_t_ini
  double precision,            public  :: tng_t_prev
  double precision,            public  :: tng_t_now
  double precision,            public  :: tng_t

  !--------------------------------------------------------------------!

  character(len=*), parameter, private :: tng_file   = 'timing.log'
  integer,          parameter, private :: u_tng      = 88
  logical,                     private :: tng_isinit = .false.
  logical,                     private :: tng_output = .false.

contains

  subroutine tng_init()

    implicit none

    open(u_tng, file=tng_file, status='replace', action='write')
    tng_output = .true.
    call tng_timing()

  end subroutine tng_init

  !--------------------------------------------------------------------!

  subroutine tng_final()

    implicit none

    call tng_timing('Timing finished.')
    close(u_tng)

  end subroutine tng_final

  !--------------------------------------------------------------------!
  !           Timing with respect to the initial time t_ini            !
  !--------------------------------------------------------------------!

  subroutine tng_timing(msg)

    implicit none

    character(len=*), optional, intent(in) :: msg
    integer                                :: cnt, cnt_rate

    call system_clock(count=cnt, count_rate=cnt_rate)
    tng_t_now = dble(cnt)/dble(cnt_rate)

    if (.not. tng_isinit) then
       tng_t_ini  = tng_t_now
       tng_t_prev = tng_t_now
       tng_isinit = .true.
       tng_t      = 0.0d0
       return
    endif

    tng_t      = tng_t + tng_t_now - tng_t_prev
    tng_t_prev = tng_t_now

    if (tng_output) then
       if (present(msg)) then
          write(u_tng, '(1x,F10.2," s",5x,A)') tng_t, msg
       else
          write(u_tng, '(1x,F10.2," s")') tng_t
       end if
    else
       if (present(msg)) then
          write(*, '(1x,F10.2," s",5x,A)') tng_t, msg
       else
          write(*, '(1x,F10.2," s")') tng_t
       end if
    end if

  end subroutine tng_timing

  !--------------------------------------------------------------------!
  !  Measure just the time passed since the last call to tng_timing()  !
  !--------------------------------------------------------------------!

  subroutine tng_timing2(msg)

    implicit none

    character(len=*), optional, intent(in) :: msg
    integer                                :: cnt, cnt_rate

    call system_clock(count=cnt, count_rate=cnt_rate)
    tng_t_now = dble(cnt)/dble(cnt_rate)

    if (.not. (tng_isinit .and. tng_output)) then
       return
    endif

    if (present(msg)) then
       write(u_tng, '(5x,F10.2," s",1x,A)') tng_t_now - tng_t_prev, msg
    else
       write(u_tng, '(5x,F10.2," s")') tng_t_now - tng_t_prev
    end if

  end subroutine tng_timing2


end module timing
