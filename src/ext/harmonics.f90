module harmonics

  use funclib, only: P_lgndr

  implicit none

  public  :: K_cubharm,    &
             K_cubharm1,   &
             K_cubharm2,   &
             rot_init,     &
             rot_final,    &
             rot_integral, &
             plot3d_open,  &
             plot3d_iter,  &
             plot3d_write, &
             rot_wigner_d

  private :: K_prefac 
             

  double precision, parameter, public :: pi          = 3.141592653589793d0
  double precision, parameter, public :: sqrt_pi     = 1.772453850905516d0
  double precision, parameter, public :: pi_inv      = 1.0d0 / pi
  double precision, parameter, public :: sqrt_pi_inv = 1.0d0 / sqrt_pi
  double precision, parameter, public :: pi2         = 2.0d0 * pi
  double precision, parameter, public :: sqrt2       = 1.4142135623731d0
  double precision, parameter, public :: sqrt2_inv   = 1.0d0/sqrt2
  double precision, parameter, public :: EPS         = 1.0d-15

  !--------------------------------------------------------------------!
  !                          cubic harmonics                           !
  !--------------------------------------------------------------------!
  
  interface K_cubharm
     module procedure K_cubharm1, K_cubharm2
  end interface

  !--------------------------------------------------------------------!
  !                         integeral rotation                         !
  !--------------------------------------------------------------------!

  logical,                                     private :: rot_is_init = .false.
  integer,                                     private :: rot_l_max
  integer,          dimension(:), allocatable, private :: rot_fact

  !--------------------------------------------------------------------!
  !                           visualization                            !
  !--------------------------------------------------------------------!

  logical,                                     private :: p3d_open = .false.
  integer,                                     private :: p3d_nx, p3d_ny, p3d_nz
  integer,                                     private :: p3d_ix, p3d_iy, p3d_iz
  double precision, dimension(:), allocatable, private :: p3d_cx, p3d_cy, p3d_cz
  integer, parameter,                          private :: p3d_unit = 97
  

contains

  !--------------------------------------------------------------------!
  !                          cubic harmonics                           !
  !                                                                    !
  ! K_{l,m}(\theta,\phi) = C \cdot P_l^{|m|}(cos\theta) \cdot f(m\phi) !
  !                                                                    !
  !   for m >= 0 : f(m\phi) =  \cos(m\phi)                             !
  !       m <  0 : f(m\phi) = -\sin(m\phi)                             !
  !                                                                    !
  !   C = \sqrt{ \frac{2l + 1}{2\pi(1 + \delta_{m,0})}                 !
  !              \cdot \frac{(l - |m|)!}{(l + |m|)!}    }              !
  !                                                                    !
  !--------------------------------------------------------------------!

  function K_cubharm1(l, m, theta, phi) result(K_lm)

    implicit none

    integer,          intent(in) :: l, m
    double precision, intent(in) :: theta, phi
    double precision             :: K_lm

    if ( (theta < 0) .or.(theta > pi) .or. &
         (phi < 0) .or. (phi > pi2) ) then
       write(0,*) 'Error: invalid coordinates in K_cubharm1'
       stop
    end if

    if ( (l<0) .or. (abs(m) > l)) then
       write(0,*) 'Error: invalid parameters in K_cubharm1'
       stop
    end if

    K_lm = K_prefac(l, abs(m)) * P_lgndr(l,abs(m),cos(theta))
    if (m < 0) then
       K_lm = K_lm * sin(abs(m)*phi)
    else
       K_lm = K_lm * cos(m*phi)
    end if

  end function K_cubharm1

  !--------------------------------------------------------------------!

  function K_cubharm2(l, m, er) result(K_lm)

    implicit none

    integer,                        intent(in) :: l, m
    double precision, dimension(3), intent(in) :: er
    double precision             :: K_lm
    double precision             :: theta, phi

    if (abs(er(1)*er(1)+er(2)*er(2)+er(3)*er(3) - 1.0d0) > EPS) then
       write(0,*) 'Error: K_cubharm2 expects nomalized vector (ex,ey,ez).'
       stop
    end if

    theta = acos(er(3))
    phi   = atan2(er(2),er(1))
    if (phi < 0.0d0) phi = phi + pi2

    K_lm  = K_cubharm1(l, m, theta, phi)
    
  end function K_cubharm2

  !--------------------------------------------------------------------!

  function K_prefac(l, m) result(pref)

    implicit none

    integer, intent(in) :: l, m
    double precision    :: pref

    double precision    :: tmp
    integer             :: im, fact

    if ((m<0) .or. (m>l)) then
       write(0,*) 'Error: invalid parameters in K_prefac.'
       stop
    end if

    tmp = 2.0d0*l + 1.0d0
    if (m == 0) then
       pref = sqrt_pi_inv*sqrt(tmp/(tmp+1.0d0))
    else if (m == l) then
       pref = sqrt_pi_inv*sqrt(tmp/(tmp-1.0d0))
       fact = 1
       do im = 2, l+m
          fact = fact * im
       end do
       pref = pref / dble(fact)
    else
       pref = sqrt_pi_inv*sqrt(tmp/(tmp-1.0d0))
       fact = 1
       do im = l-m, l+m
          fact = fact * im
       end do
       pref = pref / dble(fact)
    end if

  end function K_prefac

  !--------------------------------------------------------------------!
  !                Rotation of Slater-Koster integrals                 !
  !--------------------------------------------------------------------!

  subroutine rot_init(l_max)
    !! tabulate factorials

    implicit none

    integer, intent(in) :: l_max
    integer             :: l

    if (rot_is_init) then
       write(0,*) 'Error: double call of rot_init.'
       return
    end if

    ! remember max. l:
    rot_l_max = l_max

    ! tabulate necessary factorials:
    allocate(rot_fact(0:2*l_max))
    rot_fact(0) = 1
    do l = 1, 2*l_max
       rot_fact(l) = rot_fact(l-1)*l
    end do

    rot_is_init = .true.

  end subroutine rot_init

  !--------------------------------------------------------------------!

  subroutine rot_final()

    implicit none

    if (.not. rot_is_init) then
       write(0,*) 'Warning: nothing to finalize in rot_final.'
       return
    end if

    deallocate( rot_fact )

    rot_is_init = .false.

  end subroutine rot_final

  !--------------------------------------------------------------------!
  ! <l1m1|H|l2m2> = sum_{m=1}^{min(l1,l2)} c_m * (l1,l2,|m|)           !
  !                                                                    !
  ! Ref.: PRB 72 (2005) 165107 and PRB 69 (2004) 233101                !
  !--------------------------------------------------------------------!

  subroutine rot_integral(alpha, cosbeta, l1, m1, l2, m2, n_m, c_m)

    implicit none

    double precision,                   parameter :: EPS_here = 1.0d-6

    double precision,                 intent(in)  :: alpha, cosbeta
    integer,                          intent(in)  :: l1, l2, m1, m2
    integer,                          intent(in)  :: n_m
    double precision, dimension(n_m), intent(out) :: c_m

    integer :: l_min, mp

    if ( (abs(m1)>l1) .or. (abs(m2)>l2) .or. (n_m /= min(l1,l2)+1) &
       .or. (alpha<0.0d0) .or. (alpha>=pi2) .or. (abs(cosbeta)>1.0d0) ) then
       write(0,*) 'Error: inadequate parameters in rot_integral.'
       write(0,*) 'alpha   = ', alpha, alpha-pi2
       write(0,*) 'cosbeta = ', cosbeta
       write(0,*) 'l1, m1  = ', l1, m1
       write(0,*) 'l2, m2  = ', l2, m2
       write(0,*) 'n_m     = ', n_m
       stop
    end if

    l_min = min(l1,l2)
    c_m(1:l_min+1) = 0.0d0

    ! mp == 0
    c_m(1) = 2.0d0*rot_A(m1,alpha)*rot_A(m2,alpha)    &
                  *rot_wigner_d(l1,abs(m1),0,cosbeta) &
                  *rot_wigner_d(l2,abs(m2),0,cosbeta)
    ! mb > 0
    do mp = 1, l_min
       c_m(mp+1) = rot_S(l1,m1,mp,alpha,cosbeta) &
                 * rot_S(l2,m2,mp,alpha,cosbeta) &
                 + rot_T(l1,m1,mp,alpha,cosbeta) &
                 * rot_T(l2,m2,mp,alpha,cosbeta)
    end do

  end subroutine rot_integral

  !--------------------------------------------------------------------!

  function rot_S(l, m, mp, alpha, cosbeta) result(S)

    implicit none

    integer,          intent(in) :: l, m, mp
    double precision, intent(in) :: alpha, cosbeta
    double precision             :: S

    S = rot_wigner_d(l, abs(m), mp, cosbeta)
    if (mod(mp,2) == 1) S = -S
    S = S + rot_wigner_d(l, abs(m), -mp, cosbeta)
    S = S*rot_A(m,alpha)

  end function rot_S

  !--------------------------------------------------------------------!

  function rot_T(l, m, mp, alpha, cosbeta) result(T)

    implicit none

    integer,          intent(in) :: l, m, mp
    double precision, intent(in) :: alpha, cosbeta
    double precision             :: T

    T = rot_wigner_d(l, abs(m), mp, cosbeta)
    if (mod(mp,2) == 1) T = -T
    T = T - rot_wigner_d(l, abs(m), -mp, cosbeta)
    T = T*rot_B(m,alpha)

  end function rot_T
  
  !--------------------------------------------------------------------!

  function rot_A(m, alpha) result(A)

    implicit none

    integer,          intent(in) :: m
    double precision, intent(in) :: alpha
    double precision             :: A

    double precision, parameter  :: SQRT2_inv = 0.707106781186547d0

    if (m < 0) then
       A = sin(dble(m)*alpha)
       if (abs(mod(m,2)) == 0) A = -A
    else if (m > 0) then
       A = cos(dble(m)*alpha)
       if (abs(mod(m,2)) == 1) A = -A
    else
       A = SQRT2_inv
    end if

  end function rot_A

  !--------------------------------------------------------------------!

  function rot_B(m, alpha) result(B)

    implicit none

    integer,          intent(in) :: m
    double precision, intent(in) :: alpha
    double precision             :: B

    if (m < 0) then
       B = cos(dble(m)*alpha)
       if (abs(mod(m,2)) == 1) B = -B
    else if (m > 0) then
       B = sin(dble(m)*alpha)
       if (abs(mod(m,2)) == 0) B = -B
    else
       B = 0.0d0
    end if

  end function rot_B

  !--------------------------------------------------------------------!
  ! Call this procedure for values of the Wigner small d function.     !
  ! Tabulated values will be used where available; an analytic         !
  ! function is called else.                                           !
  !--------------------------------------------------------------------!

  function rot_wigner_d(l, m1_in, m2_in, cosbeta) result(d_wigner)

    implicit none

    integer,          intent(in) :: l, m1_in, m2_in
    double precision, intent(in) :: cosbeta
    double precision             :: d_wigner  

    integer                      :: m1, m2
    double precision             :: pre, rtmp

    ! if cos(beta) == 1
    if (abs(cosbeta - 1.0d0) < EPS) then
       if (m1_in == m2_in) then
          d_wigner = 1.0d0
       else
          d_wigner = 0.0d0
       end if
       return
    end if

    ! (1)  d(l,m2,m1) = (-1)**(m2-m1) * d(l,m1,m2)
    pre = 1.0d0
    if (abs(m2_in) > abs(m1_in)) then
       m1 = m2_in
       m2 = m1_in
       if (mod(m2-m1, 2) /= 0) pre = -pre
    else
       m1 = m1_in
       m2 = m2_in
    end if

    ! (2)  d(l,-m1,-m2) = (-1)**(m2-m1) * d(l,m1,m2)
    if (m1 < 0) then
       m1 = -m1
       m2 = -m2
       if (mod(m2-m1, 2) /= 0) pre = -pre
    end if

    if ( (l < 0) .or. (m1 < 0) .or. (m1 > l) .or. (abs(m2) > abs(m1)) &
         .or. (abs(cosbeta) > 1.0d0) ) then
       write(0,*) "Error: invalid query in wigner_d table:"
       write(0,*) "       l,m,m'    = ", l, m1, m2
       write(0,*) "       cos(beta) = ", cosbeta
       stop
    end if

                   !  l m1 m2
    l_sel : select case(l)
    case(0)        !  0  0  0
       d_wigner = 1.0d0
    case(1)
       select case(m1)
       case(0)     !  1  0  0
          d_wigner = cosbeta
       case(1)
          select case(m2)
          case(0)  !  1  1  0
             ! note: beta in [0, pi] --> sin(beta) > 0
             d_wigner = -sqrt(0.5d0*(1.0d0 - cosbeta*cosbeta))
          case(-1) !  1  1 -1
             d_wigner = 0.5d0*(1.0d0 - cosbeta)
          case(1)  !  1  1  1
             d_wigner = 0.5d0*(1.0d0 + cosbeta)
          case default 
             d_wigner = 0.0d0
             write(0,*) l, m1, m2
             stop 999
          end select ! m2
       case default 
          d_wigner = 0.0d0
          write(0,*) l, m1, m2
          stop 998
       end select ! m1
    case(2)
       select case(m1)
       case(0)     !  2  0  0
          d_wigner = 1.5d0*cosbeta*cosbeta - 0.5d0
       case(1)
          select case(m2)
          case(0)  !  2  1  0
             d_wigner = -sqrt(1.5d0*(1.0d0 - cosbeta*cosbeta))*cosbeta
          case(-1) !  2  1 -1
             d_wigner = 0.5d0*(1.0d0 - cosbeta)*(2.0d0*cosbeta + 1.0d0)
          case(1)  !  2  1  1
             d_wigner = 0.5d0*(1.0d0 + cosbeta)*(2.0d0*cosbeta - 1.0d0)
          case default 
             d_wigner = 0.0d0
             write(0,*) l, m1, m2
             stop 997
          end select ! m2
       case(2)
          select case(m2)
          case(0)  !  2  2  0
             d_wigner = sqrt(0.375d0)*(1.0d0 - cosbeta*cosbeta)
          case(-1) !  2  2 -1
             d_wigner = -0.5d0*sqrt(1.0d0 - cosbeta*cosbeta)*(1.0d0 - cosbeta)
          case(1)  !  2  2  1
             d_wigner = -0.5d0*sqrt(1.0d0 - cosbeta*cosbeta)*(1.0d0 + cosbeta)
          case(-2) !  2  2 -2
             rtmp     = 0.5d0*(1.0d0 - cosbeta)
             d_wigner = rtmp*rtmp
          case(2)  !  2  2  2
             rtmp     = 0.5d0*(1.0d0 + cosbeta)
             d_wigner = rtmp*rtmp
          case default 
             d_wigner = 0.0d0
             write(0,*) l, m1, m2
             stop 996
          end select ! m2
       case default 
          d_wigner = 0.0d0
          write(0,*) l, m1, m2
          stop 995
       end select ! m1
    case default
       d_wigner = rot_wigner_d_func(l,m1,m2,cosbeta)
    end select l_sel

    if (abs(d_wigner) > EPS) then
       d_wigner = pre*d_wigner
    else       
       d_wigner = 0.0d0
    end if

  end function rot_wigner_d

  !--------------------------------------------------------------------!
  ! implementation of an analytic form of the Wigner small d function  !
  !--------------------------------------------------------------------!

  function rot_wigner_d_func(l, m1_in, m2_in, cosbeta) result(d_wigner)

    ! FIXME !!! 
    ! carefully check this function for l > 1 !!! 
    ! --> maybe better implement a recursion based function anyway:
    !     see e.g. Blanco et al., J. Mol. Struc. 419 (1997) 19-27
    ! FIXME !!!

    implicit none

    integer,          intent(in) :: l, m1_in, m2_in
    double precision, intent(in) :: cosbeta
    double precision             :: d_wigner

    double precision :: d_pre
    integer          :: itmp, m1, m2
    integer          :: k, k_i, k_f
    double precision :: cos_bo2
    double precision :: sin_bo2

    d_pre = 1.0d0
    if (m1_in < 0) then
       m1 = -m1_in
       m2 = -m2_in
       if (mod(m1-m2, 2) /= 0) d_pre = -d_pre
    else
       m1 = m1_in
       m2 = m2_in
    end if

    if ( (l<0) .or. (m1<0) .or. (m1>l) .or. (abs(m2)>l) &
         .or. (abs(cosbeta)>1.0d0) ) then
       write(0,*) 'Error: invalid arguments in rot_wigner_d'
       print *, 'l,m1,m2', l, m1, m2
       stop
    end if

    ! if cos(beta) == 1
    if (abs(cosbeta - 1.0d0) < EPS) then
       if (m1 == m2) then
          d_wigner = 1.0d0
       else
          d_wigner = 0.0d0
       end if
       return
    end if

    cos_bo2 = 0.5d0*(cosbeta + 1.0d0)
    sin_bo2 = sqrt(1.0d0 - cos_bo2)
    cos_bo2 = sqrt(cos_bo2)

    ! prefactor:
    itmp = rot_fact(l+m1)*rot_fact(l-m1)*rot_fact(l+m2)*rot_fact(l-m2)
    d_pre = d_pre*sqrt(dble(itmp))

    ! TODO:
    ! optimize the following piece of code!!!

    k_i = max(0, m2-m1)
    k_f = min(l-m1, l+m2)
    d_wigner = 0.0d0
    do k = k_i, k_f
       itmp = rot_fact(k)*rot_fact(l+m2-k)*rot_fact(m1-m2+k)*rot_fact(l-m1-k)
       if (mod(m1-m2+k,2)/=0) then
          itmp = -itmp
       end if

       d_wigner = d_wigner + ( cos_bo2**(2*l+m2-m1-2*k) &
                             * sin_bo2**(m1-m2+2*k)     &
                            ) / dble(itmp)

    end do
    if (k_f >= k_i) then
       d_wigner = d_pre * d_wigner
    else
       d_wigner = d_pre
    end if

  end function rot_wigner_d_func

  !--------------------------------------------------------------------!
  !                          3D visualization                          !
  !                                                                    !
  ! usage:    call plot3d_open(...)                                    !
  !           do                                                       !
  !              call plot3d_iter(x, y, z, stat)                       !
  !              call plot3d_write(f(x,y,z), stat)                     !
  !              if (stat == 1) exit                                   !
  !           end do                                                   !
  !                                                                    !
  ! Will create a cube file with a grid of function values of function !
  ! f.  The grid is specified when calling plot3d_open().              !
  !--------------------------------------------------------------------!

  subroutine plot3d_open(outfile, xmin, xmax, ymin, ymax, zmin, zmax, dr)

    implicit none

    character(len=*),   intent(in)  :: outfile
    double precision,   intent(in)  :: xmin, xmax, ymin, ymax, zmin, zmax
    double precision,   intent(in)  :: dr
    integer                         :: i
    double precision                :: x0, y0, z0

    if (p3d_open) then
       write(0,*) 'Error: plot stream already open in plot3d_open.'
       stop
    end if

    p3d_nx = ceiling((xmax-xmin)/dr) + 1
    p3d_ny = ceiling((ymax-ymin)/dr) + 1
    p3d_nz = ceiling((zmax-zmin)/dr) + 1

    x0 = 0.5d0*( xmax + xmin - p3d_nx*dr )
    y0 = 0.5d0*( ymax + ymin - p3d_ny*dr )
    z0 = 0.5d0*( zmax + zmin - p3d_nz*dr )

    open(p3d_unit, file=trim(adjustl(outfile)), action="write", &
         status="replace")

    ! write cube file header:
    write(p3d_unit, *) 'cube file with volumetric data'
    write(p3d_unit, *)
    write(p3d_unit, '(I5,3F12.6)') 1, x0, y0, z0
    write(p3d_unit, '(I5,3F12.6)') p3d_nx, dr, 0.0, 0.0
    write(p3d_unit, '(I5,3F12.6)') p3d_ny, 0.0, dr, 0.0
    write(p3d_unit, '(I5,3F12.6)') p3d_nz, 0.0, 0.0, dr
    write(p3d_unit, '(I5,4F12.6)') 1, 0.0, 0.0, 0.0, 0.0

    ! calculate x,y,z coordinates:
    allocate(p3d_cx(p3d_nx),p3d_cy(p3d_ny),p3d_cz(p3d_nz))
    p3d_cx(1) = xmin
    do i = 2, p3d_nx
       p3d_cx(i) = p3d_cx(i-1) + dr
    end do
    p3d_cy(1) = ymin
    do i = 2, p3d_ny
       p3d_cy(i) = p3d_cy(i-1) + dr
    end do
    p3d_cz(1) = xmin
    do i = 2, p3d_nz
       p3d_cz(i) = p3d_cz(i-1) + dr
    end do

    p3d_open = .true.
    p3d_ix   = 1
    p3d_iy   = 1
    p3d_iz   = 1

  end subroutine plot3d_open

  !--------------------------------------------------------------------!

  subroutine plot3d_iter(x, y, z, stat)
    ! stat = 0 --- continue
    !        1 --- finish
    !        2 --- insert line break

    implicit none

    double precision, intent(out) :: x, y, z
    integer,          intent(out) :: stat

    if (.not. p3d_open) then
       write(0,*) 'Error: no plot stream open in plot3d_iter.'
       stop
    end if

    stat = 0

    x = p3d_cx(p3d_ix)
    y = p3d_cy(p3d_iy)
    z = p3d_cz(p3d_iz)

    if (mod(p3d_iz,6) == 0) stat = 2

    p3d_iz = mod(p3d_iz, p3d_nz) + 1
    if (p3d_iz == 1) then
       stat = 2
       p3d_iy = mod(p3d_iy, p3d_ny) + 1
       if (p3d_iy == 1) then
          p3d_ix = mod(p3d_ix, p3d_nx) + 1
       end if
    end if

    if ((p3d_ix == 1) .and. (p3d_iy == 1) .and. (p3d_iz == 1)) then
       stat = 1
    end if

  end subroutine plot3d_iter

  !--------------------------------------------------------------------!

  subroutine plot3d_write(val, stat)

    implicit none

    double precision, intent(in) :: val
    integer,          intent(in) :: stat

    if (.not. p3d_open) then
       write(0,*) 'Error: no plot stream open in plot3d_iter.'
       stop
    end if

    if (stat == 1) then
       write(p3d_unit, '(E12.5,1x)') val
       deallocate(p3d_cx, p3d_cy, p3d_cz)
       p3d_open = .false.
    else
       write(p3d_unit, '(E12.5,1x)', advance='no') val
       if (stat == 2) write(p3d_unit, *)
    end if

  end subroutine plot3d_write




end module harmonics











!!$  subroutine rot_integral_old(alpha, cosbeta, l1, m1, l2, m2, n_m, c_m)
!!$
!!$    implicit none
!!$
!!$    double precision, parameter :: eps_here = 1.0d-6
!!$
!!$    double precision,                 intent(in)  :: alpha, cosbeta
!!$    integer,                          intent(in)  :: l1, l2, m1, m2
!!$    integer,                          intent(in)  :: n_m
!!$    double precision, dimension(n_m), intent(out) :: c_m
!!$
!!$    double precision :: Am1, Am2, Bm1, Bm2, S1, S2, T1, T2
!!$    double precision :: m_alpha, d_wigner
!!$    integer          :: mp
!!$
!!$    if ( (abs(m1)>l1) .or. (abs(m2)>l2) .or. (n_m /= min(l1,l2)+1) &
!!$       .or. (alpha<0.0d0) .or. (alpha>=pi2) .or. (abs(cosbeta)>1.0d0) ) then
!!$       write(0,*) 'Error: inadequate parameters in rot_integral.'
!!$       write(0,*) 'alpha   = ', alpha, alpha-pi2
!!$       write(0,*) 'cosbeta = ', cosbeta
!!$       write(0,*) 'l1, m1  = ', l1, m1
!!$       write(0,*) 'l2, m2  = ', l2, m2
!!$       write(0,*) 'n_m     = ', n_m
!!$       stop
!!$    end if
!!$
!!$    if (.not. rot_is_init) then
!!$       write(0,*) 'Error: rotation procedures not initialized!'
!!$       stop
!!$    end if
!!$
!!$    !-----------------------------------------------!
!!$    ! A and B factors for (m1,alpha) and (m2,alpha) !
!!$    !-----------------------------------------------!
!!$
!!$    if (abs(alpha) < eps_here) then
!!$
!!$       ! m1:
!!$       if (m1 == 0) then
!!$          Am1 = sqrt2_inv
!!$          Bm1 = 0.0d0
!!$       else if (m1 > 0) then
!!$          Am1      = 1.0d0
!!$          Bm1      = 0.0d0
!!$          if (mod(m1,2) /= 0) Am1 = -Am1
!!$       else ! m1 < 0
!!$          Am1      = 0.0d0
!!$          Bm1      = 1.0d0
!!$          if (mod(m1,2) /= 0) Bm1 = -Bm1
!!$       end if
!!$
!!$       ! m2:
!!$       if (m2 == 0) then
!!$          Am2 = sqrt2_inv
!!$          Bm2 = 0.0d0
!!$       else if (m2 > 0) then
!!$          Am2      = 1.0d0
!!$          Bm2      = 0.0d0
!!$          if (mod(m2,2) /= 0) Am2 = -Am2
!!$       else ! m2 < 0
!!$          Am2      = 0.0d0
!!$          Bm2      = 1.0d0
!!$          if (mod(m2,2) /= 0) Bm2 = -Bm2
!!$       end if
!!$
!!$    else ! if alpha > 0
!!$
!!$       ! m1:
!!$       if (m1 == 0) then
!!$          Am1 = sqrt2_inv
!!$          Bm1 = 0.0d0
!!$       else if (m1 > 0) then
!!$          m_alpha = dble(m1)*alpha
!!$          Am1      = cos(m_alpha)
!!$          Bm1      = sin(m_alpha)
!!$          if (mod(m1,2) /= 0) then
!!$             Am1 = -Am1
!!$             Bm1 = -Bm1
!!$          end if
!!$       else ! m1 < 0
!!$          m_alpha = dble(m1)*alpha
!!$          Am1      = sin(m_alpha)
!!$          Bm1      = cos(m_alpha)
!!$          if (mod(m1,2) /= 0) then
!!$             Am1 = -Am1
!!$             Bm1 = -Bm1
!!$          end if
!!$       end if
!!$
!!$       ! m2:
!!$       if (m2 == 0) then
!!$          Am2 = sqrt2_inv
!!$          Bm2 = 0.0d0
!!$       else if (m2 > 0) then
!!$          m_alpha = dble(m2)*alpha
!!$          Am2      = cos(m_alpha)
!!$          Bm2      = sin(m_alpha)
!!$          if (mod(m2,2) /= 0) then
!!$             Am2 = -Am2
!!$             Bm2 = -Bm2
!!$          end if
!!$       else ! m2 < 0
!!$          m_alpha = dble(m2)*alpha
!!$          Am2      = sin(m_alpha)
!!$          Bm2      = cos(m_alpha)
!!$          if (mod(m2,2) /= 0) then
!!$             Am2 = -Am2
!!$             Bm2 = -Bm2
!!$          end if
!!$       end if
!!$       
!!$    end if
!!$
!!$    !---------------------------!
!!$    ! coefficient for (l1,l2,0) !
!!$    !---------------------------!
!!$
!!$    c_m(1) = 2.0d0*Am1*Am2*rot_wigner_d(l1,abs(m1),0,cosbeta) &
!!$           * rot_wigner_d(l2,abs(m2),0,cosbeta)
!!$
!!$    ! loop over `m prime' mp:
!!$    mprime : do mp = 1, min(l1, l2)
!!$
!!$       !---------------------------------------------------!
!!$       ! S and T factors for (l1,m1,mp) and for (l2,m2,mp) !
!!$       !---------------------------------------------------!
!!$
!!$       ! l1, m1:
!!$       d_wigner = rot_wigner_d(l1,abs(m1),mp,cosbeta)
!!$       if (mod(mp,2) == 0) then
!!$          S1 = d_wigner
!!$       else
!!$          S1 = -d_wigner
!!$       end if
!!$       d_wigner = rot_wigner_d(l1,abs(m1),-mp,cosbeta)
!!$       if (m1 == 0) then
!!$          T1 = 0.0d0
!!$       else
!!$          T1 = Bm1 * (S1 - d_wigner)      
!!$       end if
!!$       S1 = Am1 * (S1 + d_wigner)
!!$
!!$       ! l2, m2:
!!$       d_wigner = rot_wigner_d(l2,abs(m2),mp,cosbeta)
!!$       if (mod(mp,2) == 0) then
!!$          S2 = d_wigner
!!$       else
!!$          S2 = -d_wigner
!!$       end if
!!$       d_wigner = rot_wigner_d(l2,abs(m2),-mp,cosbeta)
!!$       if (m2 == 0) then
!!$          T2 = 0.0d0
!!$       else
!!$          T2 = Bm2 * (S2 - d_wigner)      
!!$       end if
!!$       S2 = Am2 * (S2 + d_wigner)
!!$
!!$       !-------------------------------------------!
!!$       ! calculate coefficient for current m prime !
!!$       !-------------------------------------------!
!!$
!!$       c_m(mp+1) = (S1*S2 + T1*T2)    
!!$
!!$    end do mprime
!!$
!!$  end subroutine rot_integral_old
