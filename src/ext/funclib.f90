module funclib

  !--------------------------------------------------------------------!
  !      a library with implementations of mathemetical functions      !
  !                                                                    !
  ! currently implemented:                                             !
  !                                                                    !
  ! P_lgndr  -  associated Legendre polynomial                         !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! 2010-11-03 Alexander Urban                                         !
  !--------------------------------------------------------------------!

  implicit none

  double precision, parameter :: EPS   = epsilon(1.0d0)
  double precision, parameter :: FPMIN = tiny(1.0d0)
  double precision, parameter :: FPMAX = huge(1.0d0)

contains

  !--------------------------------------------------------------------!
  !                  associated Legendre polynomials                   !
  !                                                                    !
  !     plgndr - (C) Copr. 1986-92 Numerical Recipes Software ,4-#     !
  !                   (modified by Alexander Urban)                    !
  !--------------------------------------------------------------------!

  function P_lgndr(l,m,x) result(P_lm)

    implicit none

    integer,          intent(in) :: l, m
    double precision, intent(in) :: x
    double precision             :: P_lm
    integer                      :: i, ll
    double precision             :: fact, pll, pmm, pmmp1, somx2

    if ( (m < 0) .or. (m > l) .or. (abs(x) > 1.0d0) ) then
       write(0,*) 'Error: invalid parameters in P_lgndr'
       stop
    end if

    pmm = 1.0d0

    if ( m > 0 ) then
       somx2 = sqrt((1.0d0-x)*(1.0d0+x))
       fact  = 1.0d0
       do i = 1, m
          pmm  = -pmm*fact*somx2
          fact = fact+2.
       end do
    end if

    if ( l == m ) then
       P_lm = pmm
    else
       pmmp1 = x*(2*m + 1)*pmm
       if ( l == (m+1) ) then
          P_lm = pmmp1
       else
          pll = 0.0d0
          do ll = m+2, l
             pll   = (x*(2*ll - 1)*pmmp1 - (ll + m - 1)*pmm)/(ll - m)
             pmm   = pmmp1
             pmmp1 = pll
          end do
          P_lm = pll
       end if
    end if

  end function P_lgndr

  !--------------------------------------------------------------------!
  !                           gamma function                           !
  !                                                                    !
  ! G(z) = int_0^\infty { t^{z - 1} e^{-t} dt }                        !
  !                                                                    !
  ! G(z) becomes large very quickly, therefore gammln(z) returns the   !
  ! logarithm ln(G(z)).                                                !
  !                                                                    !
  ! adapted from: Press, et al. Numerical Recipes in F77, 1986-92      !
  !                                                                    !
  !--------------------------------------------------------------------!

  function gammln(xx) result(gln)

    implicit none

    double precision, intent(in)  :: xx
    double precision              :: gln

    double precision,               parameter :: stp = 2.5066282746310005d0
    double precision, dimension(6), parameter :: cof = &
    (/ 76.18009172947146d0, -86.50532032941677d0,    24.01409824083091d0, &
       -1.231739572450155d0,  0.1208650973866179d-2, -0.5395239384953d-5 /)

    integer          :: j
    double precision :: ser, tmp, x, y

    x   = xx
    y   = x
    tmp = x + 5.5d0
    tmp = (x + 0.5d0)*log(tmp) - tmp
    ser = 1.000000000190015d0
    do j = 1, 6
       y   = y + 1.d0
       ser = ser + cof(j)/y
    end do
    gln = tmp + log(stp*ser/x)

  end function gammln

  !--------------------------------------------------------------------!
  !                     incomplete gamma functions                     !
  !                                                                    !
  ! gammap - P(a,x) = 1/gamma(a) int_0^x { e^{it} t^{a-1} dt } (a > 0) !
  ! gammaq - Q(a,x) = 1 - P(a,x)                                       !
  !                                                                    !
  ! adapted from: Press, et al. Numerical Recipes in F77, 1986-92      !
  !                                                                    !
  !--------------------------------------------------------------------!

  function gamma_P(a,x) result(P)

    implicit none

    ! NR name: gammp

    double precision, intent(in) :: a, x
    double precision             :: P
    double precision             :: gammcf, gamser, gln

    if ((x < 0.0d0) .or. (a <= 0.0d0)) then
       write(0,*) "Error: bad arguments in `gamma_P()'."
       P = 0.0d0
       return
    end if
    
    if (x < a + 1.0d0) then
       call gamma_series(gamser,a,x,gln)
       P = gamser
    else
       call gamma_cf(gammcf,a,x,gln)
       P = 1.0d0 - gammcf
    endif

  end function gamma_P

  !--------------------------------------------------------------------!

  function gamma_Q(a,x) result(Q)

    implicit none

    ! NR name: gammq()

    double precision, intent(in) :: a, x
    double precision             :: Q
    double precision             :: gammcf, gamser, gln

    if ((x < 0.0d0) .or. (a <= 0.0d0)) then
       write(0,*) "Error: bad arguments in `gamma_Q()'."
       Q = 0.0d0
       return
    end if
    
    if (x < a + 1.0d0) then
       call gamma_series(gamser,a,x,gln)
       Q = 1.0d0 - gamser
    else
       call gamma_cf(gammcf,a,x,gln)
       Q = gammcf
    endif

  end function gamma_Q

  !--------------------------------------------------------------------!

  subroutine gamma_series(gamser, a, x, gln)

    implicit none

    ! NR name: gser()

    double precision, intent(in)  :: a, x
    double precision, intent(out) :: gamser, gln

    integer,            parameter :: ITMAX = 1000
    integer                       :: n
    double precision              :: ap, del, gsum

    gln = gammln(a)

    if (x <= 0.0d0) then
       if (x < 0.0d0) then
          write(0,*) "Error: x < 0 in `gamma_series()'"
       end if
       gamser = 0.0d0
       return
    endif

    ap  = a
    gsum = 1.0d0/a
    del = gsum
    
    n = 1
    do 
       ap  = ap + 1.0d0
       del = del*x/ap
       gsum = gsum + del
       if (abs(del) < abs(gsum)*EPS) exit
       n = n + 1
       if (n > ITMAX) then
          write(0,*) "Warning: a too large, ITMAX too small in `gamma_series()'"
          exit
       end if
    end do
    
    gamser = gsum*exp(-x + a*log(x) - gln)

  end subroutine gamma_series

  !--------------------------------------------------------------------!

  subroutine gamma_cf(gammcf, a, x, gln)

    implicit none

    ! NR name: gcf()

    double precision, intent(in)  :: a, x
    double precision, intent(out) :: gammcf, gln

    integer,            parameter :: ITMAX = 1000
    integer                       :: i
    double precision              :: an, b, c, d, del, h

    gln = gammln(a)
    b   = x + 1.0d0 - a
    c   = 1.0d0/FPMIN
    d   = 1.0d0/b
    h   = d

    i = 1
    do
       an  = -dble(i)*(dble(i) - a)
       b   = b + 2.0d0
       d   = an*d + b
       if (abs(d) < FPMIN) d = FPMIN
       c   = b + an/c
       if (abs(c) < FPMIN) c = FPMIN
       d   = 1.0d0/d
       del = d*c
       h   = h*del
       if (abs(del - 1.0d0) < EPS) exit
       i = i + 1
       if (i > ITMAX) then
          write(0,*) "Warning: a too large, ITMAX too small in `gamma_cf()'"
          exit
       end if
    end do

    gammcf = exp(-x + a*log(x) - gln)*h

  end subroutine gamma_cf

  !--------------------------------------------------------------------!
  !          error function and complementary error function           !
  !                                                                    !
  ! erf(x)  = 2/sqrt(pi) int_0^x { e^{-t^2} dt }                       !
  ! erfc(x) = 1 - erf(x) = 2/sqrt(pi) int_0^\infty { e^{-t^2} dt }     !
  !                                                                    !
  ! relation to the incomplete gamma functions (see above):            !
  !                                                                    !
  !   erf(x)  = P(1/2, x^2) ;  (x >= 0)                                !
  !   erfc(x) = Q(1/2, x^2) ;  (x >= 0)                                !
  !                                                                    !
  ! adapted from: Press, et al. Numerical Recipes in F77, 1986-92      !
  !                                                                    !
  !--------------------------------------------------------------------!

  function erf(x) result(y)

    implicit none

    double precision, intent(in) :: x
    double precision             :: y

    if (x < 0.0d0) then
       y = -gamma_P(0.5d0, x*x)
    else
       y =  gamma_P(0.5d0, x*x)
    endif
    
  end function erf

  !--------------------------------------------------------------------!

  function erfc(x) result(y)

    implicit none

    double precision, intent(in) :: x
    double precision             :: y

    if (x < 0.0d0) then
       y = 1.0d0 + gamma_P(0.5d0, x*x)
    else
       y = gamma_Q(0.5d0, x*x)
    endif
    
  end function erfc


end module funclib
