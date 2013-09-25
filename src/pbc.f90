module pbc

  !--------------------------------------------------------------------!
  ! This module provides a number of procedures for operations on a    !
  ! periodic grid/lattice.  No data is stored, so the module does not  !
  ! need to be initialized.                                            !
  !--------------------------------------------------------------------!
  ! 2010-10-29 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use sortlib,   only: argsort

  implicit none

  public  :: pbc_MP_kpoint_set,   &
             pbc_d_max_origin,    &
             pbc_number_of_tvecs, &
             pbc_compute_tvecs

contains

  !--------------------------------------------------------------------!
  !                  Monkhorst/Pack like k-point set                   !
  !                                                                    !
  ! This is not the original MP k-point set but a shifted one, so that !
  ! for (s1, s2, s3) = (0, 0, 0) the Gamma point is included.          !
  !--------------------------------------------------------------------!

  subroutine pbc_MP_kpoint_set(n, s, nkpts, kpt, w)

    implicit none

    integer,          dimension(3),         intent(in)  :: n
    double precision, dimension(3),         intent(in)  :: s
    integer,                                intent(in)  :: nkpts
    double precision, dimension(1:3,nkpts), intent(out) :: kpt
    double precision, dimension(nkpts),     intent(out) :: w
    
    integer          :: i1, i2, i3, ikpt
    double precision :: n1_inv, n2_inv, n3_inv

    if (nkpts /= n(1)*n(2)*n(3)) then
       write(0,*) "Error: ill-sized k-point array in `pbc_MP_kpoint_set()'."
       stop
    end if

    n1_inv = 1.0d0/dble(n(1))
    n2_inv = 1.0d0/dble(n(2))
    n3_inv = 1.0d0/dble(n(3))

    ikpt = 0
    do i1 = 0, (n(1) - 1)
    do i2 = 0, (n(2) - 1)
    do i3 = 0, (n(3) - 1)
       ikpt = ikpt + 1
       kpt(1,ikpt) = n1_inv*(dble(i1) + s(1))
       kpt(2,ikpt) = n2_inv*(dble(i2) + s(2))
       kpt(3,ikpt) = n3_inv*(dble(i3) + s(3))
    end do
    end do
    end do

    ! for now, all k-points are equally weighted:
    w(1:nkpts) = 1.0d0/dble(nkpts)

  end subroutine pbc_MP_kpoint_set

  !--------------------------------------------------------------------!
  !              max. distance of an atom from the origin              !
  !--------------------------------------------------------------------!

  function pbc_d_max_origin(natoms, coo, avec) result(d_max)

    implicit none

    integer,                                intent(in) :: natoms
    double precision, dimension(3, natoms), intent(in) :: coo
    double precision, dimension(3,3),       intent(in) :: avec
    double precision                                   :: d_max
    
    integer                        :: iat
    double precision, dimension(3) :: cart

    d_max = 0.0d0
    do iat = 1, natoms
       cart(1:3) = coo(1,iat)*avec(1:3,1) &
                 + coo(2,iat)*avec(1:3,2) &
                 + coo(3,iat)*avec(1:3,3)
       d_max = max(d_max, sum(cart*cart))
    end do
    d_max = sqrt(d_max)

  end function pbc_d_max_origin

  !--------------------------------------------------------------------!
  !          translation vectors within given cut-off radius           !
  !                                                                    !
  ! The T vectors are returned sorted by increasing length.            !
  ! Note: T = (0,0,0) is not considered by these subroutines.          !
  !--------------------------------------------------------------------!

  function pbc_number_of_tvecs(r_cut, d_max, avec) result(ntvecs)

    implicit none

    !------------------------------------------------------------------!
    ! r_cut    : radial cut-off                                        !
    ! d_max    : max. distance of a point/atom from the origin         !
    ! avec(i,j): i-th component of the j-th lattice vector             !
    !------------------------------------------------------------------!

    double precision,                 intent(in) :: r_cut
    double precision,                 intent(in) :: d_max
    double precision, dimension(3,3), intent(in) :: avec
    integer                                      :: ntvecs

    integer                        :: ia, ib, ic
    logical                        :: lasta, lastb, lastc
    double precision, dimension(3) :: cart
    double precision               :: r2, rcut2

    rcut2 = r_cut + 2.0d0*d_max
    rcut2 = rcut2*rcut2

    ntvecs = 0
    
    ic = 0
    c_do : do
       lastc = .false.
       ib = 0
       b_do : do
          lastb = .false.
          ia = 0
          a_do : do
             lasta = .false.
             if ((ic==0) .and. (ib==0) .and. (ia==0)) then
                ia = ia + 1
                cycle a_do
             end if
             
             ! (+ia, +ib, +ic)
             cart(1:3) = dble(ia)*avec(1:3,1) &
                       + dble(ib)*avec(1:3,2) &
                       + dble(ic)*avec(1:3,3) 
             r2 = sum(cart*cart)
             if (r2 <= rcut2) then
                lasta = .true.
                ntvecs = ntvecs + 1
             end if

             if ((ia /= 0) .and. (ib /= 0)) then

                ! (+ia, -ib, +ic)
                cart(1:3) = dble(ia)*avec(1:3,1) &
                          + dble(-ib)*avec(1:3,2) &
                          + dble(ic)*avec(1:3,3) 
                r2 = sum(cart*cart)
                if (r2 <= rcut2) then
                   lasta = .true.
                   ntvecs = ntvecs + 1
                end if

                if (ic /= 0) then

                   ! (-ia, +ib, +ic)
                   cart(1:3) = dble(-ia)*avec(1:3,1) &
                             + dble(ib)*avec(1:3,2) &
                             + dble(ic)*avec(1:3,3) 
                   r2 = sum(cart*cart)
                   if (r2 <= rcut2) then
                      lasta = .true.
                      ntvecs = ntvecs + 1
                   end if
                   
                   ! (-ia, -ib, +ic)
                   cart(1:3) = dble(-ia)*avec(1:3,1) &
                             + dble(-ib)*avec(1:3,2) &
                             + dble(ic)*avec(1:3,3) 
                   r2 = sum(cart*cart)
                   if (r2 <= rcut2) then
                      lasta = .true.
                      ntvecs = ntvecs + 1
                   end if

                end if ! c /= 0

             else if ((ia /= 0) .and. (ic /= 0)) then

                   ! (-ia, +ib, +ic)
                   cart(1:3) = dble(-ia)*avec(1:3,1) &
                             + dble(ib)*avec(1:3,2) &
                             + dble(ic)*avec(1:3,3) 
                   r2 = sum(cart*cart)
                   if (r2 <= rcut2) then
                      lasta = .true.
                      ntvecs = ntvecs + 1
                   end if
                   
             else if ((ib /= 0) .and. (ic /= 0)) then

                ! (+ia, -ib, +ic)
                cart(1:3) = dble(ia)*avec(1:3,1) &
                          + dble(-ib)*avec(1:3,2) &
                          + dble(ic)*avec(1:3,3) 
                r2 = sum(cart*cart)
                if (r2 <= rcut2) then
                   lasta = .true.
                   ntvecs = ntvecs + 1
                end if

             end if

             if (.not. lasta) exit a_do
             lastb = .true.
             ia = ia + 1
          end do a_do
          if (.not. lastb) exit b_do
          lastc = .true.
          ib = ib + 1
       end do b_do
       if (.not. lastc) exit c_do
       ic = ic + 1
    end do c_do
  
  end function pbc_number_of_tvecs

  !--------------------------------------------------------------------!
  ! Note: the norm of the T vectors is returned squared T*T            !
  !--------------------------------------------------------------------!

  subroutine pbc_compute_tvecs(r_cut, d_max, avec, ntvecs, tvec, tveclen)

    implicit none

    !------------------------------------------------------------------!
    ! r_cut      : radial cut-off                                      !
    ! d_max      : max. distance of a point/atom from the origin       !
    ! avec(i,j)  : i-th component of the j-th lattice vector           !
    ! ntvecs     : number of T vectors                                 !
    ! tvec(i,j)  : i-th component of the j-th T vector (on exit)       !
    ! tveclen(i) : norm of the i-th T vector (on exit)                 !
    !------------------------------------------------------------------!

    double precision,                      intent(in)  :: r_cut
    double precision,                      intent(in)  :: d_max
    double precision, dimension(3,3),      intent(in)  :: avec
    integer,                               intent(in)  :: ntvecs
    integer,          dimension(3,ntvecs), intent(out) :: tvec
    double precision, dimension(ntvecs),   intent(out) :: tveclen

    integer,          dimension(3,ntvecs)              :: buff_tvec
    double precision, dimension(ntvecs)                :: buff_tveclen
    integer,          dimension(ntvecs)                :: idx

    integer                                            :: ivec
    integer                                            :: ia, ib, ic
    logical                                            :: lasta, lastb, lastc
    double precision, dimension(3)                     :: cart
    double precision                                   :: r2, rcut2

    rcut2 = r_cut + 2.0d0*d_max
    rcut2 = rcut2*rcut2

    ivec = 0
    
    ic = 0
    c_do : do
       lastc = .false.
       ib = 0
       b_do : do
          lastb = .false.
          ia = 0
          a_do : do
             lasta = .false.
             if ((ic == 0) .and. (ib == 0) .and. (ia == 0)) then
                ia = ia + 1
                cycle a_do
             end if
             
             ! (+ia, +ib, +ic)
             cart(1:3) = dble(ia)*avec(1:3,1) &
                       + dble(ib)*avec(1:3,2) &
                       + dble(ic)*avec(1:3,3) 
             r2 = sum(cart*cart)
             if (r2 <= rcut2) then
                lasta = .true.
                ivec = ivec + 1
                buff_tvec(1:3,ivec) = (/ia,ib,ic/)
                buff_tveclen(ivec)  = r2
             end if

             if ((ia /= 0) .and. (ib /= 0)) then

                ! (+ia, -ib, +ic)
                cart(1:3) = dble(ia)*avec(1:3,1)  &
                          + dble(-ib)*avec(1:3,2) &
                          + dble(ic)*avec(1:3,3) 
                r2 = sum(cart*cart)
                if (r2 <= rcut2) then
                   lasta = .true.
                   ivec = ivec + 1
                   buff_tvec(1:3,ivec) = (/ia,-ib,ic/)
                   buff_tveclen(ivec)   = r2
                end if

                if (ic /= 0) then

                   ! (-ia, +ib, +ic)
                   cart(1:3) = dble(-ia)*avec(1:3,1) &
                             + dble(ib)*avec(1:3,2) &
                             + dble(ic)*avec(1:3,3) 
                   r2 = sum(cart*cart)
                   if (r2 <= rcut2) then
                      lasta = .true.
                      ivec = ivec + 1
                      buff_tvec(1:3,ivec) = (/-ia,ib,ic/)
                      buff_tveclen(ivec)  = r2
                   end if
                   
                   ! (-ia, -ib, +ic)
                   cart(1:3) = dble(-ia)*avec(1:3,1) &
                             + dble(-ib)*avec(1:3,2) &
                             + dble(ic)*avec(1:3,3) 
                   r2 = sum(cart*cart)
                   if (r2 <= rcut2) then
                      lasta = .true.
                      ivec = ivec + 1
                      buff_tvec(1:3,ivec) = (/-ia,-ib,ic/)
                      buff_tveclen(ivec)   = r2
                   end if

                end if ! c /= 0

             else if ((ia /= 0) .and. (ic /= 0)) then

                   ! (-ia, +ib, +ic)
                   cart(1:3) = dble(-ia)*avec(1:3,1) &
                             + dble(ib)*avec(1:3,2) &
                             + dble(ic)*avec(1:3,3) 
                   r2 = sum(cart*cart)
                   if (r2 <= rcut2) then
                      lasta = .true.
                      ivec = ivec + 1
                      buff_tvec(1:3,ivec) = (/-ia,ib,ic/)
                      buff_tveclen(ivec)  = r2
                   end if
                   
             else if ((ib /= 0) .and. (ic /= 0)) then

                ! (+ia, -ib, +ic)
                cart(1:3) = dble(ia)*avec(1:3,1) &
                          + dble(-ib)*avec(1:3,2) &
                          + dble(ic)*avec(1:3,3) 
                r2 = sum(cart*cart)
                if (r2 <= rcut2) then
                   lasta = .true.
                   ivec = ivec + 1
                   buff_tvec(1:3,ivec) = (/ia,-ib,ic/)
                   buff_tveclen(ivec)  = r2
                end if

             end if

             if (.not. lasta) exit a_do
             lastb = .true.
             ia = ia + 1
          end do a_do
          if (.not. lastb) exit b_do
          lastc = .true.
          ib = ib + 1
       end do b_do
       if (.not. lastc) exit c_do
       ic = ic + 1
    end do c_do

    if (ivec /= ntvecs) then
       write(0,*) "Error: inconsistent number of T vectors !"
       stop
    end if

    call argsort(buff_tveclen(1:ntvecs),idx)
    do ivec = 1, ntvecs
       tvec(1:3,ivec) = buff_tvec(1:3,idx(ivec))
       tveclen(ivec)  = sqrt(buff_tveclen(idx(ivec)))
    end do
   
  end subroutine pbc_compute_tvecs


end module pbc
