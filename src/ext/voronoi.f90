module voronoi

  use sortlib, only: argsort

  implicit none

  public  :: voro_wsc, convex_hull, circumcircle, circumsphere
  private :: voro_wscpoints, det3

  integer,                                       private :: npoints
  double precision, dimension(:,:), allocatable, private :: point
  double precision, dimension(:),   allocatable, private :: dist

contains

  !--------------------------------------------------------------------!
  !                                                                    !
  !                         Wigner-Seitz cell                          !
  !                                                                    !
  ! The following procedures implement the computation of the Wigner-  !
  ! Seitz cell in three dimensions.  In general they can be used to    !
  ! compute a single voronoi cell around the origin of a periodic grid !
  ! defined by three lattice/bais vectors (e.g. also Brillouin zone).  !
  !--------------------------------------------------------------------!

  subroutine voro_wsc(avec, wscpnt, npnts)

    implicit none

    integer,                              intent(inout) :: npnts
    double precision, dimension(3,3),     intent(in)    :: avec
    double precision, dimension(3,npnts), intent(out)   :: wscpnt

    double precision, dimension(3), parameter :: O = (/0.0d0, 0.0d0, 0.0d0/)
    double precision :: x0, y0, z0, r, r2, dx, dy, dz, dist2
    integer          :: ipt1, ipt2, ipt3, ipt4
    logical          :: accept
    double precision :: rmax

    if (.not. allocated(point)) then
       call voro_wscpoints(avec)
    end if

    if (npnts < npoints) then
       npnts = -npoints
       return
    end if

    rmax = sqrt(dist(npoints))

    npnts = 0
    do ipt1 = 1, npoints
    do ipt2 = ipt1+1, npoints
    do ipt3 = ipt2+1, npoints

       call circumsphere( point(1:3,ipt3), &
                          point(1:3,ipt2), &
                          point(1:3,ipt1), &
                          O, x0, y0, z0, r )
       if ((r < 0.0d0) .or. (r > rmax)) cycle

       ! check if we already have this point:
       accept = .true.
       check1 : do ipt4 = 1, npnts
          dx = wscpnt(1,ipt4) - x0
          dy = wscpnt(2,ipt4) - y0
          dz = wscpnt(3,ipt4) - z0
          dist2 = dx*dx + dy*dy + dz*dz
          if (dist2 < 1.0d-10) then
             accept = .false.
             exit check1
          end if
       end do check1
       if (.not. accept) cycle

       ! check if other sites are within the sphere:
       r2 = 0.99999d0*r*r
       check2 : do ipt4 = 1, npoints
          dx = point(1,ipt4) - x0
          dy = point(2,ipt4) - y0
          dz = point(3,ipt4) - z0
          dist2 = dx*dx + dy*dy + dz*dz
          if (dist2 < r2) then
             accept = .false.
             exit check2
          end if
       end do check2

       if (accept) then
          npnts = npnts + 1
          wscpnt(1, npnts) = x0
          wscpnt(2, npnts) = y0
          wscpnt(3, npnts) = z0
       end if

    end do
    end do
    end do

    if (allocated(point)) then
       deallocate(point, dist)
    end if

  end subroutine voro_wsc

  !--------------------------------------------------------------------!
  ! determine all points that might be necessary to compute the        !
  ! Wigner-Seitz cell of the grid defined by the grid vectors avec     !
  !--------------------------------------------------------------------!

  subroutine voro_wscpoints(avec)

    implicit none

    !------------------------------------------------------------------!
    ! avec(i,j)  : i-th component of the j-th lattice vector           !
    !------------------------------------------------------------------!

    double precision, dimension(3,3),           intent(in) :: avec

    integer,                                    parameter  :: memunit = 256

    integer                                                :: bufsize
    double precision, dimension(:,:), allocatable, target  :: buf1, buf2
    double precision, dimension(:,:),              pointer :: buffer
    integer,          dimension(:),   allocatable          :: idx1
    double precision, dimension(3)                         :: vec1
    double precision                                       :: r2, short, long
    integer                                                :: ia, ib, ic, ip
    logical                                                :: lasta, lastb, lastc

    !-----------------!
    ! compute cut-off !
    !-----------------!

    vec1(:) = avec(:,1) + avec(:,2) + avec(:,3)
    short = sum(vec1*vec1)
    long  = short
    do ia = -1, 1
       do ib = -1, 1
          do ic = -1, 1

             if ((ia == 0) .and. (ib == 0) .and. (ic == 0)) cycle

             vec1(:) = dble(ia)*avec(1:3,1) &
                  + dble(ib)*avec(1:3,2) &
                  + dble(ic)*avec(1:3,3)
             r2 = sum(vec1*vec1)
             if (r2 < short) then
                short = r2
             else if (r2 > long) then
                long = r2
             end if

          end do
       end do
    end do

    !-------------!
    ! open buffer !
    !-------------!

    bufsize = memunit
    allocate(buf1(4,bufsize))
    buffer => buf1

    !------------------------------------!
    ! look for points within the cut-off !
    !------------------------------------!

    npoints = 0
    ia = 0
    iado : do 
       lasta = .false.
       ib = 0
       ibdo : do
          lastb = .false.
          ic = 0
          icdo : do
             lastc = .false.

             if ((ia==0) .and. (ib==0) .and. (ic==0)) then
                ic = ic + 1
                cycle
             end if

             if (bufsize <= (npoints + 8)) then
                bufsize = bufsize + memunit
                if (allocated(buf1)) then
                   allocate(buf2(4,bufsize))
                   buf2(:,1:npoints) =  buf1(:,1:npoints)
                   buffer => buf2
                   deallocate(buf1)
                else
                   allocate(buf1(4,bufsize))
                   buf1(:,1:npoints) =  buf2(:,1:npoints)
                   buffer => buf1
                   deallocate(buf2)
                end if
             end if

             vec1(:) = dble(ia)*avec(1:3,1) &
                     + dble(ib)*avec(1:3,2) &
                     + dble(ic)*avec(1:3,3)
             r2 = sum(vec1*vec1)
             if (r2 <= long) then
                lastc = .true.
                npoints = npoints + 1
                buffer(1:3,npoints) = vec1(1:3)
                buffer(4,npoints)   = r2
                npoints = npoints + 1
                buffer(1:3,npoints) = -vec1(1:3)
                buffer(4,npoints)   = r2
             end if

             if (ia /= 0) then
                vec1(:) = dble(-ia)*avec(1:3,1) &
                        + dble(ib)*avec(1:3,2) &
                        + dble(ic)*avec(1:3,3)
                r2 = sum(vec1*vec1)
                if (r2 <= long) then
                   lastc = .true.
                   npoints = npoints + 1
                   buffer(1:3,npoints) = vec1(1:3)
                   buffer(4,npoints)   = r2
                   npoints = npoints + 1
                   buffer(1:3,npoints) = -vec1(1:3)
                   buffer(4,npoints)   = r2
                end if
             end if

             if (ib /= 0) then
                vec1(:) = dble(ia)*avec(1:3,1) &
                        + dble(-ib)*avec(1:3,2) &
                        + dble(ic)*avec(1:3,3)
                r2 = sum(vec1*vec1)
                if (r2 <= long) then
                   lastc = .true.
                   npoints = npoints + 1
                   buffer(1:3,npoints) = vec1(1:3)
                   buffer(4,npoints)   = r2
                   npoints = npoints + 1
                   buffer(1:3,npoints) = -vec1(1:3)
                   buffer(4,npoints)   = r2
                end if
             end if

             if (ic /= 0) then
                vec1(:) = dble(ia)*avec(1:3,1) &
                        + dble(ib)*avec(1:3,2) &
                        + dble(-ic)*avec(1:3,3)
                r2 = sum(vec1*vec1)
                if (r2 <= long) then
                   lastc = .true.
                   npoints = npoints + 1
                   buffer(1:3,npoints) = vec1(1:3)
                   buffer(4,npoints)   = r2
                   npoints = npoints + 1
                   buffer(1:3,npoints) = -vec1(1:3)
                   buffer(4,npoints)   = r2
                end if
             end if

             if (.not. lastc) exit icdo
             lastb = .true.
             ic = ic + 1
          end do icdo
          if (.not. lastb) exit ibdo
          lasta = .true.
          ib = ib + 1
       end do ibdo
       if (.not. lasta) exit iado
       ia = ia + 1
    end do iado

    !-----------------------------------------!
    ! sort points by distance from the origin !
    !-----------------------------------------!

    allocate(point(3,npoints), dist(npoints), idx1(npoints))
    call argsort(buffer(4,1:npoints),idx1)
    do ip = 1, npoints
       point(1:3,ip) = buffer(1:3,idx1(ip))
       dist(ip)      = buffer(4,idx1(ip))
    end do
    if (allocated(buf1)) then
       deallocate(buf1, idx1)
    else
       deallocate(buf2, idx1)
    end if

  end subroutine voro_wscpoints


  !====================================================================!
  !                                                                    !
  !                       Auxilliary procedures                        !
  !                                                                    !
  ! These procedures are used in some of the above routines.  Since    !
  ! they are quite general and might be useful for other applications  !
  ! the following routines are public (at least the main once).        !
  !                                                                    !
  !====================================================================!



  !--------------------------------------------------------------------!
  !                         convex hull in 3d                          !
  !                                                                    !
  ! All nodes of the voronoi cell are points on the convex hull.  This !
  ! general routine also works if points in the node-set are inside of !
  ! the convex hull.                                                   !
  !                                                                    !
  ! Note: this implementation is absolutely inefficient and should be  !
  !       used for small point sets only !                             !
  !--------------------------------------------------------------------!

  subroutine convex_hull(point, npoints, edge, nedges)

    implicit none

    !------------------------------------------------------------------!
    ! point(1:3,i)       : i-th point in the point set                 !
    ! npoints            : total number of points                      !
    ! edge(i,j)          : i-th vertex (i=1,2) of the j-th edge        !
    ! nedges             : in : size of the edge array                 !
    !                      out: total number of edges                  !
    !                                                                  !
    ! If the size of the edge array is unsufficient, the optimal size  !
    ! will be returned as -nedges.                                     !
    !                                                                  !
    !------------------------------------------------------------------!

    integer,                                intent(in)    :: npoints
    integer,                                intent(inout) :: nedges
    double precision, dimension(3,npoints), intent(in)    :: point
    integer,          dimension(2,nedges),  intent(out)   :: edge

    logical, dimension(npoints,npoints)  :: connect
    integer                              :: ncon

    integer                              :: ip1, ip2, ip3, ip4
    double precision, dimension(3)       :: r1, r2, r3, r4
    double precision                     :: V, V0, c
    double precision, dimension(3)       :: r1_x_r2, r1_x_r3
    logical                              :: valid
    integer, dimension(:,:), allocatable :: deledge
    integer                              :: i, ndel

    connect(:,:) = .false.
    ncon         = 0

    ! determine all edges on the convex hull:
    do ip1 = 1, npoints
       do ip2 = ip1+1, npoints
          do ip3 = ip2+1, npoints

             r1(:) = point(:,ip2) - point(:,ip1)
             r2(:) = point(:,ip3) - point(:,ip1)

             V0 = 0.0d0

             valid = .true.
             do ip4 = 1, npoints
                if (any(ip4 == (/ip1,ip2,ip3/))) cycle
                r4 = point(:,ip4) - point(:,ip1)
                V  = det3(r1,r2,r4)
                if (V0 == 0.0d0) V0 = V
                if (V*V0 < 0.0d0) then
                   valid = .false.
                   exit
                end if
             end do

             if (valid) then
                if (.not. connect(ip1, ip2)) then
                   connect(ip1, ip2) = .true.
                   connect(ip2, ip1) = .true.
                   ncon = ncon + 1
                end if
                if (.not. connect(ip1, ip3)) then
                   connect(ip1, ip3) = .true.
                   connect(ip3, ip1) = .true.
                   ncon = ncon + 1
                end if
                if (.not. connect(ip2, ip3)) then
                   connect(ip2, ip3) = .true.
                   connect(ip3, ip2) = .true.
                   ncon = ncon + 1
                end if
             end if

          end do
       end do
    end do

    ! remove redundant edges:
    allocate(deledge(2,ncon))
    ndel = 0
    do ip1 = 1, npoints
       do ip2 = ip1+1, npoints
          if (.not. connect(ip1,ip2)) cycle
          P3 : do ip3 = 1, npoints
             if (any(ip3 == (/ip1, ip2/))) cycle
             if (.not. connect(ip1,ip2)) exit
             if (.not. ( connect(ip1,ip3) &
                 .and.   connect(ip2,ip3) ) ) cycle

             do ip4 = 1, npoints
                if (any(ip4 == (/ip1, ip2, ip3/))) cycle
                if (.not. ( connect(ip1,ip4) &
                    .and.   connect(ip2,ip4) ) ) cycle

                ! At this point we have to triangles (p1,p2,p3) and
                ! (p1,p2,p4). Now check, if p3 and p4 lie in one plane:

                r1(:) = point(:,ip2) - point(:,ip1)
                r2(:) = point(:,ip3) - point(:,ip1)
                r3(:) = point(:,ip4) - point(:,ip1)
                V = det3(r1,r2,r3)
                if (abs(V) < 1.0d-10) then

                   ! now check, if p3 and p4 lie on different sides of
                   ! the (p1,p2) vector:

                   r1_x_r2(1) = r1(2)*r2(3) - r1(3)*r2(2)
                   r1_x_r2(2) = r1(3)*r2(1) - r1(1)*r2(3)
                   r1_x_r2(3) = r1(1)*r2(2) - r1(2)*r2(1)

                   r1_x_r3(1) = r1(2)*r3(3) - r1(3)*r3(2)
                   r1_x_r3(2) = r1(3)*r3(1) - r1(1)*r3(3)
                   r1_x_r3(3) = r1(1)*r3(2) - r1(2)*r3(1)

                   if (r1_x_r2(1) /= 0.0d0) then
                      c = r1_x_r2(1)/r1_x_r3(1)
                   else if (r1_x_r2(2) /= 0.0d0) then
                      c = r1_x_r2(2)/r1_x_r3(2)
                   else if (r1_x_r2(3) /= 0.0d0) then
                      c = r1_x_r2(3)/r1_x_r3(3)
                   else
                      write(0,*) "Warning: check point set in `convex_hull()'."
                      c = 0.0d0
                   end if

                   if (c < 0.0d0) then
                      ndel            = ndel + 1
                      deledge(1,ndel) = ip1
                      deledge(2,ndel) = ip2
                      exit P3
                   end if
                end if
             end do

          end do P3
       end do
    end do
    do i = 1, ndel
       if (connect(deledge(1,i),deledge(2,i))) then
          connect(deledge(1,i),deledge(2,i)) = .false.
          connect(deledge(2,i),deledge(1,i)) = .false.
          ncon = ncon - 1
       end if
    end do
    deallocate(deledge)

    ! return results:
    if (nedges < ncon) then
       nedges = -ncon
    else
       nedges = ncon
       ncon   = 0
       do ip1 = 1, npoints
       do ip2 = ip1+1, npoints
          if (connect(ip1,ip2)) then
             ncon = ncon + 1
             edge(1,ncon) = ip1
             edge(2,ncon) = ip2
          end if
       end do
       end do       
    end if

  end subroutine convex_hull

  !--------------------------------------------------------------------!

  function det3(P1, P2, P3) result(d)

    implicit none

    double precision, dimension(3), intent(in) :: P1, P2, P3
    double precision :: d

    d = P1(1)*P2(2)*P3(3) + P1(3)*P2(1)*P3(2) + P1(2)*P2(3)*P3(1) &
      - P1(3)*P2(2)*P3(1) - P1(1)*P2(3)*P3(2) - P1(2)*P2(1)*P3(3)

  end function det3

  !--------------------------------------------------------------------!
  !                    circumcircle / circumsphere                     !
  !                                                                    !
  ! The nodes of the voronoi cell are the centers of circumcircles in  !
  ! 2d, and the centers of circumspheres in 3d, of adjacent points.    !
  !--------------------------------------------------------------------!

  subroutine circumcircle(p1, p2, p3, x0, y0, r)

    implicit none

    double precision, dimension(2), intent(in)  :: p1, p2, p3
    double precision,               intent(out) :: x0, y0, r

    double precision :: xy1, xy2, xy3, d21, d31, d32, dx, dy
    double precision :: f

    f  = p1(2)*(p3(1)-p2(1)) + p2(2)*(p1(1)-p3(1)) + p3(2)*(p2(1)-p1(1))
    if (abs(f) <= 1.0d-15) then
       ! the points don't lie on a circle (since they are all on one line)
       ! --> return negative radius
       r = -1.0d0
       return
    end if

    f  = 0.5d0 / f

    xy1 = p1(1)*p1(1) + p1(2)*p1(2)
    xy2 = p2(1)*p2(1) + p2(2)*p2(2)
    xy3 = p3(1)*p3(1) + p3(2)*p3(2)

    d21 = xy2 - xy1
    d31 = xy3 - xy1
    d32 = xy3 - xy2

    x0 =  f*(p1(2)*d32 - p2(2)*d31 + p3(2)*d21)
    y0 = -f*(p1(1)*d32 - p2(1)*d31 + p3(1)*d21)

    dx = p1(1) - x0
    dy = p1(2) - y0
    r  = sqrt(dx*dx + dy*dy)

  end subroutine circumcircle

  !--------------------------------------------------------------------!

  subroutine circumsphere(p1,p2,p3,p4,x0,y0,z0,r)

    implicit none

    double precision, dimension(3), intent(in)  :: p1, p2, p3, p4
    double precision,               intent(out) :: x0, y0, z0, r

    double precision :: xyz1, xyz2, xyz3, xyz4
    double precision :: d21, d31, d41, d32, d42, d43
    double precision :: dx, dy, dz
    double precision :: f

    f = p1(3)*(p2(2)*(p4(1)-p3(1)) - p3(2)*(p4(1)-p2(1)) + p4(2)*(p3(1)-p2(1))) &
      - p2(3)*(p1(2)*(p4(1)-p3(1)) - p3(2)*(p4(1)-p1(1)) + p4(2)*(p3(1)-p1(1))) &
      + p3(3)*(p1(2)*(p4(1)-p2(1)) - p2(2)*(p4(1)-p1(1)) + p4(2)*(p2(1)-p1(1))) &
      - p4(3)*(p1(2)*(p3(1)-p2(1)) - p2(2)*(p3(1)-p1(1)) + p3(2)*(p2(1)-p1(1)))

    if (abs(f) <= 1.0d-15) then
       ! the points don't lie on a circle (since they are all on one line)
       ! --> return negative radius
       r = -1.0d0
       return
    end if

    f = 0.5d0 / f

    xyz1 = p1(1)*p1(1) + p1(2)*p1(2) + p1(3)*p1(3)
    xyz2 = p2(1)*p2(1) + p2(2)*p2(2) + p2(3)*p2(3)
    xyz3 = p3(1)*p3(1) + p3(2)*p3(2) + p3(3)*p3(3)
    xyz4 = p4(1)*p4(1) + p4(2)*p4(2) + p4(3)*p4(3)

    d21  = xyz2 - xyz1
    d31  = xyz3 - xyz1
    d41  = xyz4 - xyz1
    d32  = xyz3 - xyz2
    d42  = xyz4 - xyz2
    d43  = xyz4 - xyz3

    x0   =  f*( p1(3)*(p2(2)*d43 - p3(2)*d42 + p4(2)*d32) &
              - p2(3)*(p1(2)*d43 - p3(2)*d41 + p4(2)*d31) &
              + p3(3)*(p1(2)*d42 - p2(2)*d41 + p4(2)*d21) &
              - p4(3)*(p1(2)*d32 - p2(2)*d31 + p3(2)*d21) )

    y0   = -f*( p1(3)*(p2(1)*d43 - p3(1)*d42 + p4(1)*d32) &
              - p2(3)*(p1(1)*d43 - p3(1)*d41 + p4(1)*d31) &
              + p3(3)*(p1(1)*d42 - p2(1)*d41 + p4(1)*d21) &
              - p4(3)*(p1(1)*d32 - p2(1)*d31 + p3(1)*d21) )

    z0   =  f*( p1(2)*(p2(1)*d43 - p3(1)*d42 + p4(1)*d32) &
              - p2(2)*(p1(1)*d43 - p3(1)*d41 + p4(1)*d31) &
              + p3(2)*(p1(1)*d42 - p2(1)*d41 + p4(1)*d21) &
              - p4(2)*(p1(1)*d32 - p2(1)*d31 + p3(1)*d21) )

    dx = p1(1) - x0
    dy = p1(2) - y0
    dz = p1(3) - z0
    r = sqrt(dx*dx + dy*dy + dz*dz)

  end subroutine circumsphere

end module voronoi



