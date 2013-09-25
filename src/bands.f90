program bands

  !--------------------------------------------------------------------!
  ! bands -- calculate band structures                                 !
  !                                                                    !
  ! Usage: ./bands.x [options]                                         !
  !                                                                    !
  ! Options:                                                           !
  !                                                                    !
  ! -h  --help     show usage information                              !
  ! -d  --dir      specify directory with tight-binding parameters     !
  !                (default: ./tbparam)                                !
  ! -i  --input    specify input file name (default: INP)              !
  ! -o  --output   specify outbut file name (default: bandstruc.dat)   !
  ! -f  --format   specify format of the input file (default: mbpp)    !
  ! -p  --param    specify TB parameter set (default: proj01)          !
  ! -z  --zero     set energy zero to this value [eV]                  !
  ! --untangle     try to disentangle bands (default)                  !
  ! --untangleoff  turn of disentangling of bands                      !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! 2010-10-31 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use constants, only: CURDIR, DIRSEP,         &
                       DP,                     &
                       Ha2eV

  use geometry,  only: geo_init,               &
                       geo_final,              &
                       atomType,               &
                       atomTypeName,           &
                       bandPlotBands,          &
                       bandPlotPoints,         &
                       bandPlotPoint,          &
                       bandPlotStepSize,       &
                       cooLattice,             &
                       doBandPlot,             &
                       latticeVec,             &
                       recLattVec,             &
                       recLattVecMet,          &
                       nAtoms,                 &
                       nAtomsOfType,           &
                       nTypes

  use pbc,       only: pbc_d_max_origin,       &
                       pbc_number_of_tvecs,    &
                       pbc_compute_tvecs

  use tbio,      only: tbio_init,              &
                       tbio_final,             &
                       tbio_args
                 
  use tbmatrix,  only: mat_init,               &
                       mat_final,              &
                       mat_aoNum,              &
                       mat_r_max,              &
                       mat_get_nOrbitals,      &
                       mat_nOrbsMax,           &
                       mat_setup_z,            &
                       mat_setup_d,            &
                       mat_setup_0,            &
                       mat_setup_onsite,       &
                       mat_eigen_z,            &
                       mat_print_H,            &
                       mat_print_S,            &
                       mat_H_0, mat_S_0,       &
                       mat_H_d, mat_S_d,       &
                       mat_H_z, mat_S_z

  use tbparam,   only: tbp_init,               &
                       tbp_final,              &
                       tbp_H, tbp_S,           &
                       tbp_O, tbp_pot,         &
                       nl, l_of_il

  use voronoi,   only: voro_wsc,               &
                       convex_hull

  implicit none

  !--------------------------------------------------------------------!
  !                             constants                              !
  !                                                                    !
  ! NLMAX     : max. number of different angular momentum channels per !
  !             atom type                                              !
  !--------------------------------------------------------------------!

  integer,   parameter :: NLMAX  = 3
  integer,   parameter :: u_band = 20
  integer,   parameter :: u_gpl  = 30

  !--------------------------------------------------------------------!
  ! r_max          : max. radial cut-off for all TB matrix elements    !
  ! d_max          : max. cart. distance of any atom from the origin   !
  ! nTvecs         : number of real space translation vectors T within !
  !                  radial cut-off r_max                              !
  ! Tvec(i,ivec)   : i-th component of the ivec-th T vector            !
  ! TvecLen(ivec)  : norm of the ivec-th T vector                      !
  !                                                                    !
  !----------------------- bandstructure points -----------------------!
  ! kpt(1:3)       : k-point (rec. coordinates) for eigenval. calc.    !
  ! kvec(1:3)      : vector in k-space to the next k-point             !
  ! knorm          : ||kvec||                                          !
  ! dk             : fraction of knorm associated with one step        !
  ! kout           : measure for distance (x value) in rec. coo.       !
  ! runit          : knorm for the first line in the bandplot          !
  !                                                                    !
  !---------------------- eigenvalues / -vectors ----------------------!
  ! nev            : number of eigenvalues                             !
  ! eval(i)        : i-th eigenvalue                                   !
  ! evec_z(1:nev,i): i-th eigenvector                                  !
  ! evprev(1:2,i)  : i-th eigenvalues from the two previous k-points   !
  ! evidx(i)       : eigenvalue index ---> use: eval(evidx(i))         !
  ! ev_min, ev_max : min/max eigenvalue (over all k-points)            !
  !                                                                    !
  !----------------------------- options ------------------------------!
  ! infile         : name of the input file                            !
  ! filetype       : type/format of the input file                     !
  ! paramset       : name of the TB parameters set                     !
  ! paramdir       : directory with the TB parameters files            !
  ! gnuplot        : .true., if gnuplot script shall be created        !
  ! plotbz         : .true., if Brillouin zone shall be plotted        !
  ! untangle       : .true., if bands shall be desentangled            !
  ! ezero          : zero energy value                                 !
  !--------------------------------------------------------------------!
  
  double precision                              :: r_max, d_max
  integer                                       :: nTvecs
  integer,          dimension(:,:), allocatable :: Tvec
  double precision, dimension(:),   allocatable :: TvecLen

  double precision, dimension(3)                :: kvec, kpt
  double precision                              :: dk, knorm, kout, runit
  integer                                       :: ikpt1, ikpt2
  integer                                       :: istep, nsteps, i, k, l

  integer                                       :: nev
  double precision, dimension(:),   allocatable :: eval
  complex(kind=DP), dimension(:,:), allocatable :: evec_z
  double precision, dimension(:,:), allocatable :: evprev
  integer,          dimension(:),   allocatable :: evidx
  double precision                              :: ev_min, ev_max

  character(len=100) :: infile, filetype, paramset, paramdir
  character(len=100) :: outfile, gplfile
  logical            :: gnuplot, plotbz, nountangle
  double precision   :: ezero
  

  !----------------------------- options ------------------------------!

  call tbio_init(name  = 'bands.x',                                   &
                 descr = 'band structure calculations',               & 
                 args  = (/ 'in        ', 'out       ', 'format    ', &
                            'param     ', 'dir       ', 'gnuplot   ', &
                            'plotbz    ', 'nountangle', 'ezero     '  &
                         /),                                          &
                 extra = '')

  call tbio_args('in',         infile)
  call tbio_args('out',        outfile)
  call tbio_args('format',     filetype)
  call tbio_args('param',      paramset)
  call tbio_args('dir',        paramdir)
  call tbio_args('gnuplot',    gnuplot)
  call tbio_args('plotbz',     plotbz)
  call tbio_args('nountangle', nountangle)
  call tbio_args('ezero',      ezero)

  call tbio_final()

  gplfile = trim(adjustl(outfile)) // '.gpl'

  !-------------------------- initialization --------------------------!

  call geo_init(infile, filetype)
  call tbp_init(NLMAX, nTypes, atomTypeName, paramset, paramdir)

  call mat_init(nTypes, nAtoms, atomType, nl, NLMAX, l_of_il, tbp_S, 'complex')

  r_max = mat_r_max
  d_max = pbc_d_max_origin(nAtoms, cooLattice, latticeVec)
  nev = mat_get_nOrbitals()
  allocate(eval(nev), evec_z(nev,nev), evprev(2,nev), evidx(nev))
  do i = 1, nev
     evidx(i) = i
  end do

  ! translation vectors within the interaction radius:
  nTvecs = pbc_number_of_tvecs(r_max, d_max, latticeVec)
  allocate(Tvec(3,nTvecs), TvecLen(nTvecs))
  call pbc_compute_tvecs(r_max, d_max, latticeVec, nTvecs, Tvec, Tveclen)


  !---------------------- real space Hamiltonian ----------------------!

  ! set-up real-space Hamilton matrix:
  call mat_setup_0(nTypes, nAtoms, atomType, mat_aoNum, cooLattice,   &
                   latticeVec, nl, NLMAX, l_of_il, mat_nOrbsMax, nev, &
                   mat_H_0, mat_S_0)

  ! add distance dependent on-site contributions:
  call mat_setup_onsite(nTypes, nAtoms, atomType, mat_aoNum,   &
                        cooLattice, latticeVec, nTvecs, Tvec,  &
                        nl, NLMAX, l_of_il, mat_nOrbsMax, nev, &
                        tbp_O, mat_H_0, mat_S_0)


  !-------------------- band structure calculation --------------------!

  runit = 0.0d0
  band_structure_plot : if (doBandPlot) then

     open(u_band, file=trim(outfile), status='replace', action='write')

     ! BZ plot:
     if (plotbz) call bz_plot(recLattVec, bandPlotPoint, bandPlotPoints)

     if (gnuplot) then
        open(u_gpl,  file=trim(gplfile), status='replace', action='write')
        call gnuplot_header(u_gpl, bandPlotPoints, trim(adjustl(outfile)))
     end if

     ikpt1 = 1
     kout  = 0.0d0
     kpt   = bandPlotPoint(1:3,ikpt1) 

     ! starting point of the band structure:
     call mat_setup_z(nTypes, nAtoms, atomType, mat_aoNum,          &
                      cooLattice, latticeVec, nTvecs, Tvec, nl,     &
                      NLMAX, l_of_il, kpt(1:3), mat_nOrbsMax, nev,  &
                      mat_H_0, mat_S_0, mat_H_z, mat_S_z)

     call mat_eigen_z(nev, eval, evec_z)
     evprev(1,1:nev) = eval(1:nev)
     evprev(2,1:nev) = eval(1:nev)
     ev_min = eval(1)
     ev_max = eval(nev)
     call print_eigenvalues(u_band, nev, eval, bandPlotbands, kout, &
                            evidx, ezero)

     ! now loop over the other points:
     band_structure_points : do ikpt2 = 2, bandPlotPoints
   
        kvec(1:3) = bandPlotPoint(1:3,ikpt2) - bandPlotPoint(1:3,ikpt1)

        ! norm from reciprocal coordinates
        knorm = 0.0d0
        do k = 1, 3
        do l = 1, 3
           knorm = knorm + kvec(l)*recLattVecMet(l,k)*kvec(k)
        end do
        end do
        knorm = sqrt(knorm)

        if (runit == 0.0d0) then
           runit = knorm
        end if

        if (gnuplot) call gnuplot_addpoint(u_gpl, ikpt1, kout/runit)

        nsteps    = nint(knorm/bandPlotStepSize)
        kvec(1:3) = kvec(1:3)/dble(nsteps)
        dk        = knorm/dble(nsteps)

        steps : do istep = 1, nsteps

           kout     = kout + dk
           kpt(1:3) = bandPlotPoint(1:3, ikpt1) &
                    + dble(istep)*kvec(1:3)

           call mat_setup_z(nTypes, nAtoms, atomType, mat_aoNum,          &
                            cooLattice, latticeVec, nTvecs, Tvec, nl,     &
                            NLMAX, l_of_il, kpt(1:3), mat_nOrbsMax, nev,  &
                            mat_H_0, mat_S_0, mat_H_z, mat_S_z)

           call mat_eigen_z(nev, eval, evec_z)
           if (.not. nountangle) then
              call update_eval_index(nev, eval, evprev, evidx)
           end if
           ev_min = min(ev_min, eval(1))
           ev_max = max(ev_max, eval(nev))
           call print_eigenvalues(u_band, nev, eval, bandPlotBands, &
                                  kout/runit, evidx, ezero)
   
        end do steps

        ikpt1 = ikpt2

     end do band_structure_points

     if (gnuplot) then
        call gnuplot_addpoint(u_gpl, ikpt1, kout/runit)
        ev_min = ev_min*Ha2eV - ezero
        ev_max = ev_max*Ha2eV - ezero
        call gnuplot_footer(u_gpl, kout/runit, ev_min, ev_max, &
             min(nev,bandPlotBands), bandPlotPoints, trim(adjustl(outfile)))
     end if

  end if band_structure_plot


  !--------------------------- finalization ---------------------------!

  call mat_final()
  call tbp_final()
  call geo_final()

  if (gnuplot) close(u_gpl)
  close(u_band)
  deallocate(Tvec, TvecLen, eval, evec_z, evprev, evidx)

contains

  !--------------------------------------------------------------------!
  !            print eigenvalues to the band structure file            !
  !--------------------------------------------------------------------!

  subroutine print_eigenvalues(u_out, nev, eval, nbands, kout, evidx, ezero)

    implicit none

    integer,                          intent(in) :: u_out
    integer,                          intent(in) :: nev
    double precision, dimension(nev), intent(in) :: eval
    integer,                          intent(in) :: nbands
    double precision,                 intent(in) :: kout
    integer,          dimension(nev), intent(in) :: evidx
    double precision,                 intent(in) :: ezero

    integer                     :: ib, nb

    nb = min(nev, nbands)

    write(u_out, '(1x,F12.8,1x)', advance='no') kout
    do ib = 1, nb
       write(u_out, '(1x,F12.8,1x)', advance='no') &
            (Ha2eV*eval(evidx(ib)) - ezero)
    end do

    write(u_out,*)
    
  end subroutine print_eigenvalues
  
  !--------------------------------------------------------------------!
  !                     try to keep bands together                     !
  !--------------------------------------------------------------------!

  subroutine update_eval_index(nev, eval, evprev, evidx)

    implicit none

    integer,                            intent(in)    :: nev
    double precision, dimension(nev),   intent(in)    :: eval
    double precision, dimension(2,nev), intent(inout) :: evprev
    integer,          dimension(nev),   intent(inout) :: evidx

    double precision,      parameter :: dmin = 1.0d-3

    double precision, dimension(nev) :: evidx_new
    logical,          dimension(nev) :: set
    integer                          :: iev1, iev2, ev
    double precision                 :: evex, diff1, diff2

    set(1:nev) = .false.

    do iev1 = 1, nev
       evex  = evprev(1, iev1) + (evprev(1, iev1) - evprev(2, iev1))
       ev    = evidx(iev1)
       if (set(ev)) then
          diff1 = 100.0d0
       else
          diff1 = abs(evex - eval(ev))
       end if
       do iev2 = 1, nev
          if (ev == iev2) cycle
          if (set(iev2))  cycle
          diff2 = abs(evex - eval(iev2))
          if (diff1 - diff2 > dmin) then
             ev = iev2
             diff1 = diff2
          end if
       end do
       set(ev) = .true.
       evidx_new(iev1) = ev
    end do

    evidx(1:nev) = evidx_new(1:nev)
    evprev(2,1:nev) = evprev(1,1:nev)
    do iev1 = 1, nev
       evprev(1,iev1) = eval(evidx(iev1))
    end do

  end subroutine update_eval_index

  !--------------------------------------------------------------------!
  !                        write gnuplot script                        !
  !--------------------------------------------------------------------!

  subroutine gnuplot_header(u_gpl, npoints, bandsfile)

    implicit none

    integer,          intent(in) :: u_gpl
    integer,          intent(in) :: npoints
    character(len=*), intent(in) :: bandsfile
    
    integer            :: ipt
    character(len=100) :: str

    str = 'set style line 1 lc rgb "black" lt 1 lw 2.5'
    write(u_gpl, '(A)') trim(str)
    str = 'set style line 2 lc rgb "black" lt 1 lw 1.5'
    write(u_gpl, '(A)') trim(str)
    write(u_gpl, *)
    write(u_gpl, '("set border ls 2")')
    write(u_gpl, *)
    write(u_gpl, '("set terminal postscript eps color enhanced rounded \")')
    write(u_gpl, '("    lw 2.0 dl 3.0 font 30 size 8.0, 6.0")')
    write(u_gpl, *)

    write(u_gpl, '("# edit k-point labels here:")')
    do ipt = 1, npoints
       write(str, *) ipt
       str = 'P' // trim(adjustl(str))
       str = trim(adjustl(str)) // ' = "' // trim(adjustl(str)) // '"'
       write(u_gpl, '(A)') trim(str)
    end do
    write(u_gpl, *)

    write(u_gpl, '("set title ",A)')   '"band structure plot"'
    write(u_gpl, '("set ylabel ", A)') '"eigenvalue [eV]"'
    write(u_gpl, '("unset xlabel")')
    write(u_gpl, '("unset key")')
    str = '"' // trim(adjustl(bandsfile)) // '.eps' // '"'
    write(u_gpl, '("set output ",A)')  trim(str)
    write(u_gpl, *)

  end subroutine gnuplot_header

  !--------------------------------------------------------------------!

  subroutine gnuplot_addpoint(u_gpl, ipt, x)

    implicit none

    integer,          intent(in) :: u_gpl, ipt
    double precision, intent(in) :: x

    character(len=50) :: str1, str2

    write(str1, *) ipt
    write(str2, *) x
    str1 = 'x' // trim(adjustl(str1)) // ' = ' // trim(adjustl(str2))
    write(u_gpl, '(A)') trim(str1)

  end subroutine gnuplot_addpoint

  !--------------------------------------------------------------------!

  subroutine gnuplot_footer(u_gpl, x_max, y_min, y_max, nlines, npoints, &
                            bandsfile)

    implicit none

    integer,          intent(in) :: u_gpl
    double precision, intent(in) :: x_max, y_min, y_max
    integer,          intent(in) :: nlines, npoints
    character(len=*), intent(in) :: bandsfile

    integer            :: iline, ipt
    character(len=100) :: str, pi, xi, frmt

    write(u_gpl, *)
    write(u_gpl, '("x0 = 0.0")')
    write(str, *) x_max
    write(u_gpl, '("xf = ", A)') trim(adjustl(str))
    write(str, *) y_min
    write(u_gpl, '("y0 = ", A)') trim(adjustl(str))
    write(str, *) y_max
    write(u_gpl, '("yf = ", A)') trim(adjustl(str))
    write(u_gpl, *)    
    write(u_gpl, '("set xrange [x0:xf]")')
    write(u_gpl, '("set yrange [y0:yf]")')
    write(u_gpl, *)    

    write(u_gpl, '("set xtics (P1 x1)")')
    do ipt = 2, npoints
       write(str, *) ipt
       xi = 'x' // trim(adjustl(str))
       pi = 'P' // trim(adjustl(str))
       write(u_gpl, '("set xtics add (",A," ",A,")")') trim(pi), trim(xi)
    end do
    write(u_gpl, '("set grid xtics front ls 2")')
    write(u_gpl, *)

    frmt = "'" // '"' // trim(adjustl(bandsfile)) // '"' // " u 1:',A,' w l ls 1'"
    write(u_gpl, '("plot ",'// frmt //',", \"'//')') '2'
    do iline = 2, nlines-1
       write(str, *) iline+1
       write(u_gpl, '("     ",'// frmt //',", \"'//')') trim(adjustl(str))
    end do
    write(str, *) nlines+1
    write(u_gpl, '("     ",'// frmt // ')') trim(adjustl(str))
    write(u_gpl, *)

    write(u_gpl, '("! ps2pdf -dEPSCrop ", A)') trim(adjustl(bandsfile)) // '.eps'
    write(u_gpl, *)
    write(u_gpl, '("exit 0")')

  end subroutine gnuplot_footer

  !--------------------------------------------------------------------!
  !                        Brillouin zone plot                         !
  !--------------------------------------------------------------------!

  subroutine bz_plot(recLattVec, kpt, nkpts)

    implicit none

    double precision, dimension(3,3),     intent(in) :: recLattVec
    integer,                              intent(in) :: nkpts
    double precision, dimension(3,nkpts), intent(in) :: kpt

    integer,          parameter :: u_bz     = 50
    integer,          parameter :: u_plt    = 51
    character(len=*), parameter :: bz_file  = 'bz.dat'
    character(len=*), parameter :: plt_file = 'bz.gpl'

    integer                                       :: nBZpoints
    double precision, dimension(:,:), allocatable :: BZpoint
    integer                                       :: nBZedges
    integer,          dimension(:,:), allocatable :: BZedge

    double precision               :: rmax
    double precision, dimension(3) :: K1, K2
    character(len=200)             :: frmt
    integer                        :: iedge, ikpt

    !-------------------------------------------!
    ! calculate BZ edges and store them to file !
    !-------------------------------------------!

    nBZpoints = 0
    call voro_wsc(recLattVec, BZpoint, nBZpoints)
    if (nBZpoints < 0) then
       nBZpoints = -nBZpoints
       allocate(BZpoint(3,nBZpoints))
       call voro_wsc(recLattVec, BZpoint, nBZpoints)
    end if

    nBZedges = 0
    call convex_hull(BZpoint, nBZpoints, BZedge, nBZedges)
    if (nBZedges < 0) then
       nBZedges = -nBZedges
       allocate(BZedge(2,nBZedges))
       call convex_hull(BZpoint, nBZpoints, BZedge, nBZedges)
    end if
    
    open(u_bz, file=trim(bz_file), status='replace', action='write')
    do iedge = 1, nBZedges
       write(u_bz,'(1x,3(F12.8,2x))') BZpoint(1:3,BZedge(1,iedge))
       write(u_bz,'(1x,3(F12.8,2x))') BZpoint(1:3,BZedge(2,iedge))
       write(u_bz,*)
       write(u_bz,*)
    end do
    close(u_bz)

    deallocate(BZpoint, BZedge)

    !-------------------------!
    ! generate gnuplot script !
    !-------------------------!

    rmax = maxval(recLattVec(:,:))

    open(u_plt, file=trim(plt_file), action='write', status='replace')

    write(u_plt, '(A)') 'set style line 1 lt 1 lw 2.0 lc rgb "black" pt 7 ps 2.5'
    write(u_plt, '(A)') 'set style line 2 lt 1 lw 2.0 lc rgb "red"   pt 7 ps 2.5'
    write(u_plt, '(A)') 'set style line 3 lt 1 lw 1.0 lc rgb "black" pt 7 ps 2.5'
    write(u_plt, *)
    frmt = '("set arrow ",I2," from ",F10.6,", "F10.6,", "F10.6," to '
    frmt = trim(frmt) // '",F10.6,", "F10.6,", "F10.6," ls ",I2)'
    write(u_plt, frmt) 1, 0.0,0.0,0.0, recLattVec(1:3,1), 1
    write(u_plt, frmt) 2, 0.0,0.0,0.0, recLattVec(1:3,2), 1
    write(u_plt, frmt) 3, 0.0,0.0,0.0, recLattVec(1:3,3), 1
    write(u_plt, *)
    do ikpt = 2, nkpts
       K1(1:3) = kpt(1,ikpt-1)*recLattVec(1:3,1) &
               + kpt(2,ikpt-1)*recLattVec(1:3,2) &
               + kpt(3,ikpt-1)*recLattVec(1:3,3)
       K2(1:3) = kpt(1,ikpt)*recLattVec(1:3,1) &
               + kpt(2,ikpt)*recLattVec(1:3,2) &
               + kpt(3,ikpt)*recLattVec(1:3,3)
       write(u_plt, frmt) ikpt+2, K1, K2, 2
    end do
    write(u_plt, *)
    frmt = '("set label ",I2," ",A," at ",F10.6,", "F10.6,", "F10.6,'
    frmt = trim(frmt) // '" nopoint offset 0.5")'
    write(u_plt, frmt) 1, '"b1"', 0.8d0*recLattVec(1:3,1)
    write(u_plt, frmt) 2, '"b2"', 0.8d0*recLattVec(1:3,2)
    write(u_plt, frmt) 3, '"b3"', 0.8d0*recLattVec(1:3,3)
    write(u_plt, *)
    write(u_plt, '(A)') 'set terminal x11 enhanced'
    write(u_plt, '(A)') 'set size square'
    write(u_plt, '(A)') 'set ticslevel 0'
    write(u_plt, '(A)') 'unset xtics'
    write(u_plt, '(A)') 'unset ytics'
    write(u_plt, '(A)') 'unset ztics'
    write(u_plt, '(A)') 'unset border'
    write(u_plt, '(A)') 'unset key'
    write(frmt,*) rmax
    write(u_plt, '("set xrange [-",A,":",A,"]")') trim(adjustl(frmt)), trim(adjustl(frmt))
    write(u_plt, '("set yrange [-",A,":",A,"]")') trim(adjustl(frmt)), trim(adjustl(frmt))
    write(u_plt, '("set zrange [-",A,":",A,"]")') trim(adjustl(frmt)), trim(adjustl(frmt))
    write(u_plt, *)
    write(u_plt, '(A)') 'splot "' // trim(bz_file) // '" u 1:2:3 w lp ls 3'
    write(u_plt, *)
    write(u_plt, '(A)') 'pause -1'
    write(u_plt, *)

    close(u_plt)
    

  end subroutine bz_plot

end program bands
