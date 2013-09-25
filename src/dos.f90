program dos

  !--------------------------------------------------------------------!
  ! dos.x --- calculate the density of states                          !
  !--------------------------------------------------------------------!
  ! 2011-06-01 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use constants, only: DP,                    &
                       Ha2eV, eV2Ha

  use energy,    only: nrg_occupancy, nrg_gaussian

  use geometry,  only: geo_init,              &
                       geo_final,             &
                       doDOSPlot,             &
                       DOSPoints,             &
                       DOSEmin,               &
                       DOSEmax,               &
                       DOSWidth

  use tbio,      only: tbio_init,             &
                       tbio_final,            &
                       tbio_args,             &
                       tbio_filecheck

  implicit none

  integer, parameter :: u_eig = 30
  integer, parameter :: u_out = 31

  !------------------------- k-point sampling -------------------------!
  ! nkpts          : number of k-points                                !
  ! wkpt(ikpt)     : weight of the ikpt-th k-point                     !
  ! kpt(1:3)       : coordinates (reciprocal units) of current k-point !
  ! gammaPoint     : .true., if gamma point calculation                !
  !                                                                    !
  !--------------------- eigenvalues and -vectors ---------------------!
  ! nev            : number of eigenvalues                             !
  ! eval(iev,ikpt) : iev-th eigenvalue at the ikpt-th k-point          !
  ! ev             : the current eigenvalue                            !
  !                                                                    !
  !----------------------------- DOS plot -----------------------------!
  ! Emin, Emax     : initial and final energy for the DOS plot         !
  ! npts           : number of data points in the DOS plot             !
  ! n, E           : denisty of states n(E)                            !
  !                                                                    !
  !----------------------------- options ------------------------------!
  ! inFile         : name of the input file with eigenvalues           !
  ! outFile        : name of the output file                           !
  ! eigFile        : name of the file with eigenvalues and -vectors    !
  ! fileTypee      : input file format                                 !
  ! Ezero          : zero energy value                                 !
  !--------------------------------------------------------------------!

  integer                                       :: ikpt, nkpts
  double precision, dimension(:),   allocatable :: wkpt
  double precision, dimension(3)                :: kpt
  logical                                       :: gammaPoint

  integer                                       :: iev, nev
  double precision, dimension(:,:), allocatable :: eval
  double precision                              :: ev, ev_min, ev_max

  double precision                              :: Emin, Emax
  integer                                       :: ipt, npts
  double precision                              :: n, dn
  double precision                              :: E, dE

  character(len=100)                            :: inFile, outFile
  character(len=100)                            :: eigFile, fileType
  double precision                              :: Ezero

  !----------------------------- options ------------------------------!

  call tbio_init(name  = 'dos.x',                                   &
                 descr = 'calculate the density of states',         & 
                 args  = (/ 'in    ', 'out   ', 'eig   ', 'ezero ', &
                            'format' /),                            &
                 extra = "" )

  call tbio_args('in',     inFile)
  call tbio_args('out',    outFile,  default="output.dos")
  call tbio_args('eig',    eigFile)
  call tbio_args('format', fileType)
  call tbio_args('ezero',  Ezero)

  call tbio_filecheck(eigFile)

  call tbio_final()


  !-------------------------- initialization --------------------------!

  call geo_init(inFile,fileType)

  if (.not. doDOSPlot) then
     write(0,*) "No DOS plot options in input file."
     call geo_final()
     stop
  end if


  !------------------------- read eigenvalues -------------------------!

  open(u_eig, file=trim(eigFile), status='old', form='unformatted', &
             action='read')

  read(u_eig) gammaPoint
  read(u_eig) nkpts, nev
  allocate(eval(nev,nkpts), wkpt(nkpts))

  do ikpt = 1, nkpts
     read(u_eig) kpt(1:3), wkpt(ikpt)
     read(u_eig) eval(1:nev,ikpt)
     read(u_eig)
  end do

  close(u_eig)


  !-------------------------- calculate DOS ---------------------------!

  ev_min = minval(eval(1:nev,1:nkpts))
  ev_max = maxval(eval(1:nev,1:nkpts))
  if (DOSEmin < DOSEmax) then
     Emin = max(DOSEmin + Ezero*eV2Ha, ev_min)
     Emax = min(DOSEmax + Ezero*eV2Ha, ev_max)
  else
     Emin = ev_min
     Emax = ev_max
  end if
  Emin = Emin - 2.0d0*DOSWidth
  Emax = Emax + 2.0d0*DOSWidth

  npts = DOSPoints - 1
  dE   = (Emax - Emin)/dble(npts + 1)

  open(u_out, file=trim(outFile), status='replace', action='write')

  write(u_out,'(1x,2(ES15.6,2x))') Emin*Ha2eV - Ezero, 0.0d0

  E  = Emin + dE
  do ipt = 1, npts

     n = 0.0d0
     do ikpt = 1, nkpts
        dn = 0.0d0
        do iev  = 1, nev
           ev = eval(iev, ikpt)
           dn = dn + nrg_gaussian(E, ev, DOSWidth)
        end do
        dn = wkpt(ikpt)*dn
        n = n + dn
     end do

     write(u_out,'(1x,2(ES15.6,2x))') E*Ha2eV - Ezero, 2.0d0*n*eV2Ha

     E  = E + dE

  end do

  close(u_out)

  !----------------------------- finalize -----------------------------!

  deallocate(eval, wkpt)
  call geo_final()

end program dos
