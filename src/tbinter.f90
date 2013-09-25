module tbinter

  !--------------------------------------------------------------------!
  ! generalized Tight-Binding interface                                !
  !--------------------------------------------------------------------!
  ! This module offers interfaces to the tight-binding energy routines.!
  !--------------------------------------------------------------------!
  ! 2011-03-18 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use constants, only: DP, Ha2eV

  use energy,    only: nrg_atomic,          &
                       nrg_fermi_level,     &
                       nrg_occupancy

  use geometry,  only: geo_init,              &
                       geo_final,             &
                       geo_print_status,      &
                       geo_print_lattice,     &
                       geo_print_coordinates, &
                       atomType,              &
                       atomTypeName,          &
                       cooLattice,            &
                       kptMultiple,           &
                       kptShift,              &
                       latticeVec,            &
                       recLattVec,            &
                       metallic,              &
                       nAtoms,                &
                       nAtomsOfType,          &
                       nTypes,                &
                       doForces,              &
                       useSymmetry,           &
                       smearWidth

  use io,        only: io_adjustl,            &
                       io_unit

  use pbc,       only: pbc_d_max_origin,      &
                       pbc_MP_kpoint_set,     &
                       pbc_number_of_tvecs,   &
                       pbc_compute_tvecs

  use potential, only: pot_init,              &
                       pot_final,             &
                       pot_pairpot,           &
                       pot_r_max

  use symmetry,  only: sym_ir_kpoint_set,     &
                       sym_get_symmetry,      &
                       sym_print_kpoints
                 
  use tbmatrix,  only: mat_init,              &
                       mat_final,             &
                       mat_aoNum,             &
                       mat_elem_full_d,       &
                       mat_elem_full_d_alt,   &
                       mat_elem_full_z,       &
                       mat_delem_full_d,      &
                       mat_delem_full_z,      &
                       mat_eigen_d,           &
                       mat_eigen_z,           &
                       mat_get_nOrbitals,     &
                       mat_nOrbs,             &
                       mat_nOrbsMax,          &
                       mat_r_max,             &
                       mat_setup_onsite,      &
                       mat_setup_0,           &
                       mat_setup_d,           &
                       mat_setup_z,           &
                       mat_setup_full_d_alt,  &
                       mat_H_0, mat_S_0,      &
                       mat_H_d,  mat_S_d,     &
                       mat_H_z,  mat_S_z,     &
                       mat_print_matrix

  use tbparam,   only: tbp_init,              &
                       tbp_final,             &
                       tbp_level,             &
                       tbp_nElec,             &
                       tbp_H, tbp_S,          &
                       tbp_O, tbp_pot,        &
                       nl, l_of_il

  use timing,    only: tng_init,              &
                       tng_final,             &
                       tng_timing,            &
                       tng_timing2

  !--------------------------------------------------------------------!

  implicit none

  public  :: tbi_init,                 &
             tbi_final,                &
             tbi_reset,                &
             tbi_energies,             &
             tbi_energy_analysis,      &
             tbi_forces,               &
             tbi_kpoint_sampling,      &
             tbi_print_eigenvalues,    &
             tbi_print_energies,       &
             tbi_print_forces,         &
             tbi_print_mulliken,       &
             tbi_print_occupation,     &
             tbi_print_system,         &
             tbi_setup_T0
             
  private :: tbi_band_sum,             &
             tbi_calc_forces_gamma,    &
             tbi_calc_forces_kpt,      &
             tbi_energy_partitioning,  &
             tbi_Epolar_Ecov_gamma,    &
             tbi_Epolar_Ecov_kpt,      &
             tbi_Eprom_Ecf,            &
             tbi_fermi_level,          &
             tbi_init_calc,            &
             tbi_save_density_matrix
             

  integer, parameter :: LENPARA = 100
  integer, parameter :: NLMAX   = 3

  integer            :: u_eig
  integer            :: u_rho
  integer            :: u_occ
  integer            :: u_hmat
  integer            :: u_smat
  integer, parameter :: u_dbg   = 99

  !------------------------------ files -------------------------------!
  ! inFile         : input file name                                   !
  ! outFile        : output file name                                  !
  ! eigFile        : binary output file for eigenvalues/-vectors       !
  ! rhoFile        : binary output file for the density matrix         !
  ! occupFile      : file with occupation numbers for every band on    !
  !                  separate lines                                    !
  !                                                                    !
  !----------------------------- options ------------------------------!
  ! fileType       : input file type                                   !
  ! readOccup      : if .true. --> read occupation numbers from file   !
  ! doTiming       : if .true. --> log timing information              !
  ! prtMatrices    : if .true. --> print H and S matrices              !
  !                                                                    !
  !------------------------- k-point sampling -------------------------!
  ! nkpts          : number of k-points                                !
  ! kpt(1:3,ikpt)  : ikpt-th k-point in set                            !
  ! wkpt(ikpt)     : weight of the ikpt-th k-point                     !
  ! gammaPoint     : .true. if the calculation is Gamma point only     !
  !                                                                    !
  !--------------------- eigenvalues and -vectors ---------------------!
  ! nev            : number of eigenvalues                             !
  ! nao            : number of solved bands (= nev)                    !
  ! nOccup         : number of occupied bands                          !
  ! fOccup(i,ikpt) : occupation number of band i at k-point ikpt       !
  ! eval(iev,ikpt) : iev-th eigenvalue at k-point ikpt                 !
  ! evec_?(i,iev)  : i-th component of the iev-th eigenvector for the  !
  !                  last computed k-point; ?=d -> dble, =z -> cmplx   !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! r_max          : upper radial cut-off for all interactions         !
  ! d_max          : max. cart. distance of any atom from the origin   !
  ! nTvecs         : number of real space translation vectors T within !
  !                  radial cut-off r_max                              !
  ! Tvec(i,ivec)   : i-th component of the ivec-th T vector            !
  ! TvecLen(ivec)  : norm of the ivec-th T vector                      !
  ! nElecs         : total number of electrons in the system           !
  !                                                                    !
  ! H_ii(i)        : diagonal element H(i,i) of real-space H matrix    !
  ! Qatoms(il,it)  : atomic occupation of level il of type it          !
  !                                                                    !
  ! Qmull(iao)     : Mulliken charge of basis function iao             !
  !                                                                    !
  !------------------------- energy analysis --------------------------!
  ! [ Bester and Faehnle, J. Phys.: Cond. Matter 13 (2001) 11541 ]     !
  ! flevel         : the Fermi level                                   !
  ! fermIter       : number of iterations to converge the Fermi level  !
  ! Esmear         : broadening correction                             !
  ! Eband          : band energy                                       !
  ! Eatoms         : energy of the free atoms                          !
  ! Epairpot       : pair potential energy                             !
  ! EtbBand        : total TB-Band energy                              !
  ! Eprom          : promotion energy                                  !
  ! Ecf            : crystal field energy                              !
  ! Epolar         : polarization energy                               !
  ! Ecov           : covalent bond energy                              !
  ! Echeck         : = E_cov + E_polar                                 !
  ! EtbBand        : TB-Bond cohesive energy                           !
  !                                                                    !
  !------------------------------ forces ------------------------------!
  ! forCart(i,iat) : i-th component of the cartesian force vector on   !
  !                  atom iat                                          !
  !--------------------------------------------------------------------!
  
  character(len=LENPARA),                        private :: inFile, outFile
  character(len=LENPARA),                        private :: paramSet, paramDir
  character(len=LENPARA),                        private :: eigFile, rhoFile
  character(len=LENPARA),                        private :: occupFile

  character(len=LENPARA),                        private :: fileType
  logical,                                       private :: readOccup
  logical,                                       private :: doTiming
  logical,                                       private :: prtMatrices

  integer,                                       private :: nkpts
  double precision, dimension(:,:), allocatable, private :: kpt
  double precision, dimension(:),   allocatable, private :: wkpt
  logical,                                       private :: gammaPoint

  integer,                                       private :: nev
  integer,                                       private :: nao
  integer,                                       private :: nOccup
  double precision, dimension(:,:), allocatable, private :: fOccup
  double precision, dimension(:,:), allocatable, private :: eval
  double precision, dimension(:,:), allocatable, private :: evec_d
  complex(kind=DP), dimension(:,:), allocatable, private :: evec_z

  double precision,                              private :: r_max, d_max
  integer,                                       private :: nTvecs
  integer,          dimension(:,:), allocatable, private :: Tvec
  double precision, dimension(:),   allocatable, private :: TvecLen

  integer,                                       private :: nElecs

  double precision, dimension(:),   allocatable, private :: H_ii
  double precision, dimension(:,:), allocatable, private :: Qatoms

  double precision, dimension(:),   allocatable, private :: Qmull

  double precision,                              private :: flevel
  integer,                                       private :: fermIter
  double precision,                              private :: Esmear
  double precision,                              private :: Eband
  double precision,                              private :: Eatoms
  double precision,                              private :: Epairpot
  double precision,                              private :: EtbBand

  double precision,                              private :: Eprom
  double precision,                              private :: Ecf
  double precision,                              private :: Epolar
  double precision,                              private :: Ecov
  double precision,                              private :: Echeck
  double precision,                              private :: EtbBond

  double precision, dimension(:,:), allocatable, target, private :: forCart

  !-------------------------- module status ---------------------------!
  ! memSize           : currently allocated memory size (words)        !
  !                                                                    !
  ! isInit            : .true. --> module has been initialized         !
  ! hasT0Matrix       : .true. --> H/S matrices for T=(0,0,0) are done !
  ! hasEigenvalues    : .true. --> GEVPs were solved                   !
  ! hasEnergies       : .true. --> energies and occupations calculated !
  ! hasDensityMatrix  : .true. --> the density-matrix has been  stored !
  ! hasMulliken       : .true. --> Mulliken charges are in memory      !
  ! hasEnergyAnalysis : .true. --> energy partitioning was performed   !
  ! hasForces         : .true. --> atomic forces have been calculated  !
  !--------------------------------------------------------------------!

  integer, private :: memSize           = 0

  logical, private :: isInit            = .false.
  logical, private :: hasEnergyAnalysis = .false.
  logical, private :: hasT0Matrix       = .false.
  logical, private :: hasEigenvalues    = .false.
  logical, private :: hasEnergies       = .false.
  logical, private :: hasDensityMatrix  = .false.
  logical, private :: hasForces         = .false.
  logical, private :: hasMulliken       = .false.



contains !=============================================================!

  
  subroutine tbi_init(in, ftype, pset, pdir, out, occup, time, forces, &
                      matrices)

    implicit none

    character(len=*),  intent(in) :: in, ftype
    character(len=*),  intent(in) :: pset, pdir
    character(len=*),  intent(in) :: out
    character(len=*),  intent(in) :: occup
    logical,           intent(in) :: time
    logical, optional, intent(in) :: forces
    logical, optional, intent(in) :: matrices

    if (isInit) then
       write(0,*) "Warning: multiple initialization (tbinter)."
       return
    else
       isInit = .true.
    end if

    inFile   = trim(in)
    fileType = trim(ftype)
    paramSet = trim(pset)
    paramDir = trim(pdir)
    outFile  = trim(out)
    eigFile  = trim(outfile) // '.eig'
    rhoFile  = trim(outfile) // '.rho'

    if (len_trim(occup) > 0) then
       readOccup = .true.
       occupFile = occup
    else
       readOccup = .false.
    end if

    if (time) then
       doTiming = .true.
       call tng_init()
       call tng_timing('Start of timing.')
    end if
    call geo_init(inFile, fileType)
    call tbp_init(NLMAX, nTypes, atomTypeName, paramSet, paramDir)
    call pot_init()

    ! the caller can override the input file settings:
    if (present(forces)) then
       if (forces) doForces = .true.
    end if

    ! shall the matrices be saved?
    prtMatrices = .false.
    if (present(matrices)) then
       if (matrices) then
          prtMatrices = .true.
          u_hmat = io_unit()
          open(u_hmat, file=trim(outfile)//'.hmat', status='replace', &
                       action='write')
          u_smat = io_unit()
          open(u_smat, file=trim(outfile)//'.smat', status='replace', &
                       action='write')
       end if
    end if

    ! allocate memory and set-up calculation:
    call tbi_init_calc()

    if (doTiming) call tng_timing('Initialization done.')

  end subroutine tbi_init

  !--------------------------------------------------------------------!

  subroutine tbi_final()

    implicit none

    if (.not. isInit) then
       write(0,*) "Warning: nothing to finalize (tbinter)"
       return
    end if

    ! free memory:
    deallocate(kpt, wkpt, eval, Tvec, TvecLen, Qatoms, H_ii, fOccup)
    
    memSize = memSize - nev*nkpts*2 - 3*nTvecs - 2*nTvecs &
                      - NLMAX*nTypes - 2*nev - nev*nkpts*2
    memSize = memSize - 8*kptMultiple(1)*kptMultiple(2)*kptMultiple(3)
    if (gammaPoint) then
       deallocate(evec_d)
       memSize = memSize - nev*nev*2
    else
       deallocate(evec_z)
       memSize = memSize - nev*nev*4
    end if
    if (doForces) then
       deallocate(forCart)
       memSize = memSize - 3*nAtoms*2
    end if
    if (hasMulliken) then
       deallocate(Qmull)
       memSize = memSize - nev*2
    end if

    ! close files:
    if (prtMatrices) then
       close(u_hmat)
       close(u_smat)
    end if

    ! finalize modules:
    call pot_final()
    call tbp_final()
    call geo_final()
    if (doTiming) call tng_final()

    ! reset module status:
    call tbi_reset()
    isInit            = .false.

    if (memSize /= 0) then
       write(0,*) "Warning: memory leak in `tbinit': ", memSize, &
                  " words are still allocated."
    end if

  end subroutine tbi_final

  !--------------------------------------------------------------------!
  ! reset the module state - e.g. before next iteration in relaxation  !
  !--------------------------------------------------------------------!

  subroutine tbi_reset()

    implicit none

    hasEnergyAnalysis = .false.
    hasT0Matrix       = .false.
    hasEigenvalues    = .false.
    hasEnergies       = .false.
    hasDensityMatrix  = .false.
    hasForces         = .false.
    hasMulliken       = .false.

  end subroutine tbi_reset




  !============================== PUBLIC ==============================!




  !--------------------------------------------------------------------!
  !         set-up Hamilton and overlap matrices for T=(0,0,0)         !
  !--------------------------------------------------------------------!

  subroutine tbi_setup_T0()
    
    implicit none

    integer :: iev

    if (.not. isInit) then
       write(0,*) "Error: Module `tbinter' has not been initialized."
       return
    end if

    ! set-up Hamilton matrix for T = (0,0,0):
    call mat_setup_0(nTypes, nAtoms, atomType, mat_aoNum, cooLattice,   &
                     latticeVec, nl, NLMAX, l_of_il, mat_nOrbsMax, nev, &
                     mat_H_0, mat_S_0)

    ! add on-site elements:
    call mat_setup_onsite(nTypes, nAtoms, atomType, mat_aoNum,   &
                          cooLattice, latticeVec, nTvecs, Tvec,  &
                          nl, NLMAX, l_of_il, mat_nOrbsMax, nev, &
                          tbp_O, mat_H_0, mat_S_0)

    do iev = 1, nev
       H_ii(iev) = mat_H_0(iev,iev)
    end do

    if (doTiming) call tng_timing('Matrices for T = (0,0,0) done.')

    hasT0Matrix = .true.

!DEBUG
!    open(88,file='smat0.dat')
!    call mat_print_matrix(mat_S_0, 88)
!    close(88)
!    open(88,file='hmat0.dat')
!    call mat_print_matrix(mat_H_0, 88)
!    close(88)
!END DEBUG

  end subroutine tbi_setup_T0

  !--------------------------------------------------------------------!
  !                          k-point sampling                          !
  !--------------------------------------------------------------------!

  subroutine tbi_kpoint_sampling()

    implicit none

    integer :: ikpt

    if (.not. isInit) then
       write(0,*) "Error: Module `tbinter' has not been initialized."
       return
    else if (.not. hasT0Matrix) then
       write(0,*) "Error: Matrices for T = (0,0,0) are missing."
       return
    end if

    ! open eigenvalues file:
    u_eig = io_unit()
    open(u_eig, file=trim(eigfile), status='replace', action='write', &
         form='unformatted')
    write(u_eig) gammaPoint
    write(u_eig) nkpts, nev
    if (gammaPoint) then
       call mat_setup_d(nTypes, nAtoms, atomType, mat_aoNum,        &
                        cooLattice, latticeVec, nTvecs, Tvec, nl,   &
                        NLMAX, l_of_il, mat_nOrbsMax, nev, mat_H_0, &
                        mat_S_0, mat_H_d, mat_S_d)
       if (doTiming) then
          call tng_timing2('Bloch matrix done (gamma-point)')
       end if
       if (prtMatrices) then
          call mat_print_matrix(mat_H_d, u_hmat)
          call mat_print_matrix(mat_S_d, u_smat)
       end if
       call mat_eigen_d(nev, eval(1:nev,1), evec_d)
       if (doTiming) then
          call tng_timing2('GEVP done (k-point ' & 
                          // trim(io_adjustl(ikpt)) // ')')
       end if
       ! save results to eigenvalues file:
       write(u_eig) kpt(1:3,1), wkpt(1)
       write(u_eig) eval(1:nev,1)
       write(u_eig) evec_d(1:nev,1:nev)
    else
       kpoint : do ikpt = 1, nkpts
          call mat_setup_z(nTypes, nAtoms, atomType, mat_aoNum,          &
                           cooLattice, latticeVec, nTvecs, Tvec, nl,     &
                           NLMAX, l_of_il, kpt(1:3,ikpt), mat_nOrbsMax,  &
                           nev, mat_H_0, mat_S_0, mat_H_z, mat_S_z)
          if (doTiming) then
             call tng_timing2('Bloch matrix done (k-point ' & 
                              // trim(io_adjustl(ikpt)) // ')')
          end if
          if (prtMatrices) then
             write(u_hmat,'(1x,I5,1x,3(1x,ES15.8))') ikpt, kpt(1:3,ikpt)
             write(u_hmat,*)
             call mat_print_matrix(mat_H_z, u_hmat)
             write(u_hmat,*)
             write(u_smat,'(1x,I5,1x,3(1x,ES15.8))') ikpt, kpt(1:3,ikpt)
             write(u_smat,*)
             call mat_print_matrix(mat_S_z, u_smat)
             write(u_smat,*)
          end if
          call mat_eigen_z(nev, eval(1:nev,ikpt), evec_z)
          if (doTiming) then
             call tng_timing2('GEVP done (k-point ' &
                             // trim(io_adjustl(ikpt))//')')
          end if
          ! save results to eigenvalues file:
          write(u_eig) kpt(1:3,ikpt), wkpt(ikpt)
          write(u_eig) eval(1:nev,ikpt)
          write(u_eig) evec_z(1:nev,1:nev)
       end do kpoint
    end if
    close(u_eig)
    
    if (doTiming) call tng_timing('k-point sampling done.')

    hasEigenvalues = .true.

  end subroutine tbi_kpoint_sampling

  !--------------------------------------------------------------------!
  !                          energy analysis                           !
  !--------------------------------------------------------------------!

  subroutine tbi_energies(E_fermi, E_smear, E_band, E_atoms, &
                          E_pairpot, E_coh)

    implicit none

    double precision, intent(out) :: E_fermi
    double precision, intent(out) :: E_smear
    double precision, intent(out) :: E_band
    double precision, intent(out) :: E_atoms
    double precision, intent(out) :: E_pairpot
    double precision, intent(out) :: E_coh

    integer :: itype

    if (.not. isInit) then
       write(0,*) "Error: Module `tbinter' has not been initialized."
       return
    else if (.not. hasEigenvalues) then
       write(0,*) "Error: The eigenvalues/-vectors are missing."
       return
    end if

    ! number of electrons in the system:
    nElecs = 0
    do itype = 1, nTypes
       nElecs = nElecs + nAtomsOfType(itype)*tbp_nElec(itype)
    end do

    ! Fermi level and band occupations:
    call tbi_fermi_level()

    ! energy and charges of the free atoms:
    call nrg_atomic(nTypes, nAtomsOfType, tbp_nElec, NLMAX, nl, &
         tbp_level, l_of_il, Eatoms, Qatoms)

    ! band energy:
    call tbi_band_sum()

    ! pair potential energy:
    if (tbp_pot) then
       if (doForces) then
          call pot_pairpot(nAtoms, cooLattice, atomType, latticeVec, &
                           nTvecs, Tvec, TvecLen, d_max,             &
                           energy = Epairpot, forces = forCart)
       else
          call pot_pairpot(nAtoms, cooLattice, atomType, latticeVec, &
                           nTvecs, Tvec, TvecLen, d_max,             &
                           energy = Epairpot)
       end if
       if (doTiming) call tng_timing('Pair potential done.')
    else
       Epairpot = 0.0d0
    end if

    E_fermi   = flevel
    E_band    = Eband
    E_atoms   = Eatoms
    E_pairpot = Epairpot
    E_smear   = Esmear
    E_coh     = E_band + E_smear + E_pairpot - E_atoms

    hasEnergies = .true.
    
  end subroutine tbi_energies

  !--------------------------------------------------------------------!
  !                        run energy analysis                         !
  !--------------------------------------------------------------------!

  subroutine tbi_energy_analysis(E_prom, E_cf, E_polar, E_cov, E_bond)

    implicit none

    double precision, intent(out) :: E_prom 
    double precision, intent(out) :: E_cf   
    double precision, intent(out) :: E_polar
    double precision, intent(out) :: E_cov  
    double precision, intent(out) :: E_bond 

    if (.not. isInit) then
       write(0,*) "Error: Module `tbinter' has not been initialized."
       return
    else if (.not. hasEnergies) then
       write(0,*) "Error: The energies have not yet been calculated."
       return
    end if

    ! calculate the density matrix rho:
    if (.not. hasDensityMatrix) then
       call tbi_save_density_matrix()
       if (doTiming) call tng_timing('Density matrix done.')   
    end if

    ! Energy partitioning:
    call tbi_energy_partitioning()
    if (doTiming) call tng_timing('Energy partitioning done.')

    E_prom  = Eprom
    E_cf    = Ecf
    E_polar = Epolar
    E_cov   = Ecov
    E_bond  = EtbBond

    hasEnergyAnalysis = .true.

  end subroutine tbi_energy_analysis

  !--------------------------------------------------------------------!
  !                               forces                               !
  !--------------------------------------------------------------------!

  subroutine tbi_forces(F)

    implicit none

    double precision, dimension(:,:), pointer, optional, intent(out) :: F

    if (.not. isInit) then
       write(0,*) "Error: Module `tbinter' has not been initialized."
       return
    else if (.not. hasEnergies) then
       write(0,*) "Error: The energies have not yet been calculated."
       return
    end if

    if (gammaPoint) then
       call tbi_calc_forces_gamma()
    else
       call tbi_calc_forces_kpt()
    end if
    if (doTiming) call tng_timing('Forces done.')

    if (present(F)) F => forCart

    hasForces = .true.

  end subroutine tbi_forces
  


  !============================== I / O ===============================!


  !--------------------------------------------------------------------!
  !              print out system and calculation details              !
  !--------------------------------------------------------------------!

  subroutine tbi_print_system()

    implicit none

    if (.not. isInit) return

    call geo_print_status()

    write(*,*) 'System summary'
    write(*,*) '=============='
    write(*,*)
    write(*,'(1x,"Number of T vectors    : ",A)') &
         trim(io_adjustl(nTvecs+1))
    write(*,'(1x,"Size of matrices       : ",A)') &
         trim(io_adjustl(nev))//'x'//trim(io_adjustl(nev))
    write(*,'(1x,"Gamma point algorithms : ",L1)') gammaPoint
    write(*,'(1x,"Overlap matrix         : ",L1)') tbp_S
    write(*,'(1x,"On-site parameters     : ",L1)') tbp_O
    write(*,'(1x,"Pair potential         : ",L1)') tbp_pot
    write(*,*)
    write(*,*)

    call geo_print_lattice()
    call geo_print_coordinates()

    if (useSymmetry) then
       call sym_get_symmetry(nAtoms, latticeVec, cooLattice, atomType)
       call sym_print_kpoints(nkpts, kpt, wkpt)
    end if

  end subroutine tbi_print_system

  !--------------------------------------------------------------------!
  !                         print out energies                         !
  !--------------------------------------------------------------------!

  subroutine tbi_print_energies()

    implicit none

    if (.not. hasEnergies) return

    write(*,*) 'Energies'
    write(*,*) '========'
    write(*,*)

    write(*,*) "There are " // trim(io_adjustl(nElecs)) // &
               " electrons in " // trim(io_adjustl(nOccup)) // " bands."
    write(*,*) "Fermi level converged after " &
              // trim(io_adjustl(fermIter)) // " iterations."
    write(*,*)

    write(*,'(1x,"Fermi level           : ",F20.8," eV")') flevel*Ha2eV
    write(*,*)
    write(*,'(1x,"Band energy           : ",F20.8," Ha")') Eband
    write(*,'(1x,"Broadening correction : ",F20.8," Ha")') Esmear
    write(*,'(1x,"Atomic energy         : ",F20.8," Ha")') Eatoms
    write(*,'(1x,"Pair potential energy : ",F20.8," Ha")') Epairpot
    write(*,'(1x,"------------------------------------------------")')
    EtbBand = Eband + Esmear + Epairpot - Eatoms
    write(*,'(1x,"TB-Band coh. energy   : ",F20.8," Ha")') EtbBand
    write(*,*)

    if (hasEnergyAnalysis) then
       write(*,'(1x,"Crystal field energy  : ",F20.8," Ha")') Ecf
       write(*,'(1x,"------------------------------------------------")')
       EtbBond = EtbBand - Ecf
       write(*,'(1x,"TB-Bond coh. energy   : ",F20.8," Ha")') EtbBond
       write(*,*)
       write(*,'(1x,"Promotion energy      : ",F20.8," Ha")') Eprom
       write(*,'(1x,"Polarization energy   : ",F20.8," Ha")') Epolar
       write(*,'(1x,"Covalent bond energy  : ",F20.8," Ha")') Ecov
       write(*,'(1x,"E_cov + E_polar       : ",F20.8," Ha")') Echeck
       write(*,'(1x,"------------------------------------------------")')
       EtbBond = Eprom + Ecov + Epolar + Epairpot + Esmear
       write(*,'(1x,"prom+cov+polar+pot    : ",F20.8," Ha")') EtbBond
       write(*,*)
    end if
    
    write(*,*)

  end subroutine tbi_print_energies

  !--------------------------------------------------------------------!
  !                          print out forces                          !
  !--------------------------------------------------------------------!

  subroutine tbi_print_forces()

    implicit none

    integer                        :: iatom, iat, itype
    double precision, dimension(3) :: fsum

!DEBUG
!    if (.not. hasForces) return
    if (.not. doForces) return
!END DEBUG

    write(*,*) 'Atomic forces'
    write(*,*) '============='
    write(*,*)
    write(*,'(4x,8x,A,4x)') 'cartesian force components [Ha/Bohr]'
    write(*,*)
    write(*,'(4x,3(A12,4x))') 'f_x', 'f_y', 'f_z'
    iat = 0
    do itype = 1, nTypes
    do iatom = 1, nAtomsOfType(itype)
       iat = iat + 1
       write(*,'(1x,A2,1x,3(F15.8,1x))') &
            atomTypeName(itype), forCart(1:3,iat)
    end do
    write(*,*)
    end do

    ! check if forces add up to zero:
    fsum(1) = sum(forCart(1,:))
    fsum(2) = sum(forCart(2,:))
    fsum(3) = sum(forCart(3,:))
    if (any(abs(fsum) > 1.0d-6)) then
       write(0,*) "Warning: force sum rule violated."
       write(*,'(1x,A4,1x,3(F13.8,3x),"!!!")') 'sum:', fsum(:)
       write(*,*)
    end if
    write(*,*)

  end subroutine tbi_print_forces

  !--------------------------------------------------------------------!
  !                print eigenvalues (for each k-point)                !
  !--------------------------------------------------------------------!

  subroutine tbi_print_eigenvalues()

    implicit none

    integer :: ikpt, iev

    if (.not. hasEigenvalues) return

    write(*,*) 'Eigenvalues (in eV)'
    write(*,*) '==================='
    write(*,*)
    do ikpt = 1, nkpts
       write(*,'(1x,I5,1x)', advance='no') ikpt
       do iev = 1, nev
          write(*,'(1x,F10.4)', advance='no') eval(iev,ikpt)*Ha2eV
          if((mod(iev,6) == 0) .and. (iev<nev)) then
             write(*,*)
             write(*,'(7x)', advance='no')
          end if
       end do
       write(*,*)
    end do
    write(*,*)
    write(*,*)

  end subroutine tbi_print_eigenvalues

  !--------------------------------------------------------------------!
  !                    print out occupation numbers                    !
  !--------------------------------------------------------------------!

  subroutine tbi_print_occupation()

    implicit none

    integer          :: ikpt, iev

    if (.not. hasEnergies) return

    write(*,*) 'Occupation numbers'
    write(*,*) '=================='
    write(*,*)
    do ikpt = 1, nkpts
       write(*,'(1x,I5,1x)', advance='no') ikpt
       do iev = 1, nev
          write(*, '(1x,F5.3)', advance='no') fOccup(iev,ikpt)
          if ((mod(iev,11)==0) .and. (iev<nev)) then
             write(*,*)
             write(*,'(7x)',advance='no')
          end if
       end do
       write(*,*)
    end do
    write(*,*)
    write(*,*)
    
  end subroutine tbi_print_occupation

  !--------------------------------------------------------------------!
  !                       print Mulliken charges                       !
  !--------------------------------------------------------------------!

  subroutine tbi_print_mulliken()

    implicit none


    integer          :: norbs
    integer          :: itype, iatom, il, l, m, iao
    double precision :: q, q_atom, Q_sum

    if (.not. hasMulliken) return

    write(*,*) 'Mulliken analysis'
    write(*,*) '================='
    write(*,*)
    write(*,'(1x,10x,A4,A8,A4,A4,2x,A6,2x,A6)') &
         'type', 'atom', 'l', 'm', 'Q', 'delQ'
    write(*,*)    

    iao   = 0
    Q_sum = 0.0d0
    do itype = 1, nTypes
    do iatom = 1, nAtomsOfType(itype)
    do il = 1, nl(itype)
       l = l_of_il(il, itype)
       norbs = 2*l + 1
       q_atom = Qatoms(il, itype)/dble(norbs)
       do m = -l, l
          iao = iao + 1         
          ! Mulliken charges:
          q = Qmull(iao)
          Q_sum = Q_sum + q
          write(*,'(1x,I10,I4,I8,I4,I4,2x,F6.3,2x,F6.3)') &
               iao, itype, iatom, l, m, q, q - q_atom
       end do
    end do
    end do
    end do

    write(*,'(1x,46("-"))')
    write(*,'(1x,"total charge : ",11x,F12.3," e")') Q_sum
    write(*,*)
    write(*,*)

  end subroutine tbi_print_mulliken


  !============================= PRIVATE ==============================!


  !--------------------------------------------------------------------!
  !           initialize the calculation and allocate memory           !
  !--------------------------------------------------------------------!

  subroutine tbi_init_calc()

    implicit none

    ! list of k-points:
    nkpts   = kptMultiple(1)*kptMultiple(2)*kptMultiple(3)
    allocate(kpt(3,nkpts), wkpt(nkpts))
    memSize = memSize + 4*nkpts*2
    if (useSymmetry) then
       call sym_ir_kpoint_set(nAtoms, latticeVec, cooLattice, atomType, &
                              kptMultiple, kptShift, nkpts, kpt, wkpt)
    else
       call pbc_MP_kpoint_set(kptMultiple, kptShift, nkpts, kpt, wkpt)
    end if
    
    ! initialize matrices:
    if ((nkpts == 1) .and. all(kpt(1:3,1) == 0.0d0) ) then
       ! gamma point calculation:
       call mat_init(nTypes, nAtoms, atomType, nl, NLMAX, l_of_il, &
                     tbp_S, 'real')
       nev = mat_get_nOrbitals()
       allocate(eval(nev, nkpts), evec_d(nev,nev))
       memSize =  memSize + nev*nkpts*2 + nev*nev*2
       gammaPoint = .true.
    else
       ! k-point sampling:
       call mat_init(nTypes, nAtoms, atomType, nl, NLMAX, l_of_il, &
                     tbp_S, 'complex')
       nev = mat_get_nOrbitals()
       allocate(eval(nev, nkpts), evec_z(nev,nev))
       memSize =  memSize + 2*nev*nkpts + nev*nev*4
       gammaPoint = .false.
    end if

    ! number of orbitals:
    ! (for minimal basis TB this is the same as the number of eigenvalues)
    nao = nev

    d_max = pbc_d_max_origin(nAtoms, cooLattice, latticeVec)
    r_max = max(mat_r_max, pot_r_max)
  
    ! translation vectors within the interaction radius:
    nTvecs  = pbc_number_of_tvecs(r_max, d_max, latticeVec)
    allocate(Tvec(3,nTvecs), TvecLen(nTvecs))
    memSize = memSize + 3*nTvecs + 2*nTvecs
    call pbc_compute_tvecs(r_max, d_max, latticeVec, nTvecs, Tvec, Tveclen)
    
    ! memory for atomic charges, real-space on-site elements, Mulliken charges:
    allocate(Qatoms(NLMAX, nTypes), H_ii(nev))
    memSize = memSize + NLMAX*nTypes + 2*nev

    ! memory for forces:
    if (doForces) then
       allocate(forCart(3,nAtoms))
       memSize = memSize + 3*nAtoms*2
       forCart(:,:) = 0.0d0
    end if

    ! memory for occupation numbers:
    allocate(fOccup(nev,nkpts))
    memSize = memSize + nev*nkpts*2

  end subroutine tbi_init_calc

  !--------------------------------------------------------------------!
  !                  Fermi level and band occupations                  !
  !--------------------------------------------------------------------!

  subroutine tbi_fermi_level()

    implicit none

    integer          :: iev, ikpt
    double precision :: ev, esum

    nOccup = 0
    esum   = 0.0d0

    if (readOccup) then
       flevel = eval(1,1)
       ! use given occupation numbers:
       u_occ = io_unit()
       open(u_occ, file=trim(occupFile), status='old', action='read')
       rewind(u_occ)
       do ikpt = 1, nkpts
          do iev  = 1, nev
             ev = eval(iev,ikpt)
             read(u_occ, *) fOccup(iev,ikpt)
             esum = esum + wkpt(ikpt)*fOccup(iev,ikpt)
             if (fOccup(iev,ikpt) > 0.0d0) then
                flevel = max(flevel,ev)
                nOccup = max(nOccup,iev)
             end if
          end do
       end do
       close(u_occ)
       fermIter = 0
       Esmear   = 0.0d0
    else
       call nrg_fermi_level(nElecs, nkpts, wkpt, nev, eval, metallic, &
            smearWidth, flevel, Esmear, fermIter)
       do ikpt = 1, nkpts
          do iev  = 1, nev
             ev = eval(iev,ikpt)
             fOccup(iev,ikpt)  = 2.0d0*nrg_occupancy(ev, flevel, &
                                 metallic, smearWidth)
             esum = esum + wkpt(ikpt)*fOccup(iev,ikpt)
             if (fOccup(iev,ikpt) > 0.0d0) nOccup = max(nOccup,iev)
          end do
       end do
    end if

    ! check number of electrons:
    if (abs(esum - dble(nElecs)) > 1.0d-5) then
       write(0,*) "Error: wrong number of electrons in `tbinter'."
       stop
    end if

    if (doTiming) call tng_timing('Fermi level done.')

  end subroutine tbi_fermi_level

  !--------------------------------------------------------------------!
  !                              band sum                              !
  !--------------------------------------------------------------------!

  subroutine tbi_band_sum()

    implicit none

    integer :: iev, ikpt

    Eband = 0.0d0
    do ikpt = 1, nkpts
       do iev  = 1, nev
          Eband = Eband + wkpt(ikpt)*fOccup(iev,ikpt)*eval(iev,ikpt)
       end do
    end do

    if (doTiming) call tng_timing('Band sum done.')

  end subroutine tbi_band_sum

  !--------------------------------------------------------------------!
  !                    save density matrix to file                     !
  !--------------------------------------------------------------------!

  subroutine tbi_save_density_matrix()

    implicit none

    integer :: ikpt, iev, iao

    double precision, dimension(:,:), allocatable :: rho_d, work_d
    complex(kind=DP), dimension(:,:), allocatable :: rho_z, work_z

    u_rho = io_unit()
    open(u_rho, file=trim(rhofile), status='replace', action='write', &
                form='unformatted')
    write(u_rho) gammaPoint
    write(u_rho) nkpts, nev

    u_eig = io_unit()
    open(u_eig, file=trim(eigfile), status='old', action='read', &
                form='unformatted')
    read(u_eig)
    read(u_eig)

    gamma : if (gammaPoint) then

       allocate(work_d(nev,nev), rho_d(nev,nev))
       memSize = memSize + 2*nev*nev*2

       read(u_eig)
       read(u_eig) 
       read(u_eig) evec_d(1:nev,1:nev)

       do iev = 1, nev
          do iao = 1, nev
             work_d(iao, iev) = fOccup(iev,1)*evec_d(iao, iev)
          end do
       end do
       rho_d(:,:) = matmul(work_d, transpose(evec_d))
       write(u_rho) kpt(1:3,1), wkpt(1)
       write(u_rho) rho_d(1:nev, 1:nev)

       deallocate(work_d, rho_d)
       memSize = memSize - 2*nev*nev*2

    else

       allocate(work_z(nev,nev), rho_z(nev,nev))
       memSize = memSize + 2*nev*nev*4

       kpoints : do ikpt = 1, nkpts

          read(u_eig)
          read(u_eig) 
          read(u_eig) evec_z(1:nev,1:nev)

          do iev = 1, nev
             do iao = 1, nev
                work_z(iao, iev) = fOccup(iev,ikpt)*evec_z(iao, iev)
             end do
          end do
          rho_z(:,:) = matmul(work_z, transpose(conjg(evec_z)))
          write(u_rho) kpt(1:3,ikpt), wkpt(ikpt)
          write(u_rho) rho_z(1:nev, 1:nev)

       end do kpoints

       deallocate(work_z, rho_z)
       memSize = memSize - 2*nev*nev*4

    end if gamma

    close(u_eig)
    close(u_rho)

    hasDensityMatrix = .true.

  end subroutine tbi_save_density_matrix

  !--------------------------------------------------------------------!
  !                        energy partitioning                         !
  !     Bester and Faehnle, J. Phys.: Cond. Matter 13 (2001) 11541     !
  !--------------------------------------------------------------------!

  subroutine tbi_energy_partitioning()

    implicit none

    integer :: iao, ikpt

    double precision, dimension(:,:), allocatable :: rho_d, pop_d
    complex(kind=DP), dimension(:,:), allocatable :: rho_z, pop_z
    
    if (.not. allocated(Qmull)) then
       allocate(Qmull(nev))
       memSize = memSize + nev*2
    end if

    Qmull(:) = 0.0d0
    Echeck   = 0.0d0
    Eprom    = 0.0d0
    Ecf      = 0.0d0
    Epolar   = 0.0d0

    u_rho = io_unit()
    open(u_rho, file=trim(rhofile), status='old', action='read', &
                form='unformatted')
    read(u_rho)
    read(u_rho)

    gamma : if (gammaPoint) then

       allocate(rho_d(nev,nev), pop_d(nev,nev))
       memSize = memSize + 2*nev*nev*2

       ikpt = 1

       call mat_setup_d(nTypes, nAtoms, atomType, mat_aoNum, cooLattice, &
                       latticeVec, nTvecs, Tvec, nl, NLMAX, l_of_il,     &
                       mat_nOrbsMax, nev, mat_H_0, mat_S_0, mat_H_d, mat_S_d)
       if (doTiming) call tng_timing2('Matrix set-up done.')
       read(u_rho)
       read(u_rho) rho_d(:,:)
       if (doTiming) call tng_timing2('Density matrix read.')
       ! Population matrix:
       pop_d(:,:) = matmul(rho_d, mat_S_d)
       if (doTiming) call tng_timing2('Population matrix done.')
       ! Mulliken charges:
       do iao = 1, nev
          Qmull(iao) = Qmull(iao) + wkpt(ikpt)*pop_d(iao,iao)
       end do

       ! polarization energy and covalent bond energy:
       call tbi_Epolar_Ecov_gamma(rho_d)
          
       deallocate(rho_d, pop_d)
       memSize = memSize - 2*nev*nev*2

    else

       allocate(rho_z(nev,nev), pop_z(nev,nev))
       memSize = memSize + 2*nev*nev*4

       kpoints : do ikpt = 1, nkpts
          read(u_rho)

          call mat_setup_z(nTypes, nAtoms, atomType, mat_aoNum, cooLattice, &
                           latticeVec, nTvecs, Tvec, nl, NLMAX, l_of_il,    &
                           kpt(1:3,ikpt), mat_nOrbsMax, nev, mat_H_0,       &
                           mat_S_0, mat_H_z, mat_S_z)
          if (doTiming) call tng_timing2('Matrix set-up done (k-point ' &
                          // trim(io_adjustl(ikpt)) // ')')
          read(u_rho) rho_z(:,:)
          if (doTiming) call tng_timing2('Density matrix read (k-point ' &
                          // trim(io_adjustl(ikpt)) // ')')
          ! Population matrix for this k-point:
          pop_z(:,:) = matmul(rho_z, mat_S_z)
          if (doTiming) call tng_timing2('Population matrix done (k-point ' &
                          // trim(io_adjustl(ikpt)) // ')')
          ! Mulliken charges:
          do iao = 1, nev
             Qmull(iao) = Qmull(iao) + wkpt(ikpt)*real(pop_z(iao,iao))
          end do

          ! polarization energy and covalent bond energy:
          call tbi_Epolar_Ecov_kpt(ikpt, rho_z)
       end do kpoints

       deallocate(rho_z, pop_z)
       memSize = memSize - 2*nev*nev*4

    end if gamma

    close(u_rho)

    ! Mulliken analysis, promotion energy, crystal field energy:
    call tbi_Eprom_Ecf()

    ! covalent bond energy:
    Ecov = Eband - Eatoms - Eprom - Epolar - Ecf

    hasMulliken = .true.

    if (doTiming) call tng_timing('Energy partitioning done.')

  end subroutine tbi_energy_partitioning
  
  !--------------------------------------------------------------------!
  !               polarization and covalent bond energy                !
  !--------------------------------------------------------------------!

  subroutine tbi_Epolar_Ecov_kpt(ikpt, rho_z)

    implicit none

    integer,                              intent(in) :: ikpt
    complex(kind=DP), dimension(nev,nev), intent(in) :: rho_z

    complex(kind=DP) :: Ecov_z, Epol_z
    double precision :: H_ii_mean
    integer          :: itype1, iat1, iao1, norbs1, iorb1
    integer          :: itype2, iat2, iao2, norbs2, iorb2

    Ecov_z = (0.0d0, 0.0d0)
    Epol_z = (0.0d0, 0.0d0)

    atom1 : do iat1 = 1, nAtoms
       iao1   = mat_aoNum(iat1)
       itype1 = atomType(iat1)
       norbs1 = mat_nOrbs(itype1)

       ! E_polar:
       do iorb1 = iao1,    iao1 + norbs1 - 1
       do iorb2 = iorb1+1, iao1 + norbs1 - 1
          Epol_z = Epol_z + rho_z(iorb1,iorb2)*mat_H_0(iorb2,iorb1)
          Epol_z = Epol_z + rho_z(iorb2,iorb1)*mat_H_0(iorb1,iorb2)
       end do
       end do

       atom2 : do iat2 = 1, nAtoms
          iao2   = mat_aoNum(iat2)
          itype2 = atomType(iat2)
          norbs2 = mat_nOrbs(itype2)

          do iorb1 = iao1, iao1 + norbs1 - 1
          do iorb2 = iao2, iao2 + norbs2 - 1
             H_ii_mean = 0.5d0*(H_ii(iorb1) + H_ii(iorb2))
             ! E_cov and E_polar:
             Ecov_z = Ecov_z + rho_z(iorb1,iorb2)  &
                      *( mat_H_z(iorb2,iorb1)      &
                       - mat_S_z(iorb2,iorb1)*dcmplx(H_ii_mean,0.0d0) )
          end do
          end do

       end do atom2
    end do atom1

    Echeck = Echeck + wkpt(ikpt)*real(Ecov_z)
    Epolar = Epolar + wkpt(ikpt)*real(Epol_z)

    if (doTiming) call tng_timing2('E_cov done (k-point ' &
         // trim(io_adjustl(ikpt)) // ')')

  end subroutine tbi_Epolar_Ecov_kpt

  !--------------------------------------------------------------------!

  subroutine tbi_Epolar_Ecov_gamma(rho_d)

    implicit none

    double precision, dimension(nev,nev), intent(in) :: rho_d

    double precision :: Ecov_d, Epol_d
    double precision :: H_ii_mean
    integer          :: itype1, iat1, iao1, norbs1, iorb1
    integer          :: itype2, iat2, iao2, norbs2, iorb2

    Ecov_d = 0.0d0
    Epol_d = 0.0d0

    atom1 : do iat1 = 1, nAtoms
       iao1   = mat_aoNum(iat1)
       itype1 = atomType(iat1)
       norbs1 = mat_nOrbs(itype1)

       ! E_polar:
       do iorb1 = iao1,    iao1 + norbs1 - 1
       do iorb2 = iorb1+1, iao1 + norbs1 - 1
          Epol_d = Epol_d + rho_d(iorb1,iorb2)*mat_H_0(iorb2,iorb1)
          Epol_d = Epol_d + rho_d(iorb2,iorb1)*mat_H_0(iorb1,iorb2)
       end do
       end do

       atom2 : do iat2 = 1, nAtoms
          iao2   = mat_aoNum(iat2)
          itype2 = atomType(iat2)
          norbs2 = mat_nOrbs(itype2)

          do iorb1 = iao1, iao1 + norbs1 - 1
          do iorb2 = iao2, iao2 + norbs2 - 1
             H_ii_mean = 0.5d0*(H_ii(iorb1) + H_ii(iorb2))
             ! E_cov and E_polar:
             Ecov_d = Ecov_d + rho_d(iorb1,iorb2)  &
                      *( mat_H_d(iorb2,iorb1)      &
                       - mat_S_d(iorb2,iorb1)*H_ii_mean )
          end do
          end do

       end do atom2
    end do atom1

    Echeck = Echeck + Ecov_d
    Epolar = Epolar + Epol_d

    if (doTiming) call tng_timing2('E_cov done (gamma point)')

  end subroutine tbi_Epolar_Ecov_gamma

  !--------------------------------------------------------------------!
  !              promotion energy / crystal field energy               !
  !--------------------------------------------------------------------!

  subroutine tbi_Eprom_Ecf()

    implicit none

    integer          :: norbs
    integer          :: itype, iatom, il, l, m, iao
    double precision :: q, q_atom

    Eprom  = 0.0d0
    Ecf    = 0.0d0

    iao   = 0
    do itype = 1, nTypes
    do iatom = 1, nAtomsOfType(itype)
    do il = 1, nl(itype)
       l = l_of_il(il, itype)
       norbs = 2*l + 1
       q_atom = Qatoms(il, itype)/dble(norbs)
       do m = -l, l
          iao = iao + 1
          q = Qmull(iao)
          Eprom = Eprom + (q - q_atom)*tbp_level(il,itype)
!DEBUG
!          write(0,*) "***** debugging version *****"
!          Ecf   = Ecf + abs((q - q_atom)*(H_ii(1) - tbp_level(1,itype)))
!          H: Ecf   = Ecf + (q - q_atom)*(H_ii(iao) - tbp_level(il,itype))
!          G: Ecf   = Ecf - q_atom*H_ii(1)
!          F: Ecf   = Ecf - q*tbp_level(il,itype)
!          E: Ecf   = Ecf - q_atom*H_ii(iao)
!          D: Ecf   = Ecf - q*H_ii(iao)
!          C: Ecf   = Ecf - q_atom*tbp_level(il,itype)
!          B: Ecf   = Ecf   + (H_ii(1) - tbp_level(1,itype))
!
! s-level ref  Ecf = Ecf + q_atom*(H_ii(1) - tbp_level(1,itype))
! no charge    Ecf = Ecf + q_atom*(H_ii(iao) - tbp_level(il,itype))
! no E_prom    Ecf = Ecf + (q*H_ii(iao) - q_atom*tbp_level(il,itype))
!
!END DEBUG
          Ecf   = Ecf   + q*(H_ii(iao) - tbp_level(il,itype))
       end do
    end do
    end do
    end do

    if (doTiming) call tng_timing2('Promotion energy and CF energy done.')

  end subroutine tbi_Eprom_Ecf

  !--------------------------------------------------------------------!
  !                          calculate forces                          !
  !--------------------------------------------------------------------!

  subroutine tbi_calc_forces_kpt()

    implicit none

    integer                            :: ikpt
    integer                            :: iat1, iao1, iao1f, itype1, nl1, norbs1
    integer                            :: iat2, iao2, iao2f, itype2, nl2, norbs2
    integer,          dimension(NLMAX) :: lvec1, lvec2
    double precision, dimension(3)     :: coo1, coo2
    double precision, dimension(3)     :: force
    integer                            :: i1, i2
    integer                            :: dim

    complex(kind=DP), dimension(mat_nOrbsMax,mat_nOrbsMax,3) :: dH, dS
    double precision, dimension(mat_nOrbsMax,mat_nOrbsMax,3) :: dO1, dO2

    complex(kind=DP), dimension(mat_nOrbsMax,mat_nOrbsMax)   :: rho11, rho22, rho12, erho12
    complex(kind=DP), dimension(:,:),            allocatable :: work

    allocate(work(mat_nOrbsMax,nOccup))
    memSize = memSize + 4*mat_nOrbsMax*nOccup

    u_eig = io_unit()
    open(u_eig, file=trim(eigfile), status='old', action='read', &
                form='unformatted')
    read(u_eig)
    read(u_eig)

    kpoints : do ikpt = 1, nkpts

       read(u_eig)
       read(u_eig) 
       read(u_eig) evec_z(1:nev,1:nev)

       atom1 : do iat1 = 1, nAtoms
          itype1        = atomType(iat1)
          coo1          = cooLattice(1:3,iat1)
          nl1           = nl(itype1)
          lvec1(1:nl1)  = l_of_il(1:nl1, itype1)
          norbs1        = mat_nOrbs(itype1)
          iao1          = mat_aoNum(iat1)
          iao1f         = iao1 + norbs1 - 1

          if (tbp_O) then
             ! density matrix of the orbitals of atom 1:
             call tbi_local_rho_ons_kpt(iao1, iao1f, norbs1, mat_nOrbsMax,  &
                                        nOccup, nev, fOccup(1:nOccup,ikpt), &
                                        evec_z, work, rho11)
          end if

          atom2 : do iat2 = iat1+1, nAtoms
             itype2        = atomType(iat2)
             coo2          = cooLattice(1:3,iat2)
             nl2           = nl(itype2)
             lvec2(1:nl2)  = l_of_il(1:nl2, itype2)
             norbs2        = mat_nOrbs(itype2)
             iao2          = mat_aoNum(iat2)
             iao2f         = iao2 + norbs2 - 1
             
             call mat_delem_full_z(itype1, coo1, nl1, lvec1,             &
                                   itype2, coo2, nl2, lvec2,             &
                                   latticeVec, recLattVec, nTvecs, Tvec, &
                                   kpt(:,ikpt), mat_nOrbsMax, tbp_O,     &
                                   dH, dS, dO1, dO2)
             
             ! density matrix and energy weighted density matrix for the
             ! orbitals of atoms 1 and 2:
             call tbi_local_rho_erho_kpt(iao2, iao2f, norbs2, iao1, norbs1,          &
                                         mat_nOrbsMax, nOccup, nev,                  &
                                         fOccup(1:nOccup,ikpt), eval(1:nOccup,ikpt), &
                                         evec_z, work, rho12, erho12)
             force(:) = 0.0d0
             do dim = 1, 3
             do i1 = 1, norbs1
             do i2 = 1, norbs2
                force(dim) = force(dim) + real( rho12(i1,i2)*dH(i2,i1,dim)) &
                                        - real(erho12(i1,i2)*dS(i2,i1,dim))
             end do
             end do
             end do
             force(1:3) = wkpt(ikpt)*2.0d0*force(1:3)
             forCart(1:3,iat1) = forCart(1:3,iat1) - force(1:3)
             forCart(1:3,iat2) = forCart(1:3,iat2) + force(1:3)

             ! forces from the on-site parametrization:
             if (tbp_O) then 

                ! on-site elements of the first atom:
                force(:) = 0.0d0
                do dim = 1, 3
                do i1 = 1, norbs1
                do i2 = 1, norbs1
                   force(dim) = force(dim) + real(rho11(i1, i2)*dO1(i2,i1,dim))
                end do
                end do
                end do
                force(1:3) = wkpt(ikpt)*force(1:3)
                forCart(1:3,iat1) = forCart(1:3,iat1) - force(1:3)
                forCart(1:3,iat2) = forCart(1:3,iat2) + force(1:3)

                ! density matrix of the orbitals of atom 2:
                call tbi_local_rho_ons_kpt(iao2, iao2f, norbs2, mat_nOrbsMax,  &
                                           nOccup, nev, fOccup(1:nOccup,ikpt), &
                                           evec_z, work, rho22)

                ! on-site elements of the second atom:
                force(:) = 0.0d0
                do dim = 1, 3
                do i1 = 1, norbs2
                do i2 = 1, norbs2
                   force(dim) = force(dim) + real(rho22(i1,i2)*dO2(i2,i1,dim))
                end do
                end do
                end do
                force(1:3) = wkpt(ikpt)*force(1:3)
                forCart(1:3,iat1) = forCart(1:3,iat1) - force(1:3)
                forCart(1:3,iat2) = forCart(1:3,iat2) + force(1:3)

             end if         

          end do atom2
       end do atom1

    end do kpoints

    close(u_eig)

    deallocate(work)
    memSize = memSize - 4*mat_nOrbsMax*nOccup

  end subroutine tbi_calc_forces_kpt

  !--------------------------------------------------------------------!

  subroutine tbi_calc_forces_gamma()

    implicit none

    integer                            :: iat1, iao1, iao1f, itype1, nl1, norbs1
    integer                            :: iat2, iao2, iao2f, itype2, nl2, norbs2
    integer,          dimension(NLMAX) :: lvec1, lvec2
    double precision, dimension(3)     :: coo1, coo2
    double precision, dimension(3)     :: force
    integer                            :: i1, i2
    integer                            :: dim

    double precision, dimension(mat_nOrbsMax,mat_nOrbsMax,3) :: dH, dS, dO1, dO2

    double precision, dimension(mat_nOrbsMax,mat_nOrbsMax)   :: rho11, rho22, rho12, erho12
    double precision, dimension(:,:),            allocatable :: work

    allocate(work(mat_nOrbsMax,nOccup))
    memSize = memSize + 2*mat_nOrbsMax*nOccup

    !------------------------------------------------------------------!

    u_eig = io_unit()
    open(u_eig, file=trim(eigfile), status='old', action='read', &
                form='unformatted')
    read(u_eig)
    read(u_eig)

    read(u_eig)
    read(u_eig) 
    read(u_eig) evec_d(1:nev,1:nev)
    close(u_eig)

    !------------------------------------------------------------------!

    atom1 : do iat1 = 1, nAtoms
       itype1        = atomType(iat1)
       coo1          = cooLattice(1:3,iat1)
       nl1           = nl(itype1)
       lvec1(1:nl1)  = l_of_il(1:nl1, itype1)
       norbs1        = mat_nOrbs(itype1)
       iao1          = mat_aoNum(iat1)
       iao1f         = iao1 + norbs1 - 1

       if (tbp_O) then
          ! density matrix of the orbitals of atom 1:
          call tbi_local_rho_ons(iao1, iao1f, norbs1, mat_nOrbsMax, &
                                 nOccup, nev, fOccup(1:nOccup,1),   &
                                 evec_d, work, rho11)
       end if

       atom2 : do iat2 = iat1+1, nAtoms
          itype2        = atomType(iat2)
          coo2          = cooLattice(1:3,iat2)
          nl2           = nl(itype2)
          lvec2(1:nl2)  = l_of_il(1:nl2, itype2)
          norbs2        = mat_nOrbs(itype2)
          iao2          = mat_aoNum(iat2)
          iao2f         = iao2 + norbs2 - 1
          
          call mat_delem_full_d(itype1, coo1, nl1, lvec1,             &
                                itype2, coo2, nl2, lvec2,             &
                                latticeVec, recLattVec, nTvecs, Tvec, &
                                mat_nOrbsMax, tbp_O, dH, dS, dO1, dO2)

          ! density matrix and energy weighted density matrix for the
          ! orbitals of atoms 1 and 2:
          call tbi_local_rho_erho(iao1, iao1f, norbs1, iao2, norbs2,    &
                                  mat_nOrbsMax, nOccup, nev,            &
                                  fOccup(1:nOccup,1), eval(1:nOccup,1), &
                                  evec_d, work, rho12, erho12)

          force(:) = 0.0d0
          do dim = 1, 3
          do i1 = 1, norbs1
          do i2 = 1, norbs2
             ! gamma point -> symmetric density matrix
             ! (indices of rho swapped for efficiency)
             force(dim) = force(dim) +  rho12(i2,i1)*dH(i2,i1,dim) &
                                     - erho12(i2,i1)*dS(i2,i1,dim)
          end do
          end do
          end do
          ! factor two for the complex conjugated derivative:
          ! d/dR_i <mu|H|nu> = d/dR_i conjg(<nu|H|mu>) = conjg(d/dR_i <nu|H|mu>)
          force(1:3) = 2.0d0*force(1:3)
          forCart(1:3,iat1) = forCart(1:3,iat1) - force(1:3)
          forCart(1:3,iat2) = forCart(1:3,iat2) + force(1:3)

          ! forces from the on-site parametrization:
          if (tbp_O) then 

             ! on-site elements of the first atom:
             force(:) = 0.0d0
             do dim = 1, 3
             do i1 = 1, norbs1
             do i2 = 1, norbs1
                force(dim) = force(dim) + rho11(i2,i1)*dO1(i2,i1,dim)
             end do
             end do
             end do
             forCart(1:3,iat1) = forCart(1:3,iat1) - force(1:3)
             forCart(1:3,iat2) = forCart(1:3,iat2) + force(1:3)

             ! density matrix for orbitals of atom 2:
             call tbi_local_rho_ons(iao2, iao2f, norbs2, mat_nOrbsMax, &
                                    nOccup, nev, fOccup(1:nOccup,1),   &
                                    evec_d, work, rho22)

             ! on-site elements of the second atom:
             force(:) = 0.0d0
             do dim = 1, 3
             do i1 = 1, norbs2
             do i2 = 1, norbs2
                force(dim) = force(dim) + rho22(i2,i1)*dO2(i2,i1,dim)
             end do
             end do
             end do
             forCart(1:3,iat1) = forCart(1:3,iat1) - force(1:3)
             forCart(1:3,iat2) = forCart(1:3,iat2) + force(1:3)

          end if         

       end do atom2
    end do atom1

    deallocate(work)
    memSize = memSize - 2*mat_nOrbsMax*nOccup

  end subroutine tbi_calc_forces_gamma

  !--------------------------------------------------------------------!
  !         calculate two-atoms block of the density matrices          !
  !--------------------------------------------------------------------!

  subroutine tbi_local_rho_erho(iao1, iao1f, no1, iao2, no2, nomax, &
                                nocc, nev, focc, eval, evec, work,  &
                                rho, erho)

    implicit none

    integer,                                  intent(in)    :: iao1, iao1f, no1
    integer,                                  intent(in)    :: iao2, no2
    integer,                                  intent(in)    :: nomax, nocc, nev
    double precision, dimension(nocc),        intent(in)    :: focc
    double precision, dimension(nocc),        intent(in)    :: eval
    double precision, dimension(nev,nev),     intent(in)    :: evec
    double precision, dimension(nomax,nocc),  intent(inout) :: work
    double precision, dimension(nomax,nomax), intent(out)   :: rho, erho

    integer :: iev, i2

    ! density matrix of the orbitals of atoms 1 and 2:
    do iev = 1, nocc
    do i2 = 1, no2
       work(i2, iev) = focc(iev)*evec(i2+iao2-1, iev)
    end do
    end do
    rho(1:no2,1:no1) &
    = matmul(work(1:no2,1:nocc), transpose(evec(iao1:iao1f,1:nocc)))

    ! energy weighted density matrix of the orbitals of atoms 1 and 2:
    do iev = 1, nocc
    do i2 = 1, no2
       work(i2, iev) = eval(iev)*work(i2, iev)
    end do
    end do
    erho(1:no2,1:no1) &
    = matmul(work(1:no2,1:nocc), transpose(evec(iao1:iao1f,1:nocc)))
    
  end subroutine tbi_local_rho_erho

  !--------------------------------------------------------------------!

  subroutine tbi_local_rho_erho_kpt(iao1, iao1f, no1, iao2, no2, nomax, &
                                nocc, nev, focc, eval, evec, work,  &
                                rho, erho)

    implicit none

    integer,                                  intent(in)    :: iao1, iao1f, no1
    integer,                                  intent(in)    :: iao2, no2
    integer,                                  intent(in)    :: nomax, nocc, nev
    double precision, dimension(nocc),        intent(in)    :: focc
    double precision, dimension(nocc),        intent(in)    :: eval
    complex(kind=DP), dimension(nev,nev),     intent(in)    :: evec
    complex(kind=DP), dimension(nomax,nocc),  intent(inout) :: work
    complex(kind=DP), dimension(nomax,nomax), intent(out)   :: rho, erho

    integer :: iev, i2

    ! density matrix of the orbitals of atoms 1 and 2:
    do iev = 1, nocc
    do i2 = 1, no2
       work(i2, iev) = focc(iev)*evec(i2+iao2-1, iev)
    end do
    end do
    rho(1:no2,1:no1) &
    = matmul(work(1:no2,1:nocc), transpose(conjg(evec(iao1:iao1f,1:nocc))))

    ! energy weighted density matrix of the orbitals of atoms 1 and 2:
    do iev = 1, nocc
    do i2 = 1, no2
       work(i2, iev) = eval(iev)*work(i2, iev)
    end do
    end do
    erho(1:no2,1:no1) &
    = matmul(work(1:no2,1:nocc), transpose(conjg(evec(iao1:iao1f,1:nocc))))
    
  end subroutine tbi_local_rho_erho_kpt

  !--------------------------------------------------------------------!
  !           calculate one-atom block of the density matrix           !
  !--------------------------------------------------------------------!

  subroutine tbi_local_rho_ons(iao1, iao1f, no1, nomax, nocc, nev, &
                               focc, evec, work, rho)

    implicit none

    integer,                                  intent(in)    :: iao1, iao1f, no1
    integer,                                  intent(in)    :: nomax, nocc, nev
    double precision, dimension(nocc),        intent(in)    :: focc
    double precision, dimension(nev,nev),     intent(in)    :: evec
    double precision, dimension(nomax,nocc),  intent(inout) :: work
    double precision, dimension(nomax,nomax), intent(out)   :: rho

    integer :: iev, i1

    ! density matrix of the orbitals of atoms 1 and 2:
    do iev = 1, nocc
    do i1 = 1, no1
       work(i1, iev) = focc(iev)*evec(i1+iao1-1, iev)
    end do
    end do
    rho(1:no1,1:no1) &
    = matmul(work(1:no1,1:nocc), transpose(evec(iao1:iao1f,1:nocc)))

  end subroutine tbi_local_rho_ons

  !--------------------------------------------------------------------!

  subroutine tbi_local_rho_ons_kpt(iao1, iao1f, no1, nomax, nocc, nev, &
                                   focc, evec, work, rho)

    implicit none

    integer,                                  intent(in)    :: iao1, iao1f, no1
    integer,                                  intent(in)    :: nomax, nocc, nev
    double precision, dimension(nocc),        intent(in)    :: focc
    complex(kind=DP), dimension(nev,nev),     intent(in)    :: evec
    complex(kind=DP), dimension(nomax,nocc),  intent(inout) :: work
    complex(kind=DP), dimension(nomax,nomax), intent(out)   :: rho

    integer :: iev, i1

    ! density matrix of the orbitals of atoms 1 and 2:
    do iev = 1, nocc
    do i1 = 1, no1
       work(i1, iev) = focc(iev)*evec(i1+iao1-1, iev)
    end do
    end do
    rho(1:no1,1:no1) &
    = matmul(work(1:no1,1:nocc), transpose(conjg(evec(iao1:iao1f,1:nocc))))

  end subroutine tbi_local_rho_ons_kpt


end module tbinter
