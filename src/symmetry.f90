module symmetry

  !--------------------------------------------------------------------!
  ! This module contains routines to determine the crystal symmetry    !
  ! and for the search of irreducible k-points.                        !
  !                                                                    !
  ! The SPGlib library is used to determine the symmetry operations.   !
  !                                                                    !
  !--------------------------------------------------------------------!
  ! 2011-09-26 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use spglib, only: spg_find_primitive,         &
                    spg_get_international,      &
                    spg_get_ir_reciprocal_mesh, &
                    spg_get_multiplicity,       &
                    spg_get_schoenflies,        &
                    spg_get_smallest_lattice,   &
                    spg_get_symmetry

  implicit none

  public :: sym_ir_kpoint_set, &
            sym_get_symmetry,  &
            sym_print_kpoints

contains

  !--------------------------------------------------------------------!
  !                    set of irreducible k-points                     !
  !--------------------------------------------------------------------!

  subroutine sym_ir_kpoint_set(natoms, lattice, positions, types, &
                               n, s, nkpts, kpt, w)

    implicit none

    integer,                                intent(in)    :: natoms
    double precision, dimension(3,3),       intent(in)    :: lattice
    double precision, dimension(3,natoms),  intent(in)    :: positions
    integer,          dimension(natoms),    intent(in)    :: types
    integer,          dimension(3),         intent(in)    :: n
    double precision, dimension(3),         intent(in)    :: s
    integer,                                intent(inout) :: nkpts
    double precision, dimension(1:3,nkpts), intent(out)   :: kpt
    double precision, dimension(nkpts),     intent(out)   :: w

    integer,          dimension(3,n(1)*n(2)*n(3)) :: grid_point
    integer,          dimension(n(1)*n(2)*n(3))   :: map
    integer,          dimension(3)                :: is_shift
    integer                                       :: is_time_reversal

    double precision,                   parameter :: symprec = 1.0d-8

    integer :: cnt, ikpt, ikpt1, ikpt2, nkpt_tot

    nkpt_tot = n(1)*n(2)*n(3)
    if (nkpts < (nkpt_tot/2 + 1)) then
       write(0,*) "Error: insufficient memory for k-points."
       stop
    end if

    is_time_reversal = 1
    is_shift(1:3)    = 0
    if (s(1) /= 0.0d0) is_shift(1) = 1
    if (s(2) /= 0.0d0) is_shift(2) = 1
    if (s(3) /= 0.0d0) is_shift(3) = 1

    call spg_get_ir_reciprocal_mesh(nkpts, grid_point, map, n,      &
         is_shift, is_time_reversal, transpose(lattice), positions, &
         types, natoms, symprec)

    ! irreducible k-points:
    ikpt = 0
    do ikpt1 = 1, nkpt_tot
       if (map(ikpt1) == ikpt1-1) then
          cnt  = 0
          ! determine multiplicity:
          do ikpt2 = 1, nkpt_tot
             if (map(ikpt2) == map(ikpt1)) cnt = cnt + 1
          end do
          ! new irreducible k-point identified:
          ikpt = ikpt + 1
          if (ikpt > nkpts) then
             write(0,*) "Error: invalid number of k-points."
             stop
          end if
          kpt(1, ikpt) = dble(grid_point(1,ikpt1))/dble(n(1))
          kpt(2, ikpt) = dble(grid_point(2,ikpt1))/dble(n(2))
          kpt(3, ikpt) = dble(grid_point(3,ikpt1))/dble(n(3))
          w(ikpt)      = dble(cnt)/dble(nkpt_tot)
       end if
    end do

  end subroutine sym_ir_kpoint_set

  !--------------------------------------------------------------------!
  !          get symmetry operations and print info to stdout          !
  !--------------------------------------------------------------------!

  subroutine sym_get_symmetry(num_atom, lattice, positions, atom_types)

    implicit none

    integer,                                 intent(in) :: num_atom
    double precision, dimension(3,3),        intent(in) :: lattice
    double precision, dimension(3,num_atom), intent(in) :: positions
    integer,          dimension(num_atom),   intent(in) :: atom_types

    integer,          dimension(:,:,:), allocatable     :: rot
    double precision, dimension(:,:),   allocatable     :: trans
    integer,          dimension(3,3)                    :: rotmat

    integer                                             :: spacegroup
    character(len=20)                                   :: international
    character(len=9)                                    :: schoenflies

    integer                                             :: num_atom_p
    double precision, dimension(3,3)                    :: lattice_p
    double precision, dimension(3,3)                    :: lattice_s
    double precision, dimension(3,num_atom)             :: positions_p
    integer,          dimension(num_atom)               :: atom_types_p

    integer,                                  parameter :: max_num_sym = 1024
    double precision,                         parameter :: symprec     = 1.d-8

    integer :: isym, nsym
    integer :: iat

    !---------------------!
    ! primitive unit cell !
    !---------------------!

    num_atom_p       = num_atom
    lattice_p(:,:)   = lattice(:,:)
    positions_p(:,:) = positions(:,:)
    atom_types_p(:)  = atom_types(:)

    call spg_find_primitive(lattice_p, positions_p, atom_types_p, &
                            num_atom_p, symprec)
    if (num_atom_p > 0 .and. num_atom_p < num_atom) then
       call spg_get_smallest_lattice(lattice_s, lattice_p, symprec)
       write(*,*) 'primitive lattice vectors:'
       write(*,*)
       write(*,'(7x,3(ES16.8,2x))') lattice_s(1:3,1)
       write(*,'(7x,3(ES16.8,2x))') lattice_s(1:3,2)
       write(*,'(7x,3(ES16.8,2x))') lattice_s(1:3,3)
       write(*,*)
       write(*,*) 'primitive atomic coordinates:'
       write(*,*)
       do iat = 1, num_atom_p
          write(*,'(5x,I2,3(ES16.8,2x))') atom_types_p(iat), &
                                          positions_p(1:3,iat)
       end do
       write(*,*)
    else
       write(*,*) 'The input cell is a primitive unit cell.'
       write(*,*)
    end if

    !-------------!
    ! space group !
    !-------------!

    call spg_get_international(spacegroup, international, lattice, &
                               positions, atom_types, num_atom, symprec)
    call spg_get_schoenflies(spacegroup, schoenflies, lattice, &
                             positions, atom_types, num_atom, symprec)

    write(*,*) 'Spacegroup'
    write(*,*) '=========='
    write(*,*)
    write(*,'(1x,"Spacegroup           : ",I4)') spacegroup
    write(*,'(1x,"Schoenflies symbol   : ",A)') trim(schoenflies)
    write(*,'(1x,"International Tables : ",A)') trim(international)
    write(*,*)

    !-------------------------------!
    ! number of symmetry operations !
    !-------------------------------!

    call spg_get_multiplicity(nsym, transpose(lattice), positions, &
                              atom_types, num_atom, symprec)

    if (nsym == 0) then
       write(*,*) 'No symmetry found.'
       return
    else
       write(*,'(1x,"Total number of symmetry operations: ",I3)') nsym
    end if
    write(*,*)

    allocate(rot(3,3,nsym), trans(3,nsym))

    !---------------------!
    ! symmetry operations !
    !---------------------!
  
    call spg_get_symmetry(nsym, rot, trans, max_num_sym, &
         transpose(lattice), positions, atom_types, num_atom, symprec)

    write(*,*) 'List of symmetry operations'
    write(*,*) '==========================='
    write(*,*)
    write(*,'(5x,A36,3(A6,2x))') '|-------- rotation matrix --------|', &
                                 't_x', 't_y', 't_z'
    do isym = 1, nsym
       rotmat = transpose(rot(1:3,1:3,isym))
       write(*,'(1x,I3,2x,3(I2,2x))', advance='no') isym, rotmat(1:3,1)
       write(*,'(3(I2,2x))', advance='no') rotmat(1:3,2)
       write(*,'(3(I2,2x))', advance='no') rotmat(1:3,3)
       write(*,'(3(F6.3,2x))') trans(1:3,isym)
    end do
    write(*,*)

  end subroutine sym_get_symmetry

  !--------------------------------------------------------------------!
  !                           I/O procedures                           !
  !--------------------------------------------------------------------!

  subroutine sym_print_kpoints(nkpts, kpt, wkpt)

    implicit none

    integer,                              intent(in) :: nkpts
    double precision, dimension(3,nkpts), intent(in) :: kpt
    double precision, dimension(nkpts),   intent(in) :: wkpt

    integer :: ikpt

    write(*,*) 'Irreducible k-points'
    write(*,*) '===================='
    write(*,*)
    write(*,'(1x,"Number of ir. k-points : ",I5)') nkpts
    write(*,*)
    write(*,'(10x,A8,4x,3(A10,2x))') ' weight ', 'k_x   ', 'k_y   ', 'k_z   '
    do ikpt = 1, nkpts
       write(*,'(1x,I5,4x,F8.6,4x,3(F10.8,2x))') ikpt, wkpt(ikpt), &
                                                 kpt(1:3,ikpt)
    end do
    write(*,*)

  end subroutine sym_print_kpoints

end module symmetry
