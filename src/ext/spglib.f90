module spglib

  !--------------------------------------------------------------------!
  ! This module simply contains interfaces to the routines in the      !
  ! SPGLIB symmetry library.  See 'spglib_f.c' for the C headers and   !
  ! spglib.c (of the SPGLIB package) for the implementations.          !
  ! In order to use the symmetry routines one has has to link to       !
  ! 'libsymspg.a/.so', e.g. using something like:                      !
  !                                                                    !
  !   $(LD) -lsymspg <objects> spglib.o spglib_f.o                     !
  !                                                                    !
  ! Where <objects> are the other object files to be linked.           !
  !--------------------------------------------------------------------!
  ! 2011-09-17 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  implicit none

  public :: spg_get_symmetry,                   &
            spg_get_international,              &
            spg_get_schoenflies,                &
            spg_find_primitive,                 &
            spg_refine_cell,                    &
            spg_get_smallest_lattice,           &
            spg_get_multiplicity,               &
            spg_get_max_multiplicity,           &
            spg_get_ir_kpoints,                 &
            spg_get_ir_reciprocal_mesh,         &
            spg_get_stabilized_reciprocal_mesh, &
            spg_get_triplets_reciprocal_mesh

  interface 

     subroutine spg_get_symmetry(nsym, rot, trans, size, lattice, &
                                 position, types, num_atom, symprec)

       implicit none

       integer,                                 intent(out) :: nsym
       integer,                                 intent(in)  :: size
       integer,          dimension(3,3,size),   intent(out) :: rot
       double precision, dimension(3,size),     intent(out) :: trans
       double precision, dimension(3,3),        intent(in)  :: lattice
       integer,                                 intent(in)  :: num_atom
       double precision, dimension(3,num_atom), intent(in)  :: position
       integer,          dimension(num_atom),   intent(in)  :: types
       double precision,                        intent(in)  :: symprec

     end subroutine spg_get_symmetry

     !-----------------------------------------------------------------!

     subroutine spg_get_international(spacegroup, symbol, lattice, &
                                      position, types, num_atom, symprec)

       implicit none

       integer,                                 intent(out) :: spacegroup
       character(len=20),                       intent(out) :: symbol
       double precision, dimension(3,3),        intent(in)  :: lattice
       integer,                                 intent(in)  :: num_atom
       double precision, dimension(3,num_atom), intent(in)  :: position
       integer,          dimension(num_atom),   intent(in)  :: types
       double precision,                        intent(in)  :: symprec

     end subroutine spg_get_international

     !-----------------------------------------------------------------!

     subroutine spg_get_schoenflies(spacegroup, symbol, lattice, &
                                    position, types, num_atom, symprec)

       implicit none

       integer,                                 intent(out) :: spacegroup
       character(len=9),                        intent(out) :: symbol
       double precision, dimension(3,3),        intent(in)  :: lattice
       integer,                                 intent(in)  :: num_atom
       double precision, dimension(3,num_atom), intent(in)  :: position
       integer,          dimension(num_atom),   intent(in)  :: types
       double precision,                        intent(in)  :: symprec

     end subroutine spg_get_schoenflies

     !-----------------------------------------------------------------!

     subroutine spg_find_primitive(lattice, position, types, num_atom, &
                                   symprec)

       implicit none

       double precision, dimension(3,3),        intent(inout) :: lattice
       integer,                                 intent(inout) :: num_atom
       double precision, dimension(3,num_atom), intent(inout) :: position
       integer,          dimension(num_atom),   intent(inout) :: types
       double precision,                        intent(in)    :: symprec

     end subroutine spg_find_primitive

     !-----------------------------------------------------------------!

     subroutine spg_refine_cell(lattice, position, types, num_atom, &
                                symprec)

       implicit none

       double precision, dimension(3,3),          intent(inout) :: lattice
       integer,                                   intent(inout) :: num_atom
       double precision, dimension(3,4*num_atom), intent(inout) :: position
       integer,          dimension(4*num_atom),   intent(inout) :: types
       double precision,                          intent(in)    :: symprec

     end subroutine spg_refine_cell

     !-----------------------------------------------------------------!

     subroutine spg_get_smallest_lattice(smallest_lattice, lattice, &
                                         symprec)

       implicit none

       double precision, dimension(3,3),        intent(out) :: smallest_lattice
       double precision, dimension(3,3),        intent(in)  :: lattice
       double precision,                        intent(in)  :: symprec

     end subroutine spg_get_smallest_lattice

     !-----------------------------------------------------------------!

     subroutine spg_get_multiplicity(size, lattice, position, types, &
                                     num_atom, symprec)

       implicit none

       integer,                                 intent(out) :: size  
       double precision, dimension(3,3),        intent(in)  :: lattice
       integer,                                 intent(in)  :: num_atom
       double precision, dimension(3,num_atom), intent(in)  :: position
       integer,          dimension(num_atom),   intent(in)  :: types
       double precision,                        intent(in)  :: symprec

     end subroutine spg_get_multiplicity

     !-----------------------------------------------------------------!

     subroutine spg_get_max_multiplicity(size, lattice, position, types, &
                                         num_atom, symprec)

       implicit none

       integer,                                 intent(out) :: size  
       double precision, dimension(3,3),        intent(in)  :: lattice
       integer,                                 intent(in)  :: num_atom
       double precision, dimension(3,num_atom), intent(in)  :: position
       integer,          dimension(num_atom),   intent(in)  :: types
       double precision,                        intent(in)  :: symprec

     end subroutine spg_get_max_multiplicity

     !-----------------------------------------------------------------!

     subroutine spg_get_ir_kpoints(num_ir_kpoints, map, kpoints,      &
                                   num_kpoints, lattice, position,    &
                                   types, num_atom, is_time_reversal, &
                                   symprec)

       implicit none

       integer,                                    intent(out) :: num_ir_kpoints
       integer,                                    intent(in)  :: num_kpoints
       integer,          dimension(num_kpoints),   intent(out) :: map
       double precision, dimension(3,num_kpoints), intent(in)  :: kpoints
       double precision, dimension(3,3),           intent(in)  :: lattice
       integer,                                    intent(in)  :: num_atom
       double precision, dimension(3,num_atom),    intent(in)  :: position
       integer,          dimension(num_atom),      intent(in)  :: types
       integer,                                    intent(in)  :: is_time_reversal
       double precision,                           intent(in)  :: symprec

     end subroutine spg_get_ir_kpoints

     !-----------------------------------------------------------------!

     subroutine spg_get_ir_reciprocal_mesh(num_ir_kpoints, grid_point, &
          map, mesh, is_shift, is_time_reversal, lattice, position,    &
          types, num_atom, symprec)

       implicit none

       integer,                                    intent(out) :: num_ir_kpoints
       integer,          dimension(3),             intent(in)  :: mesh
       integer,          dimension(3,mesh(1)*mesh(2)*mesh(3)), &
                                                   intent(out) :: grid_point
       integer,          dimension(mesh(1)*mesh(2)*mesh(3)),   &
                                                   intent(out) :: map
       integer,          dimension(3),             intent(in)  :: is_shift
       integer,                                    intent(in)  :: is_time_reversal
       double precision, dimension(3,3),           intent(in)  :: lattice
       integer,                                    intent(in)  :: num_atom
       double precision, dimension(3,num_atom),    intent(in)  :: position
       integer,          dimension(num_atom),      intent(in)  :: types
       double precision,                           intent(in)  :: symprec

     end subroutine spg_get_ir_reciprocal_mesh

     !-----------------------------------------------------------------!

     subroutine spg_get_stabilized_reciprocal_mesh(num_ir_kpoints,    &
          grid_point, map, mesh, is_shift, is_time_reversal, lattice, &
          num_rot, rotations, num_q, qpoints, symprec)

       implicit none

       integer,                                    intent(out) :: num_ir_kpoints
       integer,          dimension(:,:),           intent(out) :: grid_point
       integer,          dimension(:),             intent(out) :: map
       double precision, dimension(3),             intent(in)  :: mesh
       integer,          dimension(3),             intent(in)  :: is_shift
       integer,                                    intent(in)  :: is_time_reversal
       double precision, dimension(3,3),           intent(in)  :: lattice
       integer,                                    intent(in)  :: num_rot
       double precision, dimension(3,3,num_rot),   intent(in)  :: rotations
       integer,                                    intent(in)  :: num_q
       integer,          dimension(3,num_q),       intent(in)  :: qpoints
       double precision,                           intent(in)  :: symprec

     end subroutine spg_get_stabilized_reciprocal_mesh

     !-----------------------------------------------------------------!

     subroutine spg_get_triplets_reciprocal_mesh(num_ir_kpoints,      &
          triplets, weight_triplets, grid_point, num_triplets, mesh,  &
          is_time_reversal, lattice, num_rot, rotations, symprec)

       implicit none

       integer,                                     intent(out) :: num_ir_kpoints
       integer,                                     intent(in)  :: num_triplets
       integer,          dimension(3,num_triplets), intent(out) :: triplets
       integer,          dimension(num_triplets),   intent(out) :: weight_triplets
       integer,          dimension(:,:),            intent(out) :: grid_point
       double precision, dimension(3),              intent(in)  :: mesh
       integer,                                     intent(in)  :: is_time_reversal
       double precision, dimension(3,3),            intent(in)  :: lattice
       integer,                                     intent(in)  :: num_rot
       double precision, dimension(3,3,num_rot),    intent(in)  :: rotations
       double precision,                            intent(in)  :: symprec

     end subroutine spg_get_triplets_reciprocal_mesh

  end interface

end module spglib

