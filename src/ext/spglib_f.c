/* ------------------------------------------------------------------ */
/* This is a Fortran wrapper for some routines of the SPGLIB library. */
/* This file is not part of the SPGLIB package.  It assumes a single  */
/* underscore added by the fortran compiler to the function names.    */
/* ------------------------------------------------------------------ */
/* 2011-09-17 Alexander Urban (AU)                                    */
/* ------------------------------------------------------------------ */

#include "spglib.h"
#include <string.h>

void spg_get_symmetry_(int*    nsym, 
                       int     rot[][3][3], 
                       double  trans[][3], 
                       int*    size,
                       double  lattice[3][3], 
                       double  position[][3], 
                       int     types[],
                       int*    num_atom, 
                       double* symprec);

void spg_get_international_(int*    spacegroup, 
                            char    symbol[21], 
                            double  lattice[3][3],
                            double  position[][3], 
                            int     types[], 
                            int*    num_atom,
                            double* symprec);

void spg_get_schoenflies_(int*    spacegroup, 
                          char    symbol[10], 
                          double  lattice[3][3],
                          double  position[][3], 
                          int     types[], 
                          int*    num_atom,
                          double* symprec);

void spg_find_primitive_(double  lattice[3][3], 
                         double  position[][3],
                         int     types[], 
                         int*    num_atom, 
                         double* symprec);

void spg_refine_cell_(double  lattice[3][3],
                      double  position[][3],
                      int     types[],
                      int*    num_atom,
                      double* symprec);

void spg_get_smallest_lattice_(double  smallest_lattice[3][3],
                               double  lattice[3][3],
                               double* symprec);

void spg_get_multiplicity_(int*    size, 
                           double  lattice[3][3], 
                           double  position[][3],
                           int     types[], 
                           int*    num_atom, 
                           double* symprec);

void spg_get_max_multiplicity_(int*    size, 
                               double  lattice[3][3], 
                               double  position[][3],
                               int     types[], 
                               int*    num_atom, 
                               double* symprec);

void spg_get_ir_kpoints_(int*    num_ir_kpoints,
                         int     map[],
                         double  kpoints[][3],
                         int*    num_kpoint,
                         double  lattice[3][3],
                         double  position[][3],
                         int     types[],
                         int*    num_atom,
                         int*    is_time_reversal,
                         double* symprec);

void spg_get_ir_reciprocal_mesh_(int*    num_ir_kpoints,
                                 int     grid_point[][3],
                                 int     map[],
                                 int     mesh[3],
                                 int     is_shift[3],
                                 int*    is_time_reversal,
                                 double  lattice[3][3],
                                 double  position[][3],
                                 int     types[],
                                 int*    num_atom,
                                 double* symprec);

void spg_get_stabilized_reciprocal_mesh_(int*    num_ir_kpoints,
                                         int     grid_point[][3],
                                         int     map[],
                                         int     mesh[3],
                                         int     is_shift[3],
                                         int*    is_time_reversal,
                                         double  lattice[3][3],
                                         int*    num_rot,
                                         int     rotations[][3][3],
                                         int*    num_q,
                                         double  qpoints[][3],
                                         double* symprec);

// void spg_get_triplets_reciprocal_mesh_(int*    num_ir_kpoints,
//                                        int     triplets[][3],
//                                        int     weight_triplets[],
//                                        int     grid_point[][3],
//                                        int*    num_triplets,
//                                        int     mesh[3],
//                                        int*    is_time_reversal,
//                                        double  lattice[3][3],
//                                        int*    num_rot,
//                                        int     rotations[][3][3],
//                                        double* symprec);

/* ------------------------------------------------------------------ */
/*                      wrappers implementation                       */
/* ------------------------------------------------------------------ */

void spg_get_symmetry_(int*    nsym, 
                       int     rot[][3][3], 
                       double  trans[][3], 
                       int*    size,
                       double  lattice[3][3], 
                       double  position[][3], 
                       int     types[],
                       int*    num_atom, 
                       double* symprec) {
  
  *nsym = spg_get_symmetry(rot, trans, *size, lattice, position,
                           types, *num_atom, *symprec);

}

/* ------------------------------------------------------------------ */

void spg_get_international_(int*    spacegroup, 
                            char    symbol[21], 
                            double  lattice[3][3],
                            double  position[][3], 
                            int     types[], 
                            int*    num_atom,
                            double* symprec) {
  char symbol_c[21];
  int i, length;

  *spacegroup = spg_get_international(symbol_c, lattice, position, types,
                                      *num_atom, *symprec);
  length = strlen(symbol_c);
  strncpy(symbol, symbol_c, length);

  for (i=length; i<21; i++) { symbol[i] = ' '; }
}

/* ------------------------------------------------------------------ */

void spg_get_schoenflies_(int*    spacegroup, 
                          char    symbol[10], 
                          double  lattice[3][3],
                          double  position[][3], 
                          int     types[], 
                          int*    num_atom,
                          double* symprec) {
  char symbol_c[10];
  int i, length;
    
  *spacegroup = spg_get_schoenflies(symbol_c, lattice, position, types,
                                    *num_atom, *symprec);
  length = strlen(symbol_c);
  strncpy(symbol, symbol_c, length);
  
  for (i=length; i<10; i++) { symbol[i] = ' '; }
}

/* ------------------------------------------------------------------ */

void spg_find_primitive_(double  lattice[3][3], 
                         double  position[][3],
                         int     types[], 
                         int*    num_atom, 
                         double* symprec) {

  *num_atom = spg_find_primitive(lattice, position, types, *num_atom,
                                 *symprec);
}

/* ------------------------------------------------------------------ */

void spg_refine_cell_(double  lattice[3][3],
                      double  position[][3],
                      int     types[],
                      int*    num_atom,
                      double* symprec) {

  *num_atom = spg_refine_cell(lattice, position, types, *num_atom,
                              *symprec);
}

/* ------------------------------------------------------------------ */

void spg_get_smallest_lattice_(double  smallest_lattice[3][3],
                               double  lattice[3][3],
                               double* symprec) {

  spg_get_smallest_lattice(smallest_lattice, lattice, *symprec);
}

/* ------------------------------------------------------------------ */

void spg_get_multiplicity_(int*    size, 
                           double  lattice[3][3], 
                           double  position[][3],
                           int     types[], 
                           int*    num_atom, 
                           double* symprec) {
  
  *size = spg_get_multiplicity(lattice, position, types, *num_atom, 
                               *symprec);
}

/* ------------------------------------------------------------------ */

void spg_get_max_multiplicity_(int*    size, 
                               double  lattice[3][3], 
                               double  position[][3],
                               int     types[], 
                               int*    num_atom, 
                               double* symprec) {

  *size = spg_get_max_multiplicity(lattice, position, types, *num_atom, 
                                   *symprec);
}

/* ------------------------------------------------------------------ */

void spg_get_ir_kpoints_(int*    num_ir_kpoints,
                         int     map[],
                         double  kpoints[][3],
                         int*    num_kpoint,
                         double  lattice[3][3],
                         double  position[][3],
                         int     types[],
                         int*    num_atom,
                         int*    is_time_reversal,
                         double* symprec) {

  *num_ir_kpoints = spg_get_ir_kpoints(map, kpoints, *num_kpoint, 
                                       lattice, position, types, 
                                       *num_atom, *is_time_reversal,
                                       *symprec);
}

/* ------------------------------------------------------------------ */

void spg_get_ir_reciprocal_mesh_(int*    num_ir_kpoints,
                                 int     grid_point[][3],
                                 int     map[],
                                 int     mesh[3],
                                 int     is_shift[3],
                                 int*    is_time_reversal,
                                 double  lattice[3][3],
                                 double  position[][3],
                                 int     types[],
                                 int*    num_atom,
                                 double* symprec) {

  *num_ir_kpoints = spg_get_ir_reciprocal_mesh(grid_point, map, mesh,
                                               is_shift, *is_time_reversal,
                                               lattice, position, types,
                                               *num_atom, *symprec);
}

/* ------------------------------------------------------------------ */

void spg_get_stabilized_reciprocal_mesh_(int*    num_ir_kpoints,
                                         int     grid_point[][3],
                                         int     map[],
                                         int     mesh[3],
                                         int     is_shift[3],
                                         int*    is_time_reversal,
                                         double  lattice[3][3],
                                         int*    num_rot,
                                         int     rotations[][3][3],
                                         int*    num_q,
                                         double  qpoints[][3],
                                         double* symprec) {

  *num_ir_kpoints = spg_get_stabilized_reciprocal_mesh(grid_point, map,
                                                       mesh, is_shift,
                                                       *is_time_reversal, 
                                                       lattice, *num_rot,
                                                       rotations, *num_q,
                                                       qpoints, *symprec);
}

/* ------------------------------------------------------------------ */
// 
// void spg_get_triplets_reciprocal_mesh_(int*    num_ir_kpoints,
//                                        int     triplets[][3],
//                                        int     weight_triplets[],
//                                        int     grid_point[][3],
//                                        int*    num_triplets,
//                                        int     mesh[3],
//                                        int*    is_time_reversal,
//                                        double  lattice[3][3],
//                                        int*    num_rot,
//                                        int     rotations[][3][3],
//                                        double* symprec) {
// 
//   *num_ir_kpoints = spg_get_triplets_reciprocal_mesh(triplets, weight_triplets,
//                                                      grid_point, *num_triplets,
//                                                      mesh, *is_time_reversal,
//                                                      lattice, *num_rot,
//                                                      rotations, *symprec);
// }
