#!/bin/bash

# LAPACK
find lapack -name *.f \
  | sed 's/\.f/.o/g' \
  | awk '
{
  if (NR>1) {
    printf("\\\n            ")
  } else {
    printf("OBJLAPACK = ")
  }; 
  printf("%-25s", $1)
} 

END {
  print ""
}' > objlapack.inc


# BLAS
find blas -name *.f \
  | sed 's/\.f/.o/g' \
  | awk '
{
  if (NR>1) {
    printf("\\\n            ")
  } else {
    printf("OBJBLAS   = ")
  }; 
  printf("%-25s", $1)
} 

END {
  print ""
}' > objblas.inc

touch Makefile
