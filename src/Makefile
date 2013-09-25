
#---------------------------- GNU Fortran -----------------------------#

#FC       = gfortran -c
#CC       = gcc -c
#LD       = gfortran
#
#CHECKS   =  -Wall -pedantic -fbounds-check
#FCFLAGS  = -O2 $(CHECKS)
#CCFLAGS  = $(INCDIR)
#LDFLAGS  = -static-libgfortran $(CHECKS)

#--------------------------- Intel Fortran ----------------------------#

FC       = ifort -c
CC       = gcc -c
LD       = ifort
#
#CHECKS   = # -check all -check bounds
OPTFLAGS =
FCFLAGS  = -O2 $(OPTFLAGS) $(CHECKS)
CCFLAGS  = $(INCDIR)
LDFLAGS  = $(OPTFLAGS) -lmkl_intel_lp64  -lmkl_sequential -lmkl_core $(CHECKS)
#LDFLAGS  = $(OPTFLAGS) -i-static -lmkl -llapack $(CHECKS)

#----------------------------------------------------------------------#

BINDIR   = ../bin/
LIBDIR   = ./lib/
INCDIR   = -I./inc/
EXTDIR   = ./ext/
LIBSRC   = ./libsrc/

LIB      = arglib.o             \
           bfgs.o               \
           funclib.o            \
           harmonics.o          \
           interpol.o           \
           io.o                 \
           mbppio.o             \
           simplex.o            \
           sortlib.o            \
           timing.o             \
           voronoi.o            \
           xsflib.o             \
           spglib.o             \
           spglib_f.o           \
           $(LIBDIR)libsymspg.a
#           $(LIBDIR)liblapack.a \
#           $(LIBDIR)libblas.a   \


OBJECTS  = bondint.o   \
           constants.o \
           energy.o    \
           geometry.o  \
           pbc.o       \
           potential.o \
           symmetry.o  \
           tbinter.o   \
           tbio.o      \
           tbmatrix.o  \
           tbparam.o

TARGET1  = $(BINDIR)tb.x
TARGET2  = $(BINDIR)onsite.x
TARGET3  = $(BINDIR)bands.x
TARGET4  = $(BINDIR)fit.x
TARGET5  = $(BINDIR)potfit.x
TARGET6  = $(BINDIR)paraplot.x
TARGET7  = $(BINDIR)dos.x
TARGET8  = $(BINDIR)paramod.x

#----------------------------------------------------------------------#

.SUFFIXES   : .f90 .F90 .c .h .o .mod .x $(SUFFIXES)
.PHONY      : all clean veryclean tb onsite bands fit potfit paraplot dos

all         : $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) \
              $(TARGET6) $(TARGET7) $(TARGET8)

tb          : $(TARGET1)
onsite      : $(TARGET2)
bands       : $(TARGET3)
fit         : $(TARGET4)
potfit      : $(TARGET5)
paraplot    : $(TARGET6)
dos         : $(TARGET7)
paramod     : $(TARGET8)

clean       :
	rm *.o *.mod &1>&2 > /dev/null; exit 0
veryclean   :
	rm *.o *.mod *.x &1>&2 > /dev/null; exit 0

%.o %.mod : %.F90
	$(FC) $(FCFLAGS) $< -o $*.o

%.o %.mod : %.f90
	$(FC) $(FCFLAGS) $< -o $*.o

%.o : %.f90
	$(FC) $(FCFLAGS) $< -o $*.o

%.o : %.c
	$(CC) $(CCFLAGS) $< -o $*.o

%.o : $(EXTDIR)%.c Makefile
	$(CC) $(CCFLAGS) $< -o $*.o

%.o %.mod : $(EXTDIR)%.f90
	$(FC) $(FCFLAGS) $< -o $*.o

%.o %.mod : $(LIBDIR)%.f90
	$(FC) $(FCFLAGS) $< -o $*.o

$(BINDIR)%.x : %.o $(OBJECTS) $(LIB)
	$(LD) $(LDFLAGS) $< $(OBJECTS) $(LIB) -o $@

$(LIBDIR)%.a : $(LIBSRC)Makefile
	cd $(LIBSRC) && $(MAKE) $*

#--------------#
# dependencies #
#--------------#

tb.o       : tb.f90         $(OBJECTS)  $(LIB)
onsite.o   : onsite.f90     $(OBJECTS)  $(LIB)
bands.o    : bands.f90      $(OBJECTS)  $(LIB)
fit.o      : fit.f90        $(OBJECTS)  $(LIB)
potfit.o   : potfit.f90     $(OBJECTS)  $(LIB)
paramod.o  : paramod.f90    $(OBJECTS)  $(LIB)
paraplot.o : paraplot.f90   $(OBJECTS)  $(LIB)
dos.o      : dos.f90        $(OBJECTS)  $(LIB)

# libraries/modules belonging to this project:
bondint.o    : bondint.f90     interpol.o  io.o
constants.o  : constants.f90
energy.o     : energy.f90      constants.o  funclib.o  sortlib.o
geometry.o   : geometry.f90    constants.o  io.o  mbppio.o  xsflib.o
pbc.o        : pbc.f90         sortlib.o
potential.o  : potential.f90   bondint.o
symmetry.o   : symmetry.f90    spglib.o
tbio.o       : tbio.f90        arglib.o   constants.o  io.o
tbmatrix.o   : tbmatrix.f90    bondint.o  harmonics.o  io.o
tbparam.o    : tbparam.f90     bondint.o  io.o
tbinter.o    : tbinter.f90     constants.o energy.o geometry.o \
               io.o pbc.o potential.o tbio.o tbmatrix.o tbparam.o timing.o

# external libraries (not maintained within this repository):
arglib.o     : $(EXTDIR)arglib.f90
bfgs.o       : $(EXTDIR)bfgs.f90
funclib.o    : $(EXTDIR)funclib.f90
harmonics.o  : $(EXTDIR)harmonics.f90  funclib.o
interpol.o   : $(EXTDIR)interpol.f90
io.o         : $(EXTDIR)io.f90
mbppio.o     : $(EXTDIR)mbppio.f90     io.o
simplex.o    : $(EXTDIR)simplex.f90
sortlib.o    : $(EXTDIR)sortlib.f90
spglib.o     : $(EXTDIR)spglib.f90     spglib_f.o
spglib_f.o   : $(EXTDIR)spglib_f.c
timing.o     : $(EXTDIR)timing.f90
voronoi.o    : $(EXTDIR)voronoi.f90    sortlib.o
xsflib.o     : $(EXTDIR)xsflib.f90     io.o
