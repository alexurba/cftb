LIB       = ../lib/
INC       = ../inc/

FC        = gfortran -c
# FC        = ifort -c

AR        = ar
FCFLAGS   = -O2
ARFLAGS   = -crusv

TRGLAPACK = $(LIB)liblapack.a
include ./objlapack.inc

TRGBLAS   = $(LIB)libblas.a
include ./objblas.inc

TRGSYMSPG = $(LIB)libsymspg.a

#----------------------------------------------------------------------#
#                                Rules                                 #
#----------------------------------------------------------------------#

.SUFFIXES: .F90 .o .mod .f $(SUFFIXES)
.PHONY : clean all lapack blas liblapack libblas

all    : $(TRGLAPACK) $(TRGBLAS)
lapack : $(TRGLAPACK)
blas   : $(TRGBLAS)
symspg : $(TRGSYMSPG)
liblapack : $(TRGLAPACK)
libblas   : $(TRGBLAS)
libsymspg : $(TRGSYMSPG)

%.o %.mod    : %.F90 Makefile
	$(FC) $(FCFLAGS) $< -o $*.o

%.o %.mod    : %.f90 Makefile
	$(FC) $(FCFLAGS) $< -o $*.o

%.o %.mod    : %.f Makefile
	$(FC) $(FCFLAGS) $< -o $*.o

$(TRGLAPACK) : $(OBJLAPACK) Makefile
	$(AR) $(ARFLAGS) $@ $(OBJLAPACK)

$(TRGBLAS)   : $(OBJBLAS) Makefile
	$(AR) $(ARFLAGS) $@ $(OBJBLAS)

$(TRGSYMSPG) $(INC)spglib.h: Makefile
	cd ./spglib && ./configure
	cd ./spglib && $(MAKE)
	cp -a ./spglib/src/.libs/libsymspg.a $(TRGSYMSPG)
	cp -a ./spglib/src/spglib.h $(INC)

clean :
	if [ -e $(TRGLAPACK) ]; then rm -f $(TRGLAPACK); fi
	for f in $(OBJLAPACK); do if [ -e $$f ]; then rm -f $$f; fi; done
	if [ -e $(TRGBLAS) ]; then rm -f $(TRGBLAS); fi
	for f in $(OBJBLAS); do if [ -e $$f ]; then rm -f $$f; fi; done
	if [ -e $(TRGSYMSPG) ]; then rm -f $(TRGSYMSPG); fi
	if [ -e $(INC)spglib.h ]; then rm -f $(INC)spglib.h; fi
	cd ./spglib && $(MAKE) clean
	for f in *~; do rm -f $$f; done
