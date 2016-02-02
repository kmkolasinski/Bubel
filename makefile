all: quantulaba
debug:
		make C=ifortDEBUG
build:
		make C=ifort
zeus:
		make C=ifortZEUS
# --------------------------------------------------------------
# make C=ifort 		- to compile optimized code (default option)
# make C=ifortDEBUG - to compile DEUBG code
# make C=ifortZEUS  - my cutomized version for ZUES claster
# --------------------------------------------------------------


# Use commands to compile Quantulaba with different solvers,
# no option specified uses superLU library/
USED_LIBRARY= -DUSE_PARDISO
#-DUSE_PARDISO
#-DUSE_UMF_PACK


# -------------------------------------------------
SUPERLU_FILES=
UMFPACK_FILES=
# -------------------------------------------------


BUILD_DIR=
LIB_DIR=lib/

FOPT= -O3 -132
COPT= -03

ifeq ($(C),ifort)
else ifeq ($(C),ifortDEBUG)
FOPT= -O0 -132 -traceback -O0 -fstack-protector -check all  -assume realloc_lhs -ftrapuv -fpe0 -warn -traceback -debug extended
COPT= -00 -O0 -Wall -g
else ifeq ($(C),ifortZEUS)
endif

CC=gcc
CCFLAGS=-c $(COPT)
FF=ifort
FCFLAGS=-c  -Iinclude $(USED_LIBRARY) $(FOPT)
LDFLAGS=$(LIBS)  -mkl -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -L${MKLROOT}/lib/intel64



F90_SOURCES= modcommons.f90 \
		modinip.f90 \
		modunits.f90 \
		modalgs.f90 \
		modutils.f90 \
		modsys.f90  \
		modshape.f90  \
		modlead.f90  \
		modscatter.f90


C_SOURCES= \
		$(SUPERLU_FILES) \
		$(UMFPACK_FILES)

F90_OBJECTS=$(F90_SOURCES:%.f90=%.o)
BUILD_F90_OBJECTS=$(patsubst %.o,$(BUILD_DIR)%.o,$(F90_OBJECTS))



EXECUTABLE=quantulaba


quantulaba: main.f90 $(F90_OBJECTS)
		$(FF) main.f90  -I$(BUILD_DIR) -Iinclude $(BUILD_F90_OBJECTS) $(LDFLAGS) -o $@

.f90.o:
		$(FF) $(FCFLAGS) $< -o $@


modcommons.o: modcommons.f90
	$(FF) $(FCFLAGS) modcommons.f90 -o $(BUILD_DIR)$@

modutils.o: modutils.f90
	$(FF) $(FCFLAGS) modutils.f90 -o $(BUILD_DIR)$@

modinip.o: modinip.f90
	$(FF) $(FCFLAGS) modinip.f90 -o $(BUILD_DIR)$@

modsys.o: modsys.f90
	$(FF) $(FCFLAGS) modsys.f90 -o $(BUILD_DIR)$@

modlead.o: modlead.f90
	$(FF) $(FCFLAGS) modlead.f90 -o $(BUILD_DIR)$@

modscatter.o: modscatter.f90
	$(FF) $(FCFLAGS) modscatter.f90 -o $(BUILD_DIR)$@

modunits.o: modunits.f90
	$(FF) $(FCFLAGS) modunits.f90 -o $(BUILD_DIR)$@

modshape.o: modshape.f90
	$(FF) $(FCFLAGS) modshape.f90 -o $(BUILD_DIR)$@

modalgs.o: modalgs.f90
	$(FF) $(FCFLAGS)  modalgs.f90 -o $(BUILD_DIR)$@

zgssv.o:c_fortran_zgssv.c
	$(CC)   $(CCFLAGS) c_fortran_zgssv.c -o $(BUILD_DIR)$@

dgssv.o:c_fortran_dgssv.c
	$(CC)   $(CCFLAGS) c_fortran_dgssv.c -o $(BUILD_DIR)$@

umfpack.o:umfpack.f90
	$(FF) $(FCFLAGS) umfpack.f90 -o $(BUILD_DIR)$@

clean:
		rm -f $(BUILD_DIR)*.o $(BUILD_DIR)*.mod *.txt *.dat *.xml *.png *.pdf fort.* 2> /dev/null
		rm -f plots/*.txt plots/*.dat plots/*.xml plots/*.png plots/*.pdf plots/fort.* 2> /dev/null


dlib:
		ifort $(UMFPACK_MACRO) -shared -fpic -Iinclude  $(F90_SOURCES)  $(SUPERLU_FILES) -o $(LIB_DIR)libquantulaba.so

slib:
		ar rcs $(LIB_DIR)libquantulaba.a $(BUILD_DIR)*.o
