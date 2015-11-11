all: quantulaba

build:
	make C=ifort

debug:
	make C=ifortDEBUG

zeus:
	make C=ifortZEUS


#	ar rcs libquantulaba.a *.o
# -----------------------------------
# Uzyj polecenia
# UMFPACK_MACRO=-DUSE_UMF_PACK
# aby skompilowac z UMFPACKIEM
UMFPACK_MACRO=-DUSE_PARDISO
#-DUSE_UMF_PACK
#-DUSE_PARDISO

#UMFPACK_MACRO=



lib:
	ifort $(UMFPACK_MACRO) -shared -fpic  modcommons.f90 \
									      modunits.f90 \
										  modutils.f90 \
										  modsys.f90 \
										  modshape.f90 \
										  modlead.f90 \
										  modscatter.f90 \
										  $(SUPERLU_FILES) -o libquantulaba.so

slib:
	ar rcs libquantulaba.a *.o

OPTS= -O3

ifeq ($(C),ifort)
FC=ifort

BASEDIR=/home/mkk/libs
FBFLAGS=  $(OPTS)  -132
#-I$(BASEDIR)/XC


ifeq ($(UMFPACK_MACRO),-DUSE_UMF_PACK)
LIBS= $(BASEDIR)/libumfpack.a $(BASEDIR)/libamd.a
FCFLAGS= -c $(OPTS)  -132  $(UMFPACK_MACRO)
#-I$(BASEDIR)/XC
FCCFLAGS= -c $(OPTS)
SUPERLU_FILES=
UMFPACK_FILES=umfpack.o
else ifeq ($(UMFPACK_MACRO),-DUSE_PARDISO)
LIBS=
FCFLAGS= -c $(OPTS)  -132  $(UMFPACK_MACRO)
#-I$(BASEDIR)/XC
FCCFLAGS= -c $(OPTS)
SUPERLU_FILES=
UMFPACK_FILES=
else
LIBS= $(BASEDIR)/libsuperlu_4.3.a
FCFLAGS= -c  $(OPTS)  -132  -I$(BASEDIR)/SuperLU_4.3/SRC $(UMFPACK_MACRO)
#-I$(BASEDIR)/XC
FCCFLAGS= -c $(OPTS) -I$(BASEDIR)/SuperLU_4.3/SRC
SUPERLU_FILES=zgssv.o
UMFPACK_FILES=
endif
FLIBS=   $(LIBS) -mkl -lmkl_lapack95_lp64
#FLIBS=   $(LIBS)   -mkl -lmkl_lapack95_lp64 -lmkl_intel_lp64  -L${MKLROOT}/lib/intel64
#-static-intel $(BASEDIR)/libxc.a

else ifeq ($(C),ifortDEBUG)
FC=ifort
BASEDIR =/home/mkk/libs
FBFLAGS =  -O0 -132

ifeq ($(UMFPACK_MACRO),-DUSE_UMF_PACK)
FCFLAGS = -c -132 -traceback -O0 -fstack-protector -check all  -assume realloc_lhs -ftrapuv -fpe0 -warn -traceback -debug extended  $(UMFPACK_MACRO)
# -I$(BASEDIR)/XC
FCCFLAGS= -c -O0 -Wall -g
LIBS= $(BASEDIR)/libumfpack.a $(BASEDIR)/libamd.a
SUPERLU_FILES=
UMFPACK_FILES=umfpack.o
else ifeq ($(UMFPACK_MACRO),-DUSE_PARDISO)
LIBS=
FCFLAGS = -c -132 -traceback -O0 -fstack-protector -check all  -assume realloc_lhs  -ftrapuv -fpe0 -warn -traceback -debug extended  $(UMFPACK_MACRO)
#-I$(BASEDIR)/XC
FCCFLAGS= -c -O0 -Wall -g
SUPERLU_FILES=
UMFPACK_FILES=
else
LIBS= $(BASEDIR)/libsuperlu_4.3.a
FCFLAGS = -c -132 -traceback -O0 -fstack-protector -check all  -assume realloc_lhs -ftrapuv -fpe0 -warn -traceback -debug extended -I$(BASEDIR)/SuperLU_4.3/SRC $(UMFPACK_MACRO)
# -I$(BASEDIR)/XC
FCCFLAGS= -c -O0 -Wall -g -I$(BASEDIR)/SuperLU_4.3/SRC
SUPERLU_FILES=zgssv.o
UMFPACK_FILES=
endif
FLIBS=   $(LIBS)  -mkl -lmkl_lapack95_lp64 -lmkl_intel_lp64  -L${MKLROOT}/lib/intel64
# -static-intel  $(BASEDIR)/libxc.a

else ifeq ($(C),ifortZEUS)
FC=ifort
BASEDIR=/people/gjkolasi
FBFLAGS =  -O0 -132

ifeq ($(UMFPACK_MACRO),-DUSE_UMF_PACK)
LIBS= $(BASEDIR)/libumfpack.a $(BASEDIR)/libamd.a
FCFLAGS= -c $(OPTS)  -132  $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c $(OPTS)
SUPERLU_FILES=
UMFPACK_FILES=umfpack.o
else ifeq ($(UMFPACK_MACRO),-DUSE_PARDISO)
LIBS=
FCFLAGS= -c $(OPTS)  -132  $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c $(OPTS)
SUPERLU_FILES=
UMFPACK_FILES=
else
LIBS= $(BASEDIR)/libsuperlu_4.3.a
FCFLAGS= -c $(OPTS)  -132  -I$(BASEDIR)/SuperLU_4.3/SRC $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c $(OPTS) -I$(BASEDIR)/SuperLU_4.3/SRC -I$(BASEDIR)/XC
SUPERLU_FILES=zgssv.o
UMFPACK_FILES=
endif
FLIBS=   $(LIBS)  -mkl $(BASEDIR)/libxc.a

endif



quantulaba: main.f90 modcommons.o modinip.o modunits.o modutils.o $(UMFPACK_FILES) modsys.o modshape.o modlead.o modscatter.o $(SUPERLU_FILES)
	$(FC) $(FBFLAGS)  main.f90 *.o $(FLIBS)   -o $@

modcommons.o: modcommons.f90
	$(FC) $(FCFLAGS) modcommons.f90 -o $@

modutils.o: modutils.f90
	$(FC) $(FCFLAGS) modutils.f90 -o $@

modinip.o: modinip.f90
	$(FC) $(FCFLAGS) modinip.f90 -o $@

modsys.o: modsys.f90
	$(FC) $(FCFLAGS) modsys.f90 -o $@

modlead.o: modlead.f90
	$(FC) $(FCFLAGS) modlead.f90 -o $@

modscatter.o: modscatter.f90
	$(FC) $(FCFLAGS) modscatter.f90 -o $@

modunits.o: modunits.f90
	$(FC) $(FCFLAGS) modunits.f90 -o $@

modshape.o: modshape.f90
	$(FC) $(FCFLAGS) modshape.f90 -o $@

zgssv.o:c_fortran_zgssv.c
	gcc   $(FCCFLAGS) c_fortran_zgssv.c -o $@

umfpack.o:umfpack.f90
	$(FC) $(FCFLAGS) umfpack.f90 -o $@




clean:
	rm -f *.o *.mod *.txt *.dat *.xml *.png *.pdf fort.* 2> /dev/null
