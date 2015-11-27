all: build_simple_sqlat \
	 build_simple_graphene \
	 build_simple_graphene2 \
	 build_rashba_2deg \
	 build_carbon_nanotube \
	 build_strange_graphene \
	 build_superconductor \
	 build_supercnd_trans \
	 build_simple_galerkin1d

run:
		 ./simple_sqlat ; \
		 ./simple_graphene ; \
		 ./simple_graphene2 ; \
		 ./rashba_2deg ; \
		 ./carbon_nanotube ; \
		 ./strange_graphene ;\
		 ./superconductor; \
		 ./simple_galerkin1d


COMPILATION_TYPE=ifort
# To build simple example of effective band mass on quare lattice
build_simple_sqlat:
	make simple_sqlat C=$(COMPILATION_TYPE) P=simple_sqlat

build_simple_graphene:
	make simple_graphene C=$(COMPILATION_TYPE) P=simple_graphene

build_simple_graphene2:
	make simple_graphene2 C=$(COMPILATION_TYPE) P=simple_graphene2

build_rashba_2deg:
	make rashba_2deg C=$(COMPILATION_TYPE) P=rashba_2deg

build_carbon_nanotube:
	make carbon_nanotube C=$(COMPILATION_TYPE) P=carbon_nanotube

build_strange_graphene:
	make strange_graphene C=$(COMPILATION_TYPE) P=strange_graphene

build_superconductor:
	make superconductor C=$(COMPILATION_TYPE) P=superconductor

build_supercnd_trans:
	make supercnd_trans C=$(COMPILATION_TYPE) P=supercnd_trans

build_simple_galerkin1d:
	make simple_galerkin1d C=$(COMPILATION_TYPE) P=simple_galerkin1d

EXEC_FILE=$(P)
# -----------------------------------
# Use command:
# UMFPACK_MACRO=-DUSE_UMF_PACK
# to compile with UMFPACKIEM other macros
# leave it empty to compile it with SuperLU
#-DUSE_UMF_PACK
#-DUSE_PARDISO
UMFPACK_MACRO=-DUSE_PARDISO
#UMFPACK_MACRO=



ifeq ($(C),ifort)
FC=ifort

BASEDIR=/home/mkk/libs
FBFLAGS=  -O3  -132 -I..
#-I$(BASEDIR)/XC -I..

ifeq ($(UMFPACK_MACRO),-DUSE_UMF_PACK)
LIBS= $(BASEDIR)/libumfpack.a $(BASEDIR)/libamd.a
FCFLAGS= -c -O3  -132  $(UMFPACK_MACRO)
#-I$(BASEDIR)/XC
FCCFLAGS= -c -O3
SUPERLU_FILES=
UMFPACK_FILES=umfpack.o
else ifeq ($(UMFPACK_MACRO),-DUSE_PARDISO)
LIBS=
FCFLAGS= -c -O3  -132  $(UMFPACK_MACRO)
#-I$(BASEDIR)/XC
FCCFLAGS= -c -O3
SUPERLU_FILES=
UMFPACK_FILES=
else
LIBS= $(BASEDIR)/libsuperlu_4.3.a
FCFLAGS= -c -O3  -132  -I$(BASEDIR)/SuperLU_4.3/SRC $(UMFPACK_MACRO)
#-I$(BASEDIR)/XC
FCCFLAGS= -c -O3 -I$(BASEDIR)/SuperLU_4.3/SRC
SUPERLU_FILES=zgssv.o
UMFPACK_FILES=
endif
FLIBS=   $(LIBS)  -mkl -lmkl_lapack95_lp64 -lmkl_intel_lp64  -L${MKLROOT}/lib/intel64
#-static-intel $(BASEDIR)/libxc.a

else ifeq ($(C),ifortDEBUG)
FC=ifort
BASEDIR =/home/mkk/libs
FBFLAGS =  -O0 -132 -I..

ifeq ($(UMFPACK_MACRO),-DUSE_UMF_PACK)
FCFLAGS = -c -132 -traceback -O0 -check all -fpe0 -warn -traceback -debug extended  $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c -O0 -Wall -g
LIBS= $(BASEDIR)/libumfpack.a $(BASEDIR)/libamd.a
SUPERLU_FILES=
UMFPACK_FILES=umfpack.o
else ifeq ($(UMFPACK_MACRO),-DUSE_PARDISO)
LIBS=
FCFLAGS = -c -132 -traceback -O0 -check all -fpe0 -warn -traceback -debug extended  $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c -O0 -Wall -g
SUPERLU_FILES=
UMFPACK_FILES=
else
LIBS= $(BASEDIR)/libsuperlu_4.3.a
FCFLAGS = -c -132 -traceback -O0 -check all -fpe0 -warn -traceback -debug extended -I$(BASEDIR)/SuperLU_4.3/SRC $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c -O0 -Wall -g -I$(BASEDIR)/SuperLU_4.3/SRC
SUPERLU_FILES=zgssv.o
UMFPACK_FILES=
endif
FLIBS=   $(LIBS)  -mkl
# -static-intel  $(BASEDIR)/libxc.a

endif



$(EXEC_FILE): $(EXEC_FILE).f90  $(UMFPACK_FILES)
	$(FC) $(FBFLAGS)  $(EXEC_FILE).f90 ../libquantulaba.a $(FLIBS)   -o $@
#	$(FC) $(FBFLAGS)  $(EXEC_FILE).f90 -L..  -lquantulaba $(FLIBS)   -o $@


clean:
	rm -f */*.dat */*.pdf */*.png */*.xml */*.txt  2> /dev/null
	rm -f simple_sqlat \
		  simple_graphene \
		  simple_graphene2 \
		  carbon_nanotube \
		  rashba_2deg \
		  strange_graphene \
		  superconductor \
		  supercnd_trans \
		  simple_galerkin1d 2> /dev/null
