all: simple_sqlat

build:
	make C=ifort

debug:
	make C=ifortDEBUG

# To build simple example of effective band mass on quare lattice
build_simple_sqlat:
	make simple_sqlat C=ifort P=simple_sqlat

build_simple_graphene:
	make simple_graphene C=ifort P=simple_graphene

build_simple_graphene2:
	make simple_graphene2 C=ifort P=simple_graphene2

build_rashba_2deg:
	make rashba_2deg C=ifort P=rashba_2deg

build_carbon_nanotube:
	make carbon_nanotube C=ifort P=carbon_nanotube

EXEC_FILE=$(P)
# -----------------------------------
# Uzyj polecenia
# UMFPACK_MACRO=-DUSE_UMF_PACK
# aby skompilowac z UMFPACKIEM
UMFPACK_MACRO=-DUSE_PARDISO
#-DUSE_UMF_PACK
#-DUSE_PARDISO

#UMFPACK_MACRO=
ifeq ($(C),ifort)
FC=ifort

BASEDIR=/home/mkk/libs
FBFLAGS=  -O3  -132 -I$(BASEDIR)/XC -I..


ifeq ($(UMFPACK_MACRO),-DUSE_UMF_PACK)
LIBS= $(BASEDIR)/libumfpack.a $(BASEDIR)/libamd.a
FCFLAGS= -c -O3  -132  $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c -O3
SUPERLU_FILES=
UMFPACK_FILES=umfpack.o
else ifeq ($(UMFPACK_MACRO),-DUSE_PARDISO)
LIBS=
FCFLAGS= -c -O3  -132  $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c -O3
SUPERLU_FILES=
UMFPACK_FILES=
else
LIBS= $(BASEDIR)/libsuperlu_4.3.a
FCFLAGS= -c -O3  -132  -I$(BASEDIR)/SuperLU_4.3/SRC $(UMFPACK_MACRO) -I$(BASEDIR)/XC
FCCFLAGS= -c -O3 -I$(BASEDIR)/SuperLU_4.3/SRC
SUPERLU_FILES=zgssv.o
UMFPACK_FILES=
endif
FLIBS=   $(LIBS)  -mkl -static-intel $(BASEDIR)/libxc.a

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
FLIBS=   $(LIBS)  -mkl -static-intel  $(BASEDIR)/libxc.a

endif



$(EXEC_FILE): $(EXEC_FILE).f90  $(UMFPACK_FILES)
	$(FC) $(FBFLAGS)  $(EXEC_FILE).f90 ../libquantulaba.a $(FLIBS)   -o $@


clean:
	rm -f */*.dat */*.pdf */*.png */*.xml */*.txt  2> /dev/null
	rm -f simple_sqlat simple_graphene simple_graphene2 carbon_nanotube rashba_2deg  2> /dev/null
