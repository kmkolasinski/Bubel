all: build_simple_sqlat \
	 build_simple_graphene \
	 build_simple_graphene2 \
	 build_rashba_2deg \
	 build_carbon_nanotube \
	 build_strange_graphene \
	 build_superconductor \
	 build_supercnd_trans \
	 build_simple_galerkin1d \
	 build_recursive_wfm \
	 build_poisson2d \
	 build_benchmark

run:
		 ./simple_sqlat ; \
		 ./simple_graphene ; \
		 ./simple_graphene2 ; \
		 ./rashba_2deg ; \
		 ./carbon_nanotube ; \
		 ./strange_graphene ;\
		 ./superconductor; \
		 ./simple_galerkin1d; \
		 ./supercnd_trans; \
		 ./recursive_wfm; \
		 ./poisson2d;


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

build_recursive_wfm:
	make recursive_wfm C=$(COMPILATION_TYPE) P=recursive_wfm

build_poisson2d:
	make poisson2d C=$(COMPILATION_TYPE) P=poisson2d

build_benchmark:
	make benchmark C=$(COMPILATION_TYPE) P=benchmark


EXEC_FILE=$(P)


SUPERLU_FILES=
UMFPACK_FILES=

FOPT= -O3 -132
COPT= -03

ifeq ($(C),ifort)
else ifeq ($(C),ifortDEBUG)
FOPT= -O0 -132 -traceback -O0 -fstack-protector -check all  -assume realloc_lhs -ftrapuv -fpe0 -warn -traceback -debug extended
COPT= -00 -O0 -Wall -g
else ifeq ($(C),ifortZEUS)
endif

FF=ifort
FCFLAGS=-I../include -I..  $(FOPT)
LDFLAGS=$(LIBS)  -mkl -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -L${MKLROOT}/lib/intel64



$(EXEC_FILE): $(EXEC_FILE).f90  $(UMFPACK_FILES)
	$(FF) $(FCFLAGS)  $(EXEC_FILE).f90 ../lib/libquantulaba.a $(LDFLAGS) -o $@



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
		  simple_galerkin1d \
		  poisson2d \
		  recursive_wfm \
		  benchmark 2> /dev/null
