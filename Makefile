TARGET = spmv_main
OBJECTS = spmv_main.o mesh_gen.o testSpMV_mpi.o
FC = mpifort
#FC = mpiifort

FFLAGS = -lmetis -I/usr/local/include -L/usr/local/lib
LDFLAGS =

# for gfortran
ifeq (${FC},gfortran)
        #FFLAGS += -fimplicit-none -fbounds-check
        #LDFLAGS += -fopenmp -llapack -lblas
endif

# for ifort
ifeq (${FC},ifort)
        #MKLROOT = /opt/intel/composer_xe_2011_sp1.7.256/mkl
        #FFLAGS += -I${MKLROOT}/include/ia32 -I${MKLROOT}/include
        #LDFLAGS += -L${MKLROOT}/lib/ia32 ${MKLROOT}/lib/ia32/libmkl_blas95.a
        #LDFLAGS += ${MKLROOT}/lib/ia32/libmkl_lapack95.a
        #LDFLAGS += -lmkl_intel -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
endif

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	${FC} -c $< ${LDFLAGS}

${TARGET}:${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${FFLAGS}
clean:
	rm -f $(TARGET) $(OBJ)  *.o
