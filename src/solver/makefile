#Commons variables
HYPRE_DIR     =
METIS         =    /home/guy/soft/metis-4.0
MPI           =    /home/guy/soft/mpich2-install
FREEFEM       =    /home/guy/soft/FreeFem_dev/freefem++
MPI_INCLUDE     =  -I$(MPI)/include
METIS_INCLUDE   =  -I$(METIS)/Lib
METIS_LIB       =  $(METIS)/libmetis.a
FREEFEM_INCLUDE =  -I$(FREEFEM)/examples++-load/include
BLAS            =  /home/guy/Documents/soft/BLAS/blas_LINUX.a
# Variables needs by HIPS
HIPS_DIR      = /home/guy/soft/hips-1.2beta-rc3
HIPS_INCLUDE  = -I$(HIPS_DIR)/SRC/INCLUDE
LIB_DIR	      = -L$(HIPS_DIR)/LIB

LIBIO             = $(HIPS_DIR)/LIB/libio.a
LIBHIPSSEQUENTIAL = $(HIPS_DIR)/LIB/libhipssequential.a
LIBHIPS           = $(HIPS_DIR)/LIB/libhips.a
LIBIOMM	          = $(HIPS_DIR)/LIB/libio.a
LIBSPKIT	  = $(HIPS_DIR)/LIB/libspkit.a
#Variables needs by pARMS

PARMS_DIR     = /home/guy/soft/pARMS_2.2
PARMS_INCLUDE = -I$(PARMS_DIR)/INCLUDE
PARMS_LIB     =  -L$(PARMS_DIR)/LIB/ -lparms-O2
#Variables needs by HYPRE

HYPRE_DIR     = /home/guy/soft/hypre-2.4.0b
HYPRE_INCLUDE = -I$(HYPRE_DIR)/src/hypre/include/
HYPRE_LIB     = -L$(HYPRE_DIR)/src/lib/ -lHYPRE

#
# Compiler
#
ARCH       = -DLINUX
CC	   =  gcc       # C compiler 
MPICC      =  mpicc
FC         =  gfortran  # Fortran compiler 
MPIFC      =  mpif90
LD	   =  $(FC)     # Linker
MPILD      =  $(MPIFC)
CFLAGS	   =   -Wall -fPIC -DLINUX -DSPKIT_F  -DMPICH_IGNORE_CXX_SEEK  -DNDEBUG -rdynamic        # Additional C compiler flags
FFLAGS	   =   -shared -fPIC  -DNDEBUG -O3 -rdynamic        # Additional Fortran compiler flags
LFLAGS     =           # Additional linker flags
COPTFLAGS  =    -O3   # Optimization flags
FOPTFLAGS  =         # 


#Creation of Hips dynamic librairy
hips    : 
	$(MPICC) -c $(CFLAGS) $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(HIPS_INCLUDE) $(METIS_INCLUDE) hips_FreeFem.cpp 
	$(MPIFC)$  $(FFLAGS) hips_FreeFem.o -o hips_FreeFem.so $(LIBHIPS) $(LIBIO) $(LIBSPKIT) $(BLAS) $(METIS_LIB)

parms	:
	$(MPICC) -c $(CFLAGS) $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(PARMS_INCLUDE) $(METIS_INCLUDE) parms_FreeFem.cpp 
	$(MPIFC)$  $(FFLAGS) parms_FreeFem.o -o parms_FreeFem.so $(PARMS_LIB) $(BLAS) $(METIS_LIB)

hypre	:
	$(MPICC) -c $(CFLAGS) $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(HYPRE_INCLUDE) $(METIS_INCLUDE) hypre_FreeFem.cpp 
	$(MPIFC)$  $(FFLAGS) hypre_FreeFem.o -o hypre_FreeFem.so $(HYPRE_LIB) $(BLAS) $(METIS_LIB)

