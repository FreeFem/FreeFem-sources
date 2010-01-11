#include makefile-gene.inc
#include makefile-common.inc
##include makefile-hips.inc
include makefile-parms.inc
include makefile-hypre.inc
#include makefile-slu_dist.inc
#include makefile-mumps.inc
#include makefile-pastix.inc
include makefile-common.inc
include makefile-hips.inc

#include makefile-sparsesolver.inc


#Creation of Hips dynamic librairy
hips    : 
	$(MPIFC) -c $(CFLAGS) $(CFLAGS_HIPS) $(FREEFEM_INCLUDE) $(HIPS_INCLUDE) $(METIS_INCLUDE)  hips_FreeFem.cpp skitf.f
	$(MPIFC)$  $(FFLAGS) skitf.o hips_FreeFem.o -o hips_FreeFem.so $(HIPS_LIB)  $(BLAS_LIB) $(METIS_LIB)

parms	:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_PARMS)  $(FREEFEM_INCLUDE) $(PARMS_INCLUDE) $(METIS_INCLUDE)  parms_FreeFem.cpp 
	$(MPIFC)$  $(FFLAGS) parms_FreeFem.o -o parms_FreeFem.so $(PARMS_LIB) $(BLAS_LIB) $(METIS_LIB) $(PARMS_OTHER_LIB)

hypre	:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_HYPRE) $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(HYPRE_INCLUDE) $(METIS_INCLUDE) $(HYPRE_OTHER_INCLUDE) hypre_FreeFem.cpp 
	$(MPIFC)$  $(FFLAGS)  hypre_FreeFem.o -o hypre_FreeFem.so $(HYPRE_LIB) $(BLAS_LIB) $(METIS_LIB) $(HYPRE_OTHER_LIB)

mumps	:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_MUMPS) $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(MUMPS_OTHER_INCLUDE) $(MUMPS_INCLUDE) $(SCOTCH_INCLUDE)  MUMPS_FreeFem.cpp
	$(MPIFC)  $(FFLAGS) $(MPI_INCLUDE) $(MUMPS_OTHER_INCLUDE) $(MUMPS_INCLUDE) MUMPS_FreeFem.o -o MUMPS_FreeFem.$(SUF) $(MPI_LIB) $(BLAS_LIB) $(BLACS_LIB) $(SCALAPACK_LIB) $(METIS_LIB) $(PORD_LIB) $(SCOTCH_LIB) $(MUMPS_LIB) $(MUMPS_OTHER_LIB)
#
rpastix_dist:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_PASTIX)  $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(PTSCOTCH_INCLUDE) $(REAL_PASTIX_INCLUDE) $(PASTIX_OTHER_INCLUDE) real_dpastix_FreeFem.cpp
	$(MPIFC) $(FFLAGS) $(PTSCOTCH_INCLUDE) $(REAL_PASTIX_INCLUDE) $(PASTIX_OTHER_INCLUDE) real_dpastix_FreeFem.o -o real_dpastix_FreeFem.$(SUF) $(MPI_LIB) $(BLAS_LIB) $(PTSCOTCH_LIB) $(METIS_LIB) $(REAL_PASTIX_LIB) $(PASTIX_OTHER_LIB)

rpastix:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_PASTIX) $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(SCOTCH_INCLUDE) $(REAL_PASTIX_INCLUDE) $(PASTIX_OTHER_INCLUDE) real_pastix_FreeFem.cpp
	$(MPIFC) $(FFLAGS) $(SCOTCH_INCLUDE) $(REAL_PASTIX_INCLUDE) $(PASTIX_OTHER_INCLUDE) real_pastix_FreeFem.o -o real_pastix_FreeFem.$(SUF) $(MPI_LIB) $(BLAS_LIB) $(SCOTCH_LIB) $(METIS_LIB) $(REAL_PASTIX_LIB) $(PASTIX_OTHER_LIB)

cpastix_dist:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_PASTIX) -Wno-deprecated $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(SCOTCH_INCLUDE) $(COMPLEX_PASTIX_INCLUDE) $(PASTIX_OTHER_INCLUDE) complex_dpastix_FreeFem.cpp
	$(MPIFC) $(FFLAGS) $(PTSCOTCH_INCLUDE) $(COMPLEX_PASTIX_INCLUDE) $(PASTIX_OTHER_INCLUDE) complex_dpastix_FreeFem.o -o complex_dpastix_FreeFem.$(SUF) $(MPI_LIB) $(BLAS_LIB) $(PTSCOTCH_LIB) $(METIS_LIB) $(COMPLEX_PASTIX_LIB) $(PASTIX_OTHER_LIB)

cpastix:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_PASTIX) -Wno-deprecated $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(SCOTCH_INCLUDE) $(COMPLEX_PASTIX_INCLUDE) $(PASTIX_OTHER_INCLUDE) complex_pastix_FreeFem.cpp
	$(MPIFC) $(FFLAGS) $(SCOTCH_INCLUDE) $(COMPLEX_PASTIX_INCLUDE) $(PASTIX_OTHER_INCLUDE) complex_pastix_FreeFem.o -o complex_pastix_FreeFem.$(SUF) $(MPI_LIB) $(BLAS_LIB) $(SCOTCH_LIB) $(METIS_LIB) $(COMPLEX_PASTIX_LIB) $(PASTIX_OTHER_LIB)

rsludist:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_SLUDIST) $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(SLUDIST_INCLUDE) $(SLUDIST_OTHER_INCLUDE) real_SuperLU_DIST_FreeFem.cpp
	$(MPIFC) $(FFLAGS)  $(SLUDIST_INCLUDE) $(SLUDIST_OTHER_INCLUDE) real_SuperLU_DIST_FreeFem.o -o real_SuperLU_DIST_FreeFem.$(SUF) $(MPI_LIB) $(BLAS_LIB) $(PARMETIS_LIB) $(METIS_LIB) $(SLUDIST_LIB) $(SLUDIST_OTHER_LIB)

csludist:
	$(MPICC) -c $(CFLAGS) $(CFLAGS_SLUDIST) $(MPI_INCLUDE) $(FREEFEM_INCLUDE) $(SLUDIST_INCLUDE) $(SLUDIST_OTHER_INCLUDE) complex_SuperLU_DIST_FreeFem.cpp
	$(MPIFC) $(FFLAGS)  $(SLUDIST_INCLUDE) $(SLUDIST_OTHER_INCLUDE) complex_SuperLU_DIST_FreeFem.o -o complex_SuperLU_DIST_FreeFem.$(SUF) $(MPI_LIB) $(BLAS_LIB) $(PARMETIS_LIB) $(METIS_LIB) $(SLUDIST_LIB) $(SLUDIST_OTHER_LIB)

