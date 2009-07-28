#// ORIG-DATE: 02/2009
#// -*- Mode : c++ -*-
#//
#// SUMMARY  :  
#// USAGE    : LGPL      
#// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
#// AUTHOR   : Jacques Morice
#// E-MAIL   : jacques.morice@ann.jussieu.fr
#//
#/* 
# This file is part of Freefem++
# 
# Freefem++ is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 2.1 of the License, or
# (at your option) any later version.
# 
# Freefem++  is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with Freefem++; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Thank to the ARN ()  FF2A3 grant
# ref:ANR-07-CIS7-002-01 
# */
#####################################################
#    SECTION 1: 
#    General variables to compile sparse solver  
#####################################################
# Note ::   A remplir automatique avec automake
#           Regarder notation de FreeFem++
###################
#    Compiler 
###################
# C compiler
CC    = gcc-4.4 
MPICC = /Users/morice/librairie/openmpi/bin/mpicc --prefix /Users/morice/librairie/openmpi/

# C++ compiler
CXX    = g++4.4
MPICXX = /Users/morice/librairie/openmpi/bin/mpicxx --prefix /Users/morice/librairie/openmpi/

# Fortran compiler
FC     = f90
MPIFC  = /Users/morice/librairie/openmpi/bin/mpif90 --prefix /Users/morice/librairie/openmpi/
FC_LIB = 

# Linker
MPILD = 

# directory of include file for mpi
MPI_INCLUDE = -I/Users/morice/librairie/openmpi/include/
MPI_LIB     = -L/Users/morice/librairie/openmpi/lib/ -lmpi -lmpi_cxx -lmpi_f77 -lmpi_f90 -lopen-rte -lopen-pal

###############
#   FLAGS
###############
# 
LDFLAGS_FREEFEM = -bundle -undefined dynamic_lookup -DNDEBUG -g -DDRAWING -DBAMG_LONG_LONG -DCHECK_KN -fno-inline -fexceptions -I/usr/X11/include/
CFLAGS_FREEFEM  =  -DNDEBUG -g -DDRAWING -DBAMG_LONG_LONG -DCHECK_KN -fno-inline -fexceptions -I/usr/X11/include/
#
CFLAGS	   =  $(CFLAGS_FREEFEM) # Additional C compiler flags
FFLAGS	   =  $(LDFLAGS_FREEFEM) # Additional Fortran compiler flags
LDFLAGS    =  $(LDFLAGS_FREEFEM) # linker flags
COPTFLAGS  =    -O3       # Optimization flags
FOPTFLAGS  =    -O3       # 
LDOPTFLAGS =    -O3       #

# SUFFIX OF DYNAMIC LIBRARY
# linux   SUF = so
# windows SUF = dll
# mac osx SUF = dylib
SUF      =  dylib


##########################################################
# SECTION 2: 
# common variable for the different parallel sparse solver
#
#########################################################

# frefem++ variables
FREEFEM_DIR     = /Users/morice/work/postdoc/freefem++debug/
FREEFEM_INCLUDE = -I$(FREEFEM_DIR)examples++-load/include/ 

# directory of the different package
SOFT_DIR =  /Users/morice/librairie/

#==========
#   BLAS
#==========
BLAS_INCLUDE = 
BLAS_LIB     = -L$(SOFT_DIR)PATCHVECLIB/ -lwrapperdotblas  -framework veclib

#============
#  BLACS : this library is need by MUMPS 
#============
BLACS_DIR     =  $(SOFT_DIR)/BLACS/
BLACS_INCLUDE = 
BLACS_LIB     =  -L$(BLACS_DIR)LIB/ -lblacsCinit_MPI-BROCHET-0 -lblacs_MPI-BROCHET-0 -lblacsF77init_MPI-BROCHET-0

#=============
#  SCALAPACK : this library is need by MUMPS
#=============
SCALAPACK_DIR     = $(SOFT_DIR)/scalapack-1.8.0/
SCALAPACK_INCLUDE = -I$(SCALAPACK_DIR)SRC/ 
SCALAPACK_LIB     = -L$(SCALAPACK_DIR) -lscalapack

#==============================================
#   ordering : metis, parmetis, scotch 
#===============================================
#  pord : this ordering is given in the mumps package 
#         and it is only need by MUMPS
#PORD_DIR     =  $(SOFT_DIR)MUMPS_4.8.3/PORD/
#PORD_INCLUDE = -I$(PORD_DIR)include
#PORD_LIB     = -L$(PORD_DIR)lib/ -lpord

# scotch
SCOTCH_DIR = $(SOFT_DIR)/scotch_5.1_esmumps/
SCOTCH_INCLUDE = -I$(SCOTCH_DIR)include/
SCOTCH_LIB     = -L$(SCOTCH_DIR)lib/  -lesmumps  -lscotch -lscotcherr -lscotcherrexit -lm

# ptscotch 
PTSCOTCH_INCLUDE = -I$(SCOTCH_DIR)include/
PTSCOTCH_LIB     = -L$(SCOTCH_DIR)lib/ -lptscotch -lptscotcherrexit -lm

# metis
METIS_DIR      = $(SOFT_DIR)ParMetis-3.1.1/
METIS_INCLUDE  = -I$(METIS_DIR)METISLib/
METIS_LIB      = -L$(METIS_DIR) -lmetis

# parmetis
PARMETIS_DIR      = $(SOFT_DIR)ParMetis-3.1.1/
PARMETIS_INCLUDE  = -I$(PARMETIS_DIR)
PARMETIS_LIB      = -L$(PARMETIS_DIR) -lparmetis -lmetis

###########################################
#  SECTION 3: 
#      ITERATIVE/PRECONDITIONNER SPARSE SOLVER
#  
###########################################
#      HIPS
###########################################
#  To compile HIPS, we need the following librairies
#     -- ordering Metis or Scotch
#        comment lines corresponding to Metis or Scotch in SECTION 2
#     -- BLAS
#        comment lines corresponding to BLAS in SECTION 2
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  Variables needs by HIPS
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     directory of HIPS
HIPS_DIR           = /Users/morice/librairie/iterativesolver/hips-1.2beta-rc3/
#    include file of HIPS
HIPS_INCLUDE       = -I$(HIPS_DIR)SRC/INCLUDE
#   library of package hips
HIPS_ALL_LIBS      = -L$(HIPS_DIR)LIB -lio -lhipssequential -lhips -lspkit 

IO_LIB             = $(HIPS_DIR)LIB/libio.a
HIPSSEQUENTIAL_LIB = $(HIPS_DIR)LIB/libhipssequential.a
HIPS_LIB           = $(HIPS_DIR)LIB/libhips.a
SPKIT_LIB	   = $(HIPS_DIR)LIB/libspkit.a


# other compiler option for HIPS (see value of CFLAGS in SECTION 1)
CFLAGS_HIPS =
#  other include file and library file for interfacing
HIPS_OTHER_INCLUDE =
HIPS_OTHER_LIB     = 

###########################################
#      HYPRE
###########################################
#  To compile HYPRE, we need the following librairies
#     -- ordering Metis 
#        comment lines corresponding to Metis in SECTION 2
#     -- BLAS
#        comment lines corresponding to BLAS in SECTION 2
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    Variables needs by HYPRE
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# directory of HYPRE
HYPRE_DIR     = /Users/morice/librairie/iterativesolver/hypre-2.4.0b/
# include file of HYPRE
HYPRE_INCLUDE = -I$(HYPRE_DIR)src/hypre/include/
# library of HYPRE
HYPRE_LIB     = -L$(HYPRE_DIR)src/lib/ -lHYPRE

# other compiler option for HYPRE (see value of CFLAGS in SECTION 1)
CFLAGS_HYPRE        = -Wall -DHAVE_CONFIG_H -DHYPRE_TIMING -Wno-deprecated -DLINUX -DSPKIT_F  -DMPICH_SKIP_MPICXX
#  other directory of header file and library file for interfacing
HYPRE_OTHER_INCLUDE = 
HYPRE_OTHER_LIB     = -lm -lstdc++

###########################################
#      pARMS
###########################################
#  To compile pARMS, we need the following librairies
#     -- ordering Metis 
#        comment lines corresponding to Metis in SECTION 2
#     -- BLAS
#        comment lines corresponding to BLAS in SECTION 2
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   variable needs by pARMS
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# directory of pARMS
PARMS_DIR    = /Users/morice/librairie/iterativesolver/pARMS_2.2/
# include file of pARMS
PARMS_INCLUDE =  -I$(PARMS_DIR)INCLUDE/
# library of HYPRE
PARMS_LIB     =  -L$(PARMS_DIR)LIB/ -lparms-O3

# other compiler option for HIPS (see value of CFLAGS in )
# PARMS NEED THE ARCHITECTURE (see generaldefs.h). For Mac Osx defined the architecture as -DLINUX
CFLAGS_PARMS   = -DLINUX -DSPKIT_F 
# other include file and library file for interfacing
PARMS_OTHER_INCLUDE = 
PARMS_OTHER_LIB     = 


###########################################
#  SECTION 4: 
#      DIRECT SPARSE SOLVER
#  
###########################################
#      MUMPS
###########################################
#  To compile MUMPS, we need the following librairies
#     -- ordering metis, scotch 
#        comment lines corresponding to Metis and Scotch in SECTION 2
#     -- BLAS, BLACS, SCALAPACK
#        comment lines corresponding to BLAS, BLACS and SCALAPACK in SECTION 2
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    Variables needs by MUMPS
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MUMPS_HOST  = BROCHET
#  directory of MUMPS
MUMPS_DIR     =  /Users/morice/librairie/MUMPS_4.8.3/
# directory of Header File for MUMPS
MUMPS_INCLUDE = -I$(MUMPS_DIR)include/
# library of MUMPS
MUMPS_LIB     = -L$(MUMPS_DIR)lib/ -lmumps_common_$(MUMPS_HOST) -ldmumps_$(MUMPS_HOST) -lzmumps_$(MUMPS_HOST) 

#   pord : ordering include in MUMPS package 
PORD_DIR     =  $(MUMPS_DIR)PORD/
PORD_INCLUDE = -I$(PORD_DIR)include
PORD_LIB     = -L$(MUMPS_DIR)lib/ -lpord_$(MUMPS_HOST) 

# other compiler option for MUMPS (see value of CFLAGS in SECTION 1)
CFLAGS_MUMPS  =
#  other include file and library file for interfacing
MUMPS_OTHER_INCLUDE = -I/usr/X11/include
MUMPS_OTHER_LIB     = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1/ -lgfortran

###########################################
#      SuperLU_DIST
###########################################
#  To compile SuperLU_DIST, we need the following librairies
#     -- ordering Metis and ParMetis 
#        comment lines corresponding to Metis and ParMetis in SECTION 2
#     -- BLAS
#        comment lines corresponding to BLAS in SECTION 2
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    Variables needs by SuperLU_DIST
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# directory of superlu_dist
SLUDIST_DIR = /Users/morice/librairie/SuperLU_DIST_2.3/
# include file of superlu_dist
SLUDIST_INCLUDE = -I$(SLUDIST_DIR)SRC/ 
# librairie of superlu_dist
SLUDIST_LIB = -L$(SLUDIST_DIR)lib/ -lsuperlu_dist_2.3

# other compiler option for SuperLU (see value of CFLAGS in SECTION 1)
SLUDIST_CFLAGS        = 
#  other include file and library file for interfacing
SLUDIST_OTHER_INCLUDE = 
SLUDIST_OTHER_LIB     = 

###########################################
#      Pastix
###########################################
#  To compile pastix, we need the following librairies
#     -- ordering scotch, ptscotch, metis 
#        comment lines corresponding to metis, scotch and ptsccoth in SECTION 2
#     -- BLAS 
#        comment lines corresponding to BLAS in SECTION 2
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    Variables needs by pastix
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  To obtain real and complex  arithmetics, we need to install for real and for complex separately. 
#   Hence we need to have two differents pastix directory
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
#  real pastix
#  directory of pastix
REAL_PASTIX_DIR     = /Users/morice/librairie/pastix_release_2259/
# include directory of pastix
REAL_PASTIX_INCLUDE = -I$(REAL_PASTIX_DIR)install-double-nondist/ 
# library
REAL_PASTIX_LIB     = -L$(REAL_PASTIX_DIR)install-double-nondist/ -lpastix #-lmatrix_driver 

# complex pastix
#  directory of pstix
COMPLEX_PASTIX_DIR     = /Users/morice/librairie/pastix_release_2259/
# include directory of pastix
COMPLEX_PASTIX_INCLUDE = -I$(COMPLEX_PASTIX_DIR)install-complex-nondist/
# library
COMPLEX_PASTIX_LIB     = -L$(COMPLEX_PASTIX_DIR)install-complex-nondist/ -lpastix #-lmatrix_driver 


# other compiler option for PASTIX (see value of CFLAGS in SECTION 1)
CFLAGS_FLAGS        =   -DFORCE_LONG -DNUMA_ALLOC -DX_ARCHi686_mac -DDOF_CONSTANT  -DCOMMON_DEBUG -DSYMBOL_DEBUG -DDOF_DEBUG -DSOLVER_DEBUG -DORDER_DEBUG -DFAX_DEBUG -DSCOTCH_DEBUG -DDEBUG_BLEND -DBUBBLE_DEBUG -DPASTIX_DEBUG -Wno-depreciated
# other directory of header file and library file for interfacing
PASTIX_OTHER_INCLUDE = 
PASTIX_OTHER_LIB     = -lm -lsupc++ -lpthread