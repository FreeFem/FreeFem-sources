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
CC    = gcc 
MPICC = /usr/local/bin/mpicc

# C++ compiler
CXX    = g++
MPICXX = /usr/local/bin/mpicxx 

# Fortran compiler
FC     = fort77
MPIFC  = /usr/local/bin/mpif77 #mpifort77 #/Users/morice/librairie/openmpi/bin/mpif90 --prefix /Users/morice/librairie/openmpi/
FC_LIB = -L/usr/local/lib -lmpi_f77

# Linker
MPILD = 

# directory of include file for mpi
MPI_INCLUDE = /I/usr/local/include #-I/Users/morice/librairie/openmpi/include/
MPI_LIB     = -L/usr/local/lib -lmpi -lmpi_cxx -lmpi_f77 -lmpi_f90 -lopen-rte -lopen-pal

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
FREEFEM_DIR     = /Volumes/Ike/cazeaux/ff++/
FREEFEM_INCLUDE = -I$(FREEFEM_DIR)examples++-load/include/ 

FREEFEM_INCLUDE = -I$(FREEFEM_DIR)examples++-load/include/ 

# directory of the different package
#SOFT_DIR =  /Users/morice/librairie/
SOFT_DIR =  $(FREEFEM_DIR)/download/

#==========
#   BLAS
#==========
BLAS_INCLUDE = 
BLAS_LIB     = -L$(FREEFEM_DIR)ff++/src/solver/PATCHVECLIB/libwrapperdostblas.a  -framework veclib

#============
#  BLACS : this library is need by MUMPS 
#============
BLACS_DIR     =  $(SOFT_DIR)
BLACS_INCLUDE = 
BLACS_LIB     =  -L$(BLACS_DIR)lib/ -lblacsCinit_MPI-FREEFEM-0 -lblacs_MPI-FREEFEM-0 -lblacsF77init_MPI-FREEFEM-0

#=============
#  SCALAPACK : this library is need by MUMPS
#=============
SCALAPACK_DIR     = $(SOFT_DIR)
SCALAPACK_INCLUDE = -I$(SCALAPACK_DIR)include/ 
SCALAPACK_LIB     = -L$(SCALAPACK_DIR)lib/ -lscalapack

#==============================================
#   ordering : metis, parmetis, scotch 
#===============================================

# scotch
SCOTCH_DIR = $(SOFT_DIR)
SCOTCH_INCLUDE = -I$(SCOTCH_DIR)include/
SCOTCH_LIB     = -L$(SCOTCH_DIR)lib/  -lesmumps  -lscotch -lscotcherr -lscotcherrexit -lm

# ptscotch 
PTSCOTCH_INCLUDE = -I$(SCOTCH_DIR)include/
PTSCOTCH_LIB     = -L$(SCOTCH_DIR)lib/ -lptscotch -lptscotcherrexit -lm

# metis
METIS_DIR      = $(SOFT_DIR)
METIS_INCLUDE  = -I$(METIS_DIR)include/metis/
METIS_LIB      = -L$(METIS_DIR)lib -lmetis

# parmetis
PARMETIS_DIR      = $(SOFT_DIR)
PARMETIS_INCLUDE  = -I$(PARMETIS_DIR)include/
PARMETIS_LIB      = -L$(PARMETIS_DIR)lib/ -lparmetis -lmetis

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
HIPS_DIR           = $(SOFT_DIR)
#    include file of HIPS
HIPS_INCLUDE       = -I$(HIPS_DIR)include/hips/
#   library of package hips
HIPS_ALL_LIBS      = -L$(HIPS_DIR)lib/hips -lio -lhipssequential -lhips #-lspkit 

IO_LIB             = $(HIPS_DIR)lib/hips/libio.a
HIPSSEQUENTIAL_LIB = $(HIPS_DIR)lib/hips/libhipssequential.a
HIPS_LIB           = $(HIPS_DIR)lib/hips/libhips.a
#SPKIT_LIB	   = $(HIPS_DIR)LIB/libspkit.a


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
HYPRE_DIR     = $(SOFT_DIR)
# include file of HYPRE
HYPRE_INCLUDE = -I$(HYPRE_DIR)include/hypre
# library of HYPRE
HYPRE_LIB     = -L$(HYPRE_DIR)lib/hypre/ -lHYPRE

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
PARMS_DIR    = $(SOFT_DIR)
# include file of pARMS
PARMS_INCLUDE =  -I$(PARMS_DIR)include/parms
# library of HYPRE
PARMS_LIB     =  -L$(PARMS_DIR)lib/ -lparms-O3

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
MUMPS_DIR     =  $(SOFT_DIR)mumps/MUMPS_4.9.2/
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
SLUDIST_DIR = $(SOFT_DIR)

# include file of superlu_dist
SLUDIST_INCLUDE = -I$(SLUDIST_DIR)include 
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
REAL_PASTIX_DIR     = $(SOFT_DIR)/
# include directory of pastix
REAL_PASTIX_INCLUDE = -I$(REAL_PASTIX_DIR)include/pastix/double/ 
# library
REAL_PASTIX_LIB     = -L$(REAL_PASTIX_DIR)lib/pastix/double -lpastix #-lmatrix_driver 

# complex pastix
#  directory of pstix
COMPLEX_PASTIX_DIR     = $(SOFT_DIR)
# include directory of pastix
COMPLEX_PASTIX_INCLUDE = -I$(COMPLEX_PASTIX_DIR)include/pastix/complex/
# library
COMPLEX_PASTIX_LIB     = -L$(COMPLEX_PASTIX_DIR)lib/pastix/complex/ -lpastix #-lmatrix_driver 


# other compiler option for PASTIX (see value of CFLAGS in SECTION 1)
CFLAGS_FLAGS        =   -DFORCE_LONG -DNUMA_ALLOC -DX_ARCHi686_mac -DDOF_CONSTANT  -DCOMMON_DEBUG -DSYMBOL_DEBUG -DDOF_DEBUG -DSOLVER_DEBUG -DORDER_DEBUG -DFAX_DEBUG -DSCOTCH_DEBUG -DDEBUG_BLEND -DBUBBLE_DEBUG -DPASTIX_DEBUG -Wno-depreciated
# other directory of header file and library file for interfacing
PASTIX_OTHER_INCLUDE = 
PASTIX_OTHER_LIB     = -lm -lsupc++ -lpthread