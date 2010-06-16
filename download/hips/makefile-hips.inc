DOWNLOADFF=$(abs_top_builddir)/download/
include $(abs_top_builddir)/download/headers-sparsesolver.inc
### ADD for freefem++
###
###  HIPS specific compilation flags
###

##   Arithmetic
##    - default is -DTYPE_REAL in double precision (-DPREC_DOUBLE)
##    - use -DTYPE_COMPLEX to build Complex version of HIPS
##    - use -DPREC_SIMPLE to compute in single precision

#COEFTYPE     = 

COEFTYPE    = -DTYPE_REAL
#COEFTYPE    = -DTYPE_COMPLEX

#COEFTYPE    = -DTYPE_REAL    -DPREC_SIMPLE
#COEFTYPE    = -DTYPE_COMPLEX -DPREC_SIMPLE


##   Partitionner
##    - default partitioner is METIS
##    - use -DSCOTCH_PART to use SCOTCH

PARTITIONER  = 
#PARTITIONER = -DSCOTCH_PART

##   Integer size
##    - default int type is    : INTS = INTL = int (C type length)
##    - use -DINTSIZE32 To set : INTS = INTEGER(4) and INTL = INTEGER(4)
##    - use -DINTSIZE64 to set : INTS = INTEGER(4) and INTL = INTEGER(8)

INTSIZE      =
#INTSIZE     = -DINTSIZE64
#INTSIZE     = -DINTSIZE32


###
###  Compiler
###

ARCH       = -DLINUX

CC	   = @CC@       # C compiler 
MPICC      = @MPICC@ 
CPP	   = cpp
FC         = 
MPIFC      = @MPIFC@
LD	   = $(CC)     # Linker
MPILD      = $(MPICC)

CFLAGS	   =  -std=c99       # Additional C compiler flags
FFLAGS	   =           # Additional Fortran compiler flags
LFLAGS     =           # Additional linker flags




COPTFLAGS  = -O3 #-fno-inline      # Optimization flags
FOPTFLAGS  = -O3 #-fno-inline      # 

###
###  Library
###

IBLAS      =  @BLASINC@   # BLAS include path
LBLAS      =  @BLASLIBS@   @LIBSPTHREAD@ # BLAS linker flags

IMPI       =  $(FFMPIINCLUDE)      # Additional MPI include path
LMPI       =  $(FFMPILIB)  -lsupc++ @FLIBS@  # Additional MPI linker flags

##   METIS_DIR : path to METIS
METIS_DIR  =  $(DOWNLOADFF)#
IMETIS     = -I$(METIS_DIR)include/metis/
LMETIS     = -L$(METIS_DIR)lib/parmetis/ -lmetis

##   SCOTCH_DIR : path to SCOTCH
SCOTCH_DIR = $(DOWNLOADFF)#/Users/morice/work/postdoc/freefem++prod/download/
ISCOTCH    = -I$(SCOTCH_DIR)include/scotch
LSCOTCH    = -L$(SCOTCH_DIR)lib/ -lscotch -lscotcherr

###
###  Misc
###

MAKE	   = make
AR	   = @AR@       #//ar
ARFLAGS	   = @ARFLAGS@ #//-crs  
LN	   = ln
CP	   = cp

###

##   Uncomment that to run in DEBUG mode
#DEBUG     = -g -DDEBUG_M

WGET    = @WGET@