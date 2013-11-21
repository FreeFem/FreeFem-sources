# ALH - avoid blocking when Makefile.inc is not built yet
-include @abs_top_builddir@/download/headers-sparsesolver.inc

EXE	=
LIB	= .a
OBJ	= .o
DIRMPI  = @FFMPIDIR@
MAKE	= make
AR	= @AR@
ARFLAGS	= @ARFLAGS@
CAT	= cat
CCS     = @CC@ 
CCP	= @MPICC@ 
CCD     = @CC@ -I'@MPI_INC_DIR@'

# FFCS - 16/4/13 - some plain C files seem to need access to mpi.h too

# FFCS - 23/4/13 - remove thread and fork dependency on Windows with
# -DCOMMON_STUB_FORK without -DCOMMON_PTHREAD
# [[file:scotch_5.1_esmumps/INSTALL.TXT::COMMON_STUB_FORK]]

CFLAGS	= @CFLAGS@ -I'@MPI_INC_DIR@' $(FFINTSCOTCH) -Drestrict=__restrict -DCOMMON_STUB_FORK  -DCOMMON_RANDOM_FIXED_SEED -DCOMMON_TIMING_OLD -DSCOTCH_RENAME -DSCOTCH_DETERMINISTIC

LDFLAGS	= @LDFLAGS@ -lm @MPI_LIB@  @LIBSPTHREAD@  #-L/usr/local/lib/gcc/
CP	= cp
LEX	= @LEX@
LN	= ln
MKDIR	= mkdir
MV	= mv
RANLIB	= @RANLIB@
YACC	= @YACC@
WGET    = @WGET@
