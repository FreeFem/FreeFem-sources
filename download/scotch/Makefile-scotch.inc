# ALH - avoid blocking when Makefile.inc is not built yet
#include @abs_top_builddir@/download/headers-sparsesolver.inc

EXE	=
LIB	= .a
OBJ	= .o
#DIRMPI  = @FFMPIDIR@
#MAKE	= make
AR	= @AR@
ARFLAGS	= @ARFLAGS@
CAT	= cat
CCS     = @CC@ 
CCP	= @MPICC@ 
CCD     = @CC@ -I'@MPI_INC_DIR@'
FFINTSCOTCH = -DINTSIZE32 -DIDXSIZE32 -g
# FFCS - 16/4/13 - some plain C files seem to need access to mpi.h too

# FFCS - 23/4/13 - remove thread and fork dependency on Windows with
# -DCOMMON_STUB_FORK without -DCOMMON_PTHREAD
# [[file:scotch_5.1_esmumps/INSTALL.TXT::COMMON_STUB_FORK]]
# s/-DCOMMON_TIMING_OLD/-Du_int32_t='unsigned int'/" -e "s/-DCOMMON_RANDOM_FIXED_SEED/-DCOMMON_RANDOM_RAND  -DCOMMON_RANDOM_FIXED_SEED -D'pipe(pfds)=_pipe(pfds,1024,0x8000)'
CFLAGSWIN = -DHAVE_SYS_TIME_H -DCOMMON_STUB_FORK -Du_int32_t='unsigned int'  -DSCOTCH_RENAME -DCOMMON_RANDOM_RAND  -DCOMMON_RANDOM_FIXED_SEED -D'pipe(pfds)=_pipe(pfds,1024,0x8000)' 
CFLAGSSTD = -DCOMMON_STUB_FORK  -DCOMMON_RANDOM_FIXED_SEED -DCOMMON_TIMING_OLD -DSCOTCH_RENAME -DSCOTCH_DETERMINISTIC

ifeq (@WIN32DLLTARGET@,) 
CFLAGS	= @CFLAGS@ -I'@MPI_INC_DIR@' $(FFINTSCOTCH) -Drestrict=__restrict $(CFLAGSSTD)
else
CFLAGS	= @CFLAGS@ -I'@MPI_INC_DIR@' $(FFINTSCOTCH) -Drestrict=__restrict $(CFLAGSWIN)
endif




LDFLAGS	= @LDFLAGS@ -lm @MPI_LIBC@  @LIBSPTHREAD@  
#-L/usr/local/lib/gcc/
CP	= cp
LEX	= @LEX@
LN	= ln
MKDIR	= mkdir
MV	= mv
RANLIB	= @RANLIB@
YACC	= @YACC@
WGET    = @WGET@
