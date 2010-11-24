include @abs_top_builddir@/download/headers-sparsesolver.inc
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
CCD     = @CC@ -I@MPI_INC_DIR@
CFLAGS	= @CFLAGS@ $(FFINTSCOTCH) -Drestrict=__restrict -DCOMMON_PTHREAD -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_RANDOM_FIXED_SEED -DCOMMON_TIMING_OLD -DSCOTCH_RENAME -DSCOTCH_DETERMINISTIC #enlever -DSCOTCH_PTHREAD do not allow thread
LDFLAGS	= @LDFLAGS@ -lz -lm @LIBSPTHREAD@  #-L/usr/local/lib/gcc/
CP	= cp
LEX	= @LEX@
LN	= ln
MKDIR	= mkdir
MV	= mv
RANLIB	= @RANLIB@
YACC	= @YACC@
WGET    = @WGET@
