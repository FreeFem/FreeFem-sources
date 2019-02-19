# BLAS
#FFBLASINCLUDE = -I/Users/morice/librairie/PATCHVECLIB/
#FFBLASDIRLIBS = 
#FFBLASLIB     = -L/Users/morice/librairie/PATCHVECLIB/ -lwrapperdotblas -framework veclib
#FFBLASLIB2    = -lblas


# LAPACK
#FFLAPACKDIRLIBS = 
#FFLAPACKLIB     = -framework veclib
#FFLAPACKLIB2    = -llapack

# MPI
#FFMPIDIR      = /Users/morice/librairie/openmpi-gcc-gfortran-4.4/
#FFMPIINCLUDE  = -I/Users/morice/librairie/openmpi-gcc-gfortran-4.4/include/
#FFMPILIB      = -L/Users/morice/librairie/openmpi-gcc-gfortran-4.4/lib/ -lmpi -lmpi_cxx -lopen-pal -lopen-rte -lotf -lvt
#FFMPIDIRLIBS  = /Users/morice/librairie/openmpi-gcc-gfortran-4.4/lib/
#FFMPILIB2     = mpi mpi_cxx open-pal open-rte otf vt

# INT

#pastix
FFVERSIONINT = _long
FFCTYPESINT  = -DFORCE_LONG -DLONG

#scotch
FFINTSCOTCH  = -DLONG

#hips :: -DINTSIZE32, -DINTSIZE64 ou default
FFINTHIPS    = 

# particularite
#
FFLIBOTHERSMUMPS = -framework Carbon -framework AppKit 

#blasdef :: 
FFBLASDEF =  -DUSE_VENDOR_BLAS

#### parameter for blacs ####

#  ------------------------------------------------------------------------
#  Allows the user to vary the topologies that the BLACS default topologies
#  (TOP = ' ') correspond to.  If you wish to use a particular topology
#  (as opposed to letting the BLACS make the choice), uncomment the
#  following macros, and replace the character in single quotes with the
#  topology of your choice.
#  ------------------------------------------------------------------------
#  DEFBSTOP   = -DDefBSTop="'1'"
#  DEFCOMBTOP = -DDefCombTop="'1'"

FFDEFBSTOP   = 
FFDEFCOMBTOP =

#  -------------------------------------------------------------------
#  If your MPI_Send is locally-blocking, substitute the following line
#  for the empty macro definition below.
#  SENDIS = -DSndIsLocBlk
#  -------------------------------------------------------------------
FFSENDIS =

#  --------------------------------------------------------------------
#  If your MPI handles packing of non-contiguous messages by copying to
#  another buffer or sending extra bytes, better performance may be
#  obtained by replacing the empty macro definition below with the
#  macro definition on the following line.
#  BUFF = -DNoMpiBuff
#  --------------------------------------------------------------------
FFBUFF = 

#  -----------------------------------------------------------------------
#  If you know something about your system, you may make it easier for the
#  BLACS to translate between C and fortran communicators.  If the empty
#  macro defininition is left alone, this translation will cause the C
#  BLACS to globally block for MPI_COMM_WORLD on calls to BLACS_GRIDINIT
#  and BLACS_GRIDMAP.  If you choose one of the options for translating
#  the context, neither the C or fortran calls will globally block.
#  If you are using MPICH, or a derivitive system, you can replace the 
#  empty macro definition below with the following (note that if you let
#  MPICH do the translation between C and fortran, you must also indicate
#  here if your system has pointers that are longer than integers.  If so,
#  define -DPOINTER_64_BITS=1.)  For help on setting TRANSCOMM, you can
#  run BLACS/INSTALL/xtc_CsameF77 and BLACS/INSTALL/xtc_UseMpich as
#  explained in BLACS/INSTALL/README.
#   TRANSCOMM = -DUseMpich
#
#  If you know that your MPI uses the same handles for fortran and C
#  communicators, you can replace the empty macro definition below with
#  the macro definition on the following line.
#  TRANSCOMM = -DCSameF77
#  -----------------------------------------------------------------------
FFTRANSCOMM =  -DUseMpi2

#  --------------------------------------------------------------------------
#  You may choose to have the BLACS internally call either the C or Fortran77
#  interface to MPI by varying the following macro.  If TRANSCOMM is left
#  empty, the C interface BLACS_GRIDMAP/BLACS_GRIDINIT will globally-block if
#  you choose to use the fortran internals, and the fortran interface will
#  block if you choose to use the C internals.  It is recommended that the 
#  user leave this macro definition blank, unless there is a strong reason
#  to prefer one MPI interface over the other.
#  WHATMPI = -DUseF77Mpi
#  WHATMPI = -DUseCMpi
#  --------------------------------------------------------------------------
FFWHATMPI = 

#  ---------------------------------------------------------------------------
#  Some early versions of MPICH and its derivatives cannot handle user defined
#  zero byte data types.  If your system has this problem (compile and run
#  BLACS/INSTALL/xsyserrors to check if unsure), replace the empty macro
#  definition below with the macro definition on the following line.
#  SYSERRORS = -DZeroByteTypeBug
#  ---------------------------------------------------------------------------
FFSYSERRORS = 