macro(ff_configure_thirdparty_optional)
  
  include(ff_find_package)

  set(FF_THIRDPARTY_OPTIONAL CBLAS
                             FFTW
                             IPOPT
                             LAPACK
                             METIS
                             MPI
                             MUMPS
                             NLOPT
                             SCOTCH
                             SUPERLU
                             TETGEN
                             UMFPACK
                             Threads)

  foreach(PACKAGE ${FF_THIRDPARTY_OPTIONAL})
    ff_find_package(${PACKAGE})
  endforeach(PACKAGE ${FF_THIRDPARTY_OPTIONAL})

endmacro(ff_configure_thirdparty_optional)
