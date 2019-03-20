macro(ff_configure_thirdparty)

  include(ff_find_package)
  include(ff_install_package)

  # Required packages

  set(FF_THIRDPARTY_REQUIRED ARPACK
                             DLOPEN
                             FLEX
                             OpenGL
                             GLUT
                             GSL
                             HDF5)

  foreach(PACKAGE ${FF_THIRDPARTY_REQUIRED})
    ff_find_package(${PACKAGE} REQUIRED)
    if(NOT FF_${PACKAGE}_FOUND)
      message(SEND_ERROR "Required package ${PACKAGE} is missing")
    endif(NOT FF_${PACKAGE}_FOUND)
  endforeach(PACKAGE ${FF_THIRDPARTY_REQUIRED})
  

  # Optional packages

  set(FF_THIRDPARTY_OPTIONAL CBLAS
                             CHOLMOD
                             FFTW
                             IPOPT
                             LAPACK
                             METIS
                             MPI
                             MUMPS
                             NLOPT
                             OpenMP
                             SCOTCH
                             SUPERLU
                             TETGEN
                             UMFPACK
                             Threads)

  foreach(PACKAGE ${FF_THIRDPARTY_OPTIONAL})
    ff_find_package(${PACKAGE} REQUIRED)
  endforeach(PACKAGE ${FF_THIRDPARTY_OPTIONAL})

  


  if(WITH_PETSC)
    include(FFInstallPackage)
    FF_INSTALL_PACKAGE(PETSC)
  endif(WITH_PETSC)

  foreach(PACKAGE ${FF_THIRDPARTY})
    ff_find_package(${PACKAGE})
    if(NOT FF_${PACKAGE}_FOUND)
      list(APPEND FF_THIRDPARTY_MISSING ${PACKAGE})
    endif(NOT FF_${PACKAGE}_FOUND)
  endforeach(PACKAGE ${FF_THIRDPARTY})

  message(STATUS "The following packages are missing: ${FF_THIRDPARTY_MISSING}")


  list(APPEND MODULE_LIST CHOLMOD
                          FFTW
                          #GMM
                          IPOPT
                          METIS
                          MUMPS
                          NLOPT
                          SCOTCH
                          SUPERLU
                          TETGEN)

  # SuiteSparse is systematically downloaded on Windows
  # because the associated pacman package
  # is not compiled in the expected way
  if(NOT MINGW)
    list(APPEND MODULE_LIST SUITESPARSE)
  endif(NOT MINGW)


  foreach(MODULE ${MODULE_LIST})
#    ff_find_package(${MODULE})
    if(NOT FREEFEM_${MODULE}_INSTALLED)
      list(APPEND DOWNLOAD_LIST ${MODULE})
    endif(NOT FREEFEM_${MODULE}_INSTALLED)
  endforeach(MODULE)

  if(ENABLE_DOWNLOAD)
    message(STATUS "The following modules will be downloaded: ${DOWNLOAD_LIST}")
  endif(ENABLE_DOWNLOAD)


endmacro(ff_configure_thirdparty)
