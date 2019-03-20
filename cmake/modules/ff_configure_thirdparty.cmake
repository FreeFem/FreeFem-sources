macro(ff_configure_thirdparty)

  include(ff_find_package)


  ff_find_package(ARPACK REQUIRED)
  find_package(DLOPEN REQUIRED)
  find_package(FLEX REQUIRED)
  find_package(OpenGL)
  find_package(GLUT REQUIRED)
  find_package(GSL REQUIRED)
  find_package(HDF5 REQUIRED)
  
  find_package(CBLAS)
  find_package(LAPACK)
  find_package(MPI)
  find_package(Threads)
  find_package(UMFPACK)

  if(WITH_PETSC)
    include(FFInstallPackage)
    FF_INSTALL_PACKAGE(PETSC)
  endif(WITH_PETSC)


  list(APPEND MODULE_LIST CHOLMOD
                          FFTW
                          #GMM
                          IPOPT
                          METIS
                          MUMPS
                          NLOPT
                          SCOTCH
                          SUITESPARSE
                          SUPERLU
                          TETGEN)


  foreach(MODULE ${MODULE_LIST})
    ff_find_package(${MODULE})
    if(NOT FREEFEM_${MODULE}_INSTALLED)
      list(APPEND DOWNLOAD_LIST ${MODULE})
    endif(NOT FREEFEM_${MODULE}_INSTALLED)
  endforeach(MODULE)

#  if(ENABLE_DOWNLOAD)
#    message(STATUS "The following modules will be downloaded: ${DOWNLOAD_LIST}")
#  endif(ENABLE_DOWNLOAD)


endmacro(ff_configure_thirdparty)
