 macro(ff_configure_thirdparty_required)

  include(ff_find_package)

  set(FF_THIRDPARTY_REQUIRED AMD
                             ARPACK
                             CAMD
                             CCOLAMD
                             CHOLMOD
                             COLAMD
                             DLOPEN
                             FLEX
                             #OpenGL
                             #GLUT
                             GSL
                             HDF5
                             SUITESPARSE
                             SUITESPARSECONFIG)

  foreach(PACKAGE ${FF_THIRDPARTY_REQUIRED})
    ff_find_package(${PACKAGE} REQUIRED)
    if(NOT FF_${PACKAGE}_FOUND)
      message(SEND_ERROR "Required package ${PACKAGE} is missing")
    endif(NOT FF_${PACKAGE}_FOUND)
  endforeach(PACKAGE ${FF_THIRDPARTY_REQUIRED})

endmacro(ff_configure_thirdparty_required)
