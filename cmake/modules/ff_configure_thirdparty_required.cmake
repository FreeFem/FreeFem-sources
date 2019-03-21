 macro(ff_configure_thirdparty_required)

  include(ff_find_package)

  set(FF_THIRDPARTY_REQUIRED ARPACK
                             DLOPEN
                             FLEX
                             OpenGL
                             GLUT
                             GSL
                             HDF5)

  # SuiteSparse is systematically downloaded on Windows
  # because the associated pacman package
  # is not compiled in the expected way
  if(MINGW)
    list(APPEND FF_THIRDPARTY_DOWNLOAD SUITESPARSE)
  else()
    list(APPEND FF_THIRDPARTY_REQUIRED SUITESPARSE)
  endif(MINGW)
  
  foreach(PACKAGE ${FF_THIRDPARTY_REQUIRED})
    ff_find_package(${PACKAGE} REQUIRED)
    if(NOT FF_${PACKAGE}_FOUND)
      message(SEND_ERROR "Required package ${PACKAGE} is missing")
    endif(NOT FF_${PACKAGE}_FOUND)
  endforeach(PACKAGE ${FF_THIRDPARTY_REQUIRED})

  ff_install_package(SUITESPARSE)

endmacro(ff_configure_thirdparty_required)
