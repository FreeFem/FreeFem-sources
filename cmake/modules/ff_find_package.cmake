function(ff_find_package PACKAGE)

  find_package(${PACKAGE} ${ARGN})

  if(${PACKAGE}_FOUND)
  
    # Unifying variables names for include directory paths

    set(POSSIBLE_INCLUDE_PATH_VARIABLES "${PACKAGE}_INCLUDE"
                                        "${PACKAGE}_INCLUDES"
                                        "${PACKAGE}_INCLUDE_PATH"
                                        "${PACKAGE}_INCLUDE_PATHS"
                                        "${PACKAGE}_C_INCLUDE_PATHS")

    foreach(INCLUDE_PATH_VARIABLE IN LISTS ${POSSIBLE_INCLUDE_PATH_VARIABLES}) 
      list(APPEND INCLUDE_PATHS ${INCLUDE_PATH_VARIABLE})
    endforeach(INCLUDE_PATH_VARIABLE)

    # Unifying variables names for library directory paths

    set(POSSIBLE_LIBRARIES_VARIABLES "${PACKAGE}_LIBRARY"
                                     "${PACKAGE}_LIBRARIES"
                                     "${PACKAGE}_C_LIBRARY"
                                     "${PACKAGE}_C_LIBRARIES")

    foreach(LIBRARIES_VARIABLE IN LISTS ${POSSIBLE_LIBRARIES_VARIABLES}) 
      list(APPEND LIBRARIES ${LIBRARIES_VARIABLE})
    endforeach(LIBRARIES_VARIABLE)

    # Return unified names

    set(FF_${PACKAGE}_FOUND TRUE PARENT_SCOPE)
    set(FF_${PACKAGE}_INCLUDE_PATHS ${INCLUDE_PATHS} PARENT_SCOPE)
    set(FF_${PACKAGE}_LIBRARIES ${LIBRARIES} PARENT_SCOPE)

    add_library(FREEFEM::${PACKAGE} UNKNOWN IMPORTED)
    foreach(LIB ${LIBRARIES})
      set_target_properties(FREEFEM::${PACKAGE} PROPERTIES IMPORTED_LOCATION ${LIB})
    endforeach()
    set_target_properties(FREEFEM::${PACKAGE} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${INCLUDE_PATHS}")



  
  endif(${PACKAGE}_FOUND)

endfunction(ff_find_package)
