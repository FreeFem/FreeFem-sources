function(ff_write_config_file)

  # Filename
  set(FILENAME "cmake-config.inc")

  # Header
  file(WRITE ${CMAKE_BINARY_DIR}/${FILENAME} "# Paths found by CMake\n")

  # Compilers
  file(APPEND ${CMAKE_BINARY_DIR}/${FILENAME} "CMAKE_C_COMPILER ${CMAKE_C_COMPILER}\n")
  file(APPEND ${CMAKE_BINARY_DIR}/${FILENAME} "CMAKE_C_FLAGS ${CMAKE_C_FLAGS}\n")
  file(APPEND ${CMAKE_BINARY_DIR}/${FILENAME} "CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER}\n")
  file(APPEND ${CMAKE_BINARY_DIR}/${FILENAME} "CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}\n")

  # Third-party libraries
  foreach(PACKAGE ${MODULE_LIST})
    file(APPEND ${CMAKE_BINARY_DIR}/${FILENAME} "FREEFEM_${PACKAGE}_INCLUDE_PATHS ${FREEFEM_${PACKAGE}_INCLUDE_PATHS}\n") 
    file(APPEND ${CMAKE_BINARY_DIR}/${FILENAME} "FREEFEM_${PACKAGE}_LIBRARIES ${FREEFEM_${PACKAGE}_LIBRARIES}\n") 
  endforeach(PACKAGE)



endfunction(ff_write_config_file)

