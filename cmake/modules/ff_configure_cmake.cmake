function(ff_configure_cmake)

  include(ProcessorCount)
  
  ProcessorCount(numprocs)
  
  if(NOT numprocs EQUAL 0)
    SET(CMAKE_BUILD_FLAGS -j ${numprocs})
  ENDIF()


endfunction(ff_configure_cmake)
