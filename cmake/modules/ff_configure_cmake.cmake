function(ff_configure_cmake)

  # make sure that "make" will run in parallel 
  include(ProcessorCount)
  
  ProcessorCount(numprocs)
  
  if(NOT numprocs EQUAL 0)
    set(CMAKE_BUILD_FLAGS -j ${numprocs})
  endif()

  if(NOT ENABLE_DOWNLOAD)
    set(ENABLE_DOWNLOAD False)
  endif(NOT ENABLE_DOWNLOAD)

endfunction(ff_configure_cmake)
