function(ff_configure_cmake)

  # make sure that "make" will run in parallel 
  include(ProcessorCount)
  
  ProcessorCount(numprocs)
  
  if(NOT numprocs EQUAL 0)
    set(CMAKE_BUILD_FLAGS -j ${numprocs})
  endif()



endfunction(ff_configure_cmake)
