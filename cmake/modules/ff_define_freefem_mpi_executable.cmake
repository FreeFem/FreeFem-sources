macro(ff_define_freefem_mpi_executable)

  # Steps are the same as those required for FreeFem executable
  # (see ff_define_freefem_executable for an explanation)

  if(MPI_FOUND)
    
    add_executable(FreeFem++-mpi ${CMAKE_SOURCE_DIR}/src/mpi/parallelempi.cpp
                                 ${CMAKE_SOURCE_DIR}/src/Graphics/sansrgraph.cpp
                                 ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp)
  
    target_compile_definitions(FreeFem++-mpi PRIVATE PARALLELE)

    include_directories(${MPI_CXX_INCLUDE_PATH} 
                        ${CMAKE_SOURCE_DIR}/src/femlib)
  
    target_link_libraries(FreeFem++-mpi lglib 
                                        libff 
                                        ${MPI_CXX_LIBRARIES})

    install(TARGETS FreeFem++-mpi
            RUNTIME DESTINATION bin)

  endif()

endmacro()
