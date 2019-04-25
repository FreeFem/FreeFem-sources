macro(ff_define_libff_library)

  include(ff_define_strversionnumber_library)
  ff_define_strversionnumber_library()

  # Get paths to cpp files in femlib and fflib and put them in LIBFF_SRC
  file(GLOB FEMLIB_SRC ${CMAKE_SOURCE_DIR}/src/femlib/*.cpp)
  file(GLOB_RECURSE FFLIB_SRC ${CMAKE_SOURCE_DIR}/src/fflib/*.cpp)

  set(LIBFF_SRC ${FEMLIB_SRC} ${FFLIB_SRC})

  # Remove spurious cpp files from LIBFF_SRC
  list(REMOVE_ITEM LIBFF_SRC ${CMAKE_SOURCE_DIR}/src/bamglib/Meshgibbs.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/ConjuguedGradrientNL.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/FESpace-v0.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/glutdraw.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/InvIntFunc.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/mortar.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/Pkorder.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/P3korder.cpp
                             ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp
                             ${CMAKE_SOURCE_DIR}/src/fflib/strversionnumber.cpp)

  # Add other required cpp files to LIBFF_SRC
  list(APPEND LIBFF_SRC ${CMAKE_SOURCE_DIR}/src/Algo/lgalgo.cpp 
                        ${CMAKE_SOURCE_DIR}/src/Eigen/eigenvalue.cpp
                        ${CMAKE_SOURCE_DIR}/src/femlib/libmesh5.c
                        ${CMAKE_SOURCE_DIR}/src/Graphics/DefColor.cpp)

  # Set a required definition for C++ compilers on systems other than Windows/MINGW
  if(DLOPEN_FOUND AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows" AND NOT MINGW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D HAVE_DLFCN_H")
  endif()


  # Directories containing required headers
  include_directories(${CMAKE_SOURCE_DIR}/src/bamglib 
                      ${CMAKE_SOURCE_DIR}/src/fflib/ 
                      ${CMAKE_SOURCE_DIR}/src/Graphics/ 
                      ${CMAKE_SOURCE_DIR}/src/lglib/
                      ${CMAKE_SOURCE_DIR}/src/femlib)

  # Definition of libff library
  add_library(libff ${LIBFF_SRC})


  if(FF_SUITESPARSE_FOUND AND FF_CHOLMOD_FOUND)


    # Compilation definitions 

    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

      TARGET_COMPILE_DEFINITIONS(libff PRIVATE HAVE_LIBUMFPACK 
                                       PRIVATE HAVE_SUITESPARSE_UMFPACK_H)

    else()

      TARGET_COMPILE_DEFINITIONS(libff PRIVATE HAVE_LIBUMFPACK 
                                       PRIVATE HAVE_UMFPACK_H)

    endif()

    # Required libraries for linking libff 
    
    if(MINGW)
      
      find_package(OpenMP)
      string(REPLACE "C:/msys64" "" OpenMP_C_LIBRARIES "${OpenMP_C_LIBRARIES}")
      target_link_libraries(libff bamglib Comdlg32 strversionnumber FREEFEM::SUITESPARSE
                                                                    FREEFEM::CHOLMOD 
                                                                    FREEFEM::AMD 
                                                                    FREEFEM::CAMD
                                                                    FREEFEM::CCOLAMD
                                                                    FREEFEM::COLAMD 
                                                                    FREEFEM::METIS
                                                                    FREEFEM::SUITESPARSECONFIG 
                                                                    ${FF_LAPACK_LIBRARIES}
                                                                    ${OpenMP_C_LIBRARIES})
    else()
    
      target_link_libraries(libff bamglib dl strversionnumber FREEFEM::SUITESPARSE FREEFEM::CHOLMOD)
  
    endif(MINGW)

  else()

    if(MINGW)
      target_link_libraries(libff bamglib Comdlg32 strversionnumber)
    else()
      target_link_libraries(libff bamglib dl strversionnumber)
    endif(MINGW)

  endif() 


  # Remove lib prefix from the name of the library
  # (libff instead of liblibff)
  set_target_properties(libff PROPERTIES PREFIX "")


endmacro()
