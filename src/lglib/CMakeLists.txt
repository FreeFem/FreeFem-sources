# TODO: remove this line as soon as autoconf is removed (only here for config.h and config-wrapper.h)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}) 


FILE(GLOB_RECURSE CPPFILES *.cpp)
ADD_LIBRARY(lglib ${CPPFILES})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/bamglib
                    ${CMAKE_SOURCE_DIR}/src/femlib
                    ${CMAKE_SOURCE_DIR}/src/fflib
                    ${CMAKE_SOURCE_DIR}/src/Graphics
                    ${CMAKE_SOURCE_DIR}/src/lglib/)
SET_TARGET_PROPERTIES(lglib PROPERTIES PREFIX "")

IF(FF_ARPACK_FOUND)
  TARGET_COMPILE_DEFINITIONS(lglib PRIVATE -DHAVE_LIBARPACK)
  if(MINGW)
    TARGET_LINK_LIBRARIES(lglib FREEFEM::ARPACK FREEFEM::LAPACK gfortran)
  else()
    TARGET_LINK_LIBRARIES(lglib FREEFEM::ARPACK FREEFEM::LAPACK)
  endif()
ENDIF()


