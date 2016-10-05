FILE(GLOB CPPFILES *.cpp)
LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-load/bmo_user.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/dfft.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-Ipopt.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-NLopt.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ffrandom.cpp # config.h
                          ${CMAKE_SOURCE_DIR}/examples++-load/freeyams.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/gsl.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ilut.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MetricPkuate.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d-v4.0.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mshmet.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS_seq.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/netgen.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/NewSolver.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/PARDISO.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/scotch.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/SuperLu.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/tetgen.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/UMFPACK64.cpp)

FOREACH(CPPFILE ${CPPFILES})

  GET_FILENAME_COMPONENT(FILENAME ${CPPFILE} NAME_WE)
  ADD_LIBRARY(${FILENAME} ${CPPFILE})
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/bamglib
                      ${CMAKE_SOURCE_DIR}/src/fflib
                      ${CMAKE_SOURCE_DIR}/src/femlib
                      ${CMAKE_SOURCE_DIR}/src/Graphics
                      ${CMAKE_SOURCE_DIR}/src/libMesh)
ENDFOREACH(CPPFILE)

