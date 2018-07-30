FILE(GLOB EDP_FILES *.edp)

LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-3d/fallingspheres.edp                # depends on mmg3d-v4.0
                           ${CMAKE_SOURCE_DIR}/examples++-3d/Laplace-Adapt.edp                 # depends on freeyams
                           ${CMAKE_SOURCE_DIR}/examples++-3d/Laplace-Adapt-3d.edp              # depends on mshmet
                           ${CMAKE_SOURCE_DIR}/examples++-3d/Laplace-Adapt-aniso-3d.edp        # depends on mshmet
                           ${CMAKE_SOURCE_DIR}/examples++-3d/Period-Poisson-cube-ballon.edp    # depends on freeyams
                           ${CMAKE_SOURCE_DIR}/examples++-3d/schwarz-nm-3d.edp                 # depends on metis
                           ${CMAKE_SOURCE_DIR}/examples++-3d/sphereincube.edp)                 # depends on freeyams
                    
# tests depending on freeyams library
IF(FREEFEM_FREEYAMS_INSTALLED)
  ADD_TEST(NAME Laplace-Adapt COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${CMAKE_SOURCE_DIR}/examples++-3d/Laplace-Adapt.edp)
  ADD_TEST(NAME Period-Poisson-cube-ballon COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${CMAKE_SOURCE_DIR}/examples++-3d/Period-Poisson-cube-ballon.edp)
  ADD_TEST(NAME sphereincube COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${CMAKE_SOURCE_DIR}/examples++-3d/sphereincube.edp)
ENDIF(FREEFEM_FREEYAMS_INSTALLED)

# tests depending on metis library
IF(FREEFEM_METIS_INSTALLED)
  ADD_TEST(NAME schwarz-nm-3d COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${CMAKE_SOURCE_DIR}/examples++-3d/schwarz-nm-3d.edp)
ENDIF(FREEFEM_METIS_INSTALLED)

# tests depending on mshmet library
IF(FREEFEM_MMG3D_INSTALLED)
  ADD_TEST(NAME fallingspheres COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${CMAKE_SOURCE_DIR}/examples++-3d/fallingspheres.edp)
ENDIF(FREEFEM_MMG3D_INSTALLED)

# tests depending on mshmet library
IF(FREEFEM_MSHMET_INSTALLED)
  ADD_TEST(NAME Laplace-Adapt-3d COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${CMAKE_SOURCE_DIR}/examples++-3d/Laplace-Adapt-3d.edp)
  ADD_TEST(NAME Laplace-Adapt-aniso-3d COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${CMAKE_SOURCE_DIR}/examples++-3d/Laplace-Adapt-aniso-3d.edp)
ENDIF(FREEFEM_MSHMET_INSTALLED)

MESSAGE(STATUS "tests = ${EDP_FILES}")

# tests with without dependencies
FOREACH(EDP_FILE ${EDP_FILES})
  GET_FILENAME_COMPONENT(FILENAME ${EDP_FILE} NAME_WE)
  ADD_TEST(NAME ${FILENAME} COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${EDP_FILE})
ENDFOREACH(EDP_FILE)

