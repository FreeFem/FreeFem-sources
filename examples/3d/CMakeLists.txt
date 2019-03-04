FILE(GLOB EDP_FILES *.edp)

# remove tests depending on freeyams library if not installed
IF(NOT FREEFEM_FREEYAMS_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/3d/Laplace-Adapt.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/Period-Poisson-cube-ballon.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/sphereincube.edp)
ENDIF(NOT FREEFEM_FREEYAMS_INSTALLED)

# remove tests depending on metis library if not installed
IF(NOT FREEFEM_METIS_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/3d/schwarz-nm-3d.edp)
ENDIF(NOT FREEFEM_METIS_INSTALLED)

# remove tests depending on mshmet library if not installed
IF(NOT FREEFEM_MMG3D_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/3d/fallingspheres.edp)
ENDIF(NOT FREEFEM_MMG3D_INSTALLED)

# remove tests depending on mshmet library if not installed
IF(NOT FREEFEM_MSHMET_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/3d/schwarz-nm-3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/Laplace-Adapt-3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/Laplace-Adapt-aniso-3d.edp)
ENDIF(NOT FREEFEM_MSHMET_INSTALLED)


# tests with without dependencies
FOREACH(EDP_FILE ${EDP_FILES})
  GET_FILENAME_COMPONENT(FILENAME ${EDP_FILE} NAME_WE)
  ADD_TEST(NAME ${FILENAME} COMMAND "${CMAKE_BINARY_DIR}/src/FreeFem++" ${EDP_FILE}
                            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/3d)
ENDFOREACH(EDP_FILE)

