FILE(GLOB EDP_FILES *.edp)

IF(NOT FREEFEM_MSH3_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/3d/3d-Leman.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/EqPoisson.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/Lac.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/Laplace3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/LaplaceRT-3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/NSI3d-carac.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/NSI3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/Poisson-cube-ballon.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/Sphere-Isocahedron.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/Stokes.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/TruncLac.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/bottle.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/cone.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/convect-3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/cube-period.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/cylinder.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/cylinder-3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/extract-boundary3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/periodic-3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/pyramide.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/refinesphere.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/regtests.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/sphere2.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/sphere6.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/tetgencube.edp
                             ${CMAKE_SOURCE_DIR}/examples/3d/tetgenholeregion.edp
                             )
ENDIF(NOT FREEFEM_MSH3_INSTALLED)


IF(NOT FREEFEM_TETGEN_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/Poisson3d.edp)
ENDIF(NOT FREEFEM_TETGEN_INSTALLED)


IF(NOT FREEFEM_MEDIT_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/3d/beam-3d.edp)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/3d/intlevelset3d.edp)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/3d/meditddm.edp)
ENDIF(NOT FREEFEM_MEDIT_INSTALLED)


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

