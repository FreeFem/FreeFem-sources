FILE(GLOB EDP_FILES *.edp)

IF(NOT FF_PLUGIN_MEDIT)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/Helmholtz-2d-HPDDM-BGMRES.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/elasticity-3d-thirdlevelgeneo.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/Helmholtz-3d-simple.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/Navier-2d-marmousi2.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/diffusion-2d-thirdlevelgeneo.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/Helmholtz-3d-overthrust.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/diffusion-3d-simple.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/natural_convection.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/Richards-2d.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/Helmholtz-2d-simple.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/diffusion-3d-minimal-ddm.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/natural_convection_3D_obstacle.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/diffusion-3d-minimal-direct.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/Helmholtz-2d-marmousi.edp
                             ${CMAKE_CURRENT_SOURCE_DIR}/elasticity-3d-simple.edp)
ENDIF(NOT FF_PLUGIN_MEDIT)

IF(NOT FF_PLUGIN_ELEMENT_MIXTE3D)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/Maxwell-3d-simple.edp)
ENDIF(NOT FF_PLUGIN_ELEMENT_MIXTE3D)

IF(NOT FF_PLUGIN_METIS)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/Maxwell_Cobracavity.edp)
ENDIF(NOT FF_PLUGIN_METIS)

FOREACH(EDP_FILE ${EDP_FILES})
  GET_FILENAME_COMPONENT(FILENAME ${EDP_FILE} NAME_WE)
  ADD_TEST(NAME ${FILENAME} COMMAND ${CMAKE_BINARY_DIR}/src/FreeFem++ ${EDP_FILE}
                            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ffddm)
ENDFOREACH(EDP_FILE)

