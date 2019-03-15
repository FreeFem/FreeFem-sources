FILE(GLOB EDP_FILES *.edp)


# Tests which are known to fail are removed from test list
LIST(REMOVE_ITEM EDP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/tutorial/exception.edp
                           ${CMAKE_CURRENT_SOURCE_DIR}/tutorial/Laplace-RHS-Dirac.edp
                           ${CMAKE_CURRENT_SOURCE_DIR}/tutorial/regtests.edp)



IF(NOT FREEFEM_MEDIT_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/tutorial/Newton.edp)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/tutorial/medit.edp)
ENDIF(NOT FREEFEM_MEDIT_INSTALLED)

IF(NOT FREEFEM_MSH3_INSTALLED)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/tutorial/sphere.edp)
ENDIF(NOT FREEFEM_MSH3_INSTALLED)

FOREACH(EDP_FILE ${EDP_FILES})
  GET_FILENAME_COMPONENT(FILENAME ${EDP_FILE} NAME_WE)
  ADD_TEST(NAME ${FILENAME} COMMAND ${CMAKE_BINARY_DIR}/src/FreeFem++ ${EDP_FILE}
                            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/tutorial)
ENDFOREACH(EDP_FILE)

