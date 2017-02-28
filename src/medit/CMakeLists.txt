FILE(WRITE ${CMAKE_SOURCE_DIR}/src/medit/compil.date "\#define COMPIL  \" jeudi 4 août 2016, 13:47:25 (UTC+0200) (with ff++ 3.47)\"")

FILE(GLOB_RECURSE CFILES *.c)
ADD_EXECUTABLE(ffmedit ${CFILES})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/libMesh 
                    ${CMAKE_SOURCE_DIR}/src/medit)
TARGET_LINK_LIBRARIES(ffmedit libMesh m ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLUT_LIBRARIES} )

