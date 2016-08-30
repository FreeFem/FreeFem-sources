ADD_EXECUTABLE(ffmedit ${MEDIT_SRC})
FILE(WRITE src/medit/compil.date "\#define COMPIL  \" jeudi 4 août 2016, 13:47:25 (UTC+0200) (with ff++ 3.47)\"")
INCLUDE_DIRECTORIES(src/libMesh src/medit)
TARGET_LINK_LIBRARIES(ffmedit libMesh m ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLUT_LIBRARIES} )


