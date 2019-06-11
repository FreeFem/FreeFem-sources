macro(ff_define_ffglut_executable)

  if(GLUT_FOUND)

    # steps are the same as for FreeFEM executable
    # (see ff_define_freefem_executable for an explanation)

    add_executable(ffglut ${CMAKE_SOURCE_DIR}/src/Graphics/ffglut.cpp
                          ${CMAKE_SOURCE_DIR}/src/femlib/fem.cpp
			  ${CMAKE_SOURCE_DIR}/src/femlib/MeshSn.cpp
                          ${CMAKE_SOURCE_DIR}/src/femlib/Mesh3dn.cpp
                          ${CMAKE_SOURCE_DIR}/src/femlib/Mesh2dn.cpp
                          ${CMAKE_SOURCE_DIR}/src/femlib/Mesh1dn.cpp
                          ${CMAKE_SOURCE_DIR}/src/femlib/GQuadTree.cpp
                          ${CMAKE_SOURCE_DIR}/src/femlib/FQuadTree.cpp
                          ${CMAKE_SOURCE_DIR}/src/femlib/Drawing.cpp
                          ${CMAKE_SOURCE_DIR}/src/femlib/mshptg.cpp
                          ${CMAKE_SOURCE_DIR}/src/Graphics/gggg.cpp
                          ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp
                          ${CMAKE_SOURCE_DIR}/src/Graphics/ffthreads.cpp)
  
    include_directories(${CMAKE_SOURCE_DIR} 
                        ${CMAKE_SOURCE_DIR}/src/mpi 
                        ${CMAKE_SOURCE_DIR}/src/femlib 
                        ${CMAKE_SOURCE_DIR}/src/fflib/ 
                        ${CMAKE_SOURCE_DIR}/src/Graphics/)
  
   target_link_libraries(ffglut lglib libff
                                ${CMAKE_THREAD_LIBS_INIT} 
                                ${OPENGL_gl_LIBRARY} 
                                ${OPENGL_glu_LIBRARY} 
                                ${GLUT_LIBRARIES})
  
    install(TARGETS ffglut
            RUNTIME DESTINATION bin)

  endif()


endmacro()
