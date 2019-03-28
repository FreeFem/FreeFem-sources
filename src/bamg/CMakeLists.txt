if(MINGW)
  find_package(OpenMP)
  string(REPLACE "C:/msys64" "" OpenMP_C_LIBRARIES ${OpenMP_C_LIBRARIES})
endif()

ADD_EXECUTABLE(bamg ${CMAKE_SOURCE_DIR}/src/bamg/bamg.cpp 
                    ${CMAKE_SOURCE_DIR}/src/bamg/global.cpp) 

if(MINGW)
  TARGET_LINK_LIBRARIES(bamg libff bamglib ${OpenMP_C_LIBRARIES})
  target_compile_options(bamg ${OpenMP_C_FLAGS})
else()
  TARGET_LINK_LIBRARIES(bamg libff bamglib)
endif()

INSTALL(TARGETS bamg
        RUNTIME DESTINATION bin)  


ADD_EXECUTABLE(cvmsh2 ${CMAKE_SOURCE_DIR}/src/bamg/cvmsh2.cpp 
                      ${CMAKE_SOURCE_DIR}/src/bamg/global.cpp)

if(MINGW)
  TARGET_LINK_LIBRARIES(cvmsh2 bamglib libff ${OpenMP_C_LIBRARIES})
  target_compile_options(cvmsh2 ${OpenMP_C_FLAGS})
else()
  TARGET_LINK_LIBRARIES(cvmsh2 bamglib libff)
endif()

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} 
                    ${CMAKE_SOURCE_DIR}/src/bamg 
                    ${CMAKE_SOURCE_DIR}/src/bamglib 
                    ${CMAKE_SOURCE_DIR}/src/femlib 
                    ${CMAKE_SOURCE_DIR}/src/fflib 
                    ${CMAKE_SOURCE_DIR}/src/Graphics)

INSTALL(TARGETS cvmsh2
        RUNTIME DESTINATION bin)  

