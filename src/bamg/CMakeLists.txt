ADD_EXECUTABLE(bamg ${CMAKE_SOURCE_DIR}/src/bamg/bamg.cpp 
                    ${CMAKE_SOURCE_DIR}/src/bamg/global.cpp) 
TARGET_LINK_LIBRARIES(bamg libff bamglib ${HDF5_LIBRARIES})


ADD_EXECUTABLE(cvmsh2 ${CMAKE_SOURCE_DIR}/src/bamg/cvmsh2.cpp 
                      ${CMAKE_SOURCE_DIR}/src/bamg/global.cpp)
TARGET_LINK_LIBRARIES(cvmsh2 bamglib libff hdf5 ${HDF5_LIBRARIES})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} 
                    ${CMAKE_SOURCE_DIR}/src/bamg 
                    ${CMAKE_SOURCE_DIR}/src/bamglib 
                    ${CMAKE_SOURCE_DIR}/src/fflib/ 
                    ${CMAKE_SOURCE_DIR}/src/Graphics/)

