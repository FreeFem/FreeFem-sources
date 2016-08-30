#FILE(GLOB_RECURSE BAMGLIB_SRC src/bamglib/*.cpp)
LIST(REMOVE_ITEM BAMGLIB_SRC ${CMAKE_SOURCE_DIR}/src/bamglib/Meshgibbs.cpp)
ADD_EXECUTABLE(bamg ${CMAKE_SOURCE_DIR}/src/bamg/bamg.cpp 
                    ${CMAKE_SOURCE_DIR}/src/bamg/global.cpp 
                    ${BAMGLIB_SRC})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} src/bamg src/bamglib src/fflib/ src/Graphics/)
TARGET_LINK_LIBRARIES(bamg libff ${HDF5_LIBRARIES})

