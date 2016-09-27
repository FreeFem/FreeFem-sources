FILE(GLOB_RECURSE CPPFILES *.cpp)
ADD_LIBRARY(lglib ${CPPFILES})
INCLUDE_DIRECTORIES(src/bamglib src/femlib/ src/fflib/ src/Graphics/ src/lglib/)
SET_TARGET_PROPERTIES(lglib PROPERTIES PREFIX "")

