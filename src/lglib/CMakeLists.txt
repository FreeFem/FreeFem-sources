ADD_LIBRARY(lglib ${LGLIB_SRC})
INCLUDE_DIRECTORIES(src/bamglib src/femlib/ src/fflib/ src/Graphics/ src/lglib/)
SET_TARGET_PROPERTIES(lglib PROPERTIES PREFIX "")

