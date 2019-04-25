macro(ff_define_strversionnumber_library)

  include(ff_create_strversionnumber)

  ff_create_strversionnumber()
  add_library(strversionnumber STATIC strversionnumber.cpp)
  target_compile_definitions(strversionnumber PRIVATE VersionFreeFempp=${FREEFEM_VERSION})

endmacro()


