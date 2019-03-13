macro(ff_configure_compilers)

  set(CMAKE_CXX_STANDARD 11)

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DBAMG_LONG_LONG -DCMAKE -DHAVE_GETENV")

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_UNISTD_H")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  if(MINGW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPURE_WIN32")
  endif(MINGW)

  if(CMAKE_BUILD_TYPE EQUAL "Debug")
    set(CMAKE_C_FLAGS "${FF_C_FLAGS} -DCHECK_KN -g" )
    set(CMAKE_CXX_FLAGS "${FF_CXX_FLAGS} -DCHECK_KN -g" )
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(CMAKE_C_FLAGS "${FF_C_FLAGS} -fno-inline -fexceptions" )
        set(CMAKE_CXX_FLAGS "${FF_CXX_FLAGS} -fno-inline -fexceptions" )
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNCHECKPTR -O3" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNCHECKPTR -O3" )
  endif(CMAKE_BUILD_TYPE EQUAL "Debug")


endmacro(ff_configure_compilers)
