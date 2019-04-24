macro(ff_configure_compilers)

  # Use C++11 standard
  set(CMAKE_CXX_STANDARD 11)
  
  # Set default C and C++ compiler flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DBAMG_LONG_LONG -DCMAKE -DHAVE_GETENV")


  # Set a special flag for Mac OS
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_UNISTD_H")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  # Set a special flag for MinGW
  if(MINGW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPURE_WIN32")
  endif(MINGW)

  # Manage Debug and Release modes
  if(CMAKE_BUILD_TYPE EQUAL "Debug")

    set(CMAKE_C_FLAGS "${FF_C_FLAGS} -DCHECK_KN -g" )
    set(CMAKE_CXX_FLAGS "${FF_CXX_FLAGS} -DCHECK_KN -g" )

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(CMAKE_C_FLAGS "${FF_C_FLAGS} -fno-inline -fexceptions" )
        set(CMAKE_CXX_FLAGS "${FF_CXX_FLAGS} -fno-inline -fexceptions" )

    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  else() # Release mode

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNCHECKPTR -O3" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNCHECKPTR -O3" )

  endif(CMAKE_BUILD_TYPE EQUAL "Debug")


endmacro(ff_configure_compilers)
