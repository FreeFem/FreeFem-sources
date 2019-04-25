macro(ff_configure_thirdparty)

  include(ff_find_package)
  include(ff_install_package)
  include(ff_configure_thirdparty_required)
  include(ff_configure_thirdparty_optional)

  ff_configure_thirdparty_required() 
  ff_configure_thirdparty_optional() 
 
endmacro(ff_configure_thirdparty)
