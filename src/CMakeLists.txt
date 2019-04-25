
include(ff_define_ffglut_executable)
include(ff_define_freefem_executable)
include(ff_define_freefem_mpi_executable)

# define the three executables:
# - ffglut
# - FreeFem++
# - FreeFem++-mpi

ff_define_ffglut_executable()
ff_define_freefem_executable()
ff_define_freefem_mpi_executable()


add_subdirectory(bamg)
add_subdirectory(bamglib)
add_subdirectory(fflib)
add_subdirectory(lglib)
add_subdirectory(libMesh)
add_subdirectory(medit)
add_subdirectory(mpi)

