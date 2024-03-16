if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model strict -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fconvert=big-endian -ffree-line-length-0")
endif()

add_executable(psop_ncio.x psop_ncio.f90)
target_link_libraries(psop_ncio.x PRIVATE NetCDF::NetCDF_Fortran)
target_link_libraries(psop_ncio.x PRIVATE MPI::MPI_Fortran)
target_link_libraries(psop_ncio.x PRIVATE ncio::ncio)
target_link_libraries(psop_ncio.x PRIVATE ncdiag::ncdiag)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(psop_ncio.x PRIVATE OpenMP::OpenMP_Fortran)
endif()

add_executable(psop_ncio_ensmean.x psop_ncio_ensmean.f90)
target_link_libraries(psop_ncio_ensmean.x PRIVATE NetCDF::NetCDF_Fortran)
target_link_libraries(psop_ncio_ensmean.x PRIVATE ncio::ncio)
target_link_libraries(psop_ncio_ensmean.x PRIVATE ncdiag::ncdiag)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(psop_ncio_ensmean.x PRIVATE OpenMP::OpenMP_Fortran)
endif()

# Install executable targets
install(TARGETS psop_ncio.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS psop_ncio_ensmean.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

