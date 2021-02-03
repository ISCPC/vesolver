SET(CMAKE_SYSTEM_NAME Linux)
#SET(CMAKE_C_COMPILER /opt/nec/ve/bin/ncc)
#SET(CMAKE_CXX_COMPILER /opt/nec/ve/bin/nc++)
#SET(CMAKE_Fortran_COMPILER /opt/nec/ve/bin/nfort)
SET(CMAKE_C_COMPILER ${SXAT_NCC_PATH}/ncc)
SET(CMAKE_CXX_COMPILER ${SXAT_NCC_PATH}/nc++)
SET(CMAKE_Fortran_COMPILER ${SXAT_NCC_PATH}/nfort)

SET(CMAKE_C_FLAGS "-O2 -g -lcblas -lblas_openmp -fopenmp"  CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-O2 -g -lcblas -lblas_openmp -fopenmp" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS "-O2 -g -llapack -lblas_openmp -fopenmp" CACHE STRING "")

SET(BLAS_LIBRARIES ${SXAT_NLC_LIBRARY_PATH}/libblas_openmp.so)
SET(LAPACK_LIBRARIES ${SXAT_NLC_LIBRARY_PATH}/liblapack.so)
