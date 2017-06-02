cmake_minimum_required(VERSION 2.8.12)
project(RMINC)

FIND_PATH( MINC_TOOLKIT_DIR bin/mincstats  HINTS /opt/minc DOC "minc-toolkit install prefix")
FIND_PACKAGE( LIBMINC REQUIRED HINTS ${MINC_TOOLKIT_DIR}/lib${LIB_SUFFIX})
INCLUDE( ${LIBMINC_USE_FILE} )

set(CMAKE_BLAS_LIBS "-lblas" CACHE STRING "Blas libraries")
set(CMAKE_LAPACK_LIBS "-llapack" CACHE STRING "Lapack libraries")

include(GenerateExportHeader)

INCLUDE_DIRECTORIES(${LIBMINC_INCLUDE_DIRS})
link_directories(${LIBMINC_LIBRARY_DIRS})
include_directories("${CMAKE_CURRENT_BINARY_DIR}") 

add_library(RMINC SHARED 
    src/minc2_lib_reg.c
    src/minc2_read.c
    src/minc2_support.c
    src/minc2_write.c
    src/minc_anova.c
    src/minc_apply.cpp
    src/minc_cpp.cpp
    src/minc_reader.c
    src/modelling_functions.c
    src/slice_loop_functions.c
    src/vertex_modelling.c
)
target_link_libraries(RMINC ${LIBMINC_LIBRARIES} ${CMAKE_LAPACK_LIBS} ${CMAKE_BLAS_LIBS} )

generate_export_header(RMINC)

install(TARGETS RMINC DESTINATION lib)

# SET(CPPFLAGS  ${CMAKE_CXX_FLAGS} -I${NETCDF_INCLUDE_DIR} -I${HDF5_INCLUDE_DIR} -I${LIBMINC_INCLUDE_DIRS})
# SET(LDFLAGS  -L${LIBMINC_LIBRARY_DIRS})
# SET(LIBS     ${LIBMINC_LIBRARIES})

# configure flags
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/Makevars.in
#   ${CMAKE_CURRENT_BINARY_DIR}/Makevars @ONLY )
