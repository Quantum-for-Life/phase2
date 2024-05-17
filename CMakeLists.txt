cmake_minimum_required(VERSION 3.19.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)

add_compile_options(-Wall -Wextra -march=native)
set(CMAKE_C_FLAGS_DEBUG "-Og -g")
set(CMAKE_C_FLAGS_RELEASE "-O2")
link_libraries(m)

message(STATUS "MPI is required")
find_package(MPI REQUIRED)
link_libraries(${MPI_C_LIBRARIES})
include_directories(${MPI_C_INCLUDE_DIRS})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
if (HDF5_IS_PARALLEL)
    message(STATUS "Parallel HDF5 found")
else ()
    message(FATAL_ERROR "Parallel HDF5 not found")
endif ()
link_libraries(${HDF5_LIBRARIES})
include_directories(${HDF5_INCLUDE_DIRS})

set(PRECISION 2 CACHE STRING
        "Use single or double floating point precision in qreg. {1,2}")
if (NOT (${PRECISION} EQUAL 1) AND
        NOT (${PRECISION} EQUAL 2))
    message(FATAL_ERROR "valid options for PRECISION(=${PRECISION}) \
     are 1 for single precision, 2 for double precision")
endif ()
add_compile_definitions(QREG_PREC=${PRECISION})

set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
include_directories(${SRC_DIR})
set(PHASE2_SOURCES
        ${SRC_DIR}/data.c
        ${SRC_DIR}/qreg.c
        ${SRC_DIR}/circ.c
        ${SRC_DIR}/log.c
        ${SRC_DIR}/xoshiro256starstar.c
)

add_executable(ph2run
        ${SRC_DIR}/main.c
        ${PHASE2_SOURCES}
)

# -----------------------------------------------------------------------------
# Testing
#set(PH2_SIMUL_DATA ${PROJECT_SOURCE_DIR}/simul_data)
#add_compile_definitions(PH2_SIMUL_DATA="${PH2_SIMUL_DATA}")

#include(CTest)
#if (BUILD_TESTING)
#    enable_testing()
#    add_subdirectory(test)
#endif ()
