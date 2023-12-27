cmake_minimum_required(VERSION 3.19.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)

set(CMAKE_VERBOSE_MAKEFILE ON)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)


option(DISTRIBUTED "Enable distributed mode (MPI)" OFF)
if (DISTRIBUTED)
    message(STATUS "MPI mode enabled")
    find_package(MPI REQUIRED)
    add_definitions(-DDISTRIBUTED)

    set(HDF5_PREFER_PARALLEL TRUE)
endif ()

find_package(HDF5 REQUIRED)
if (HDF5_IS_PARALLEL)
    message(STATUS "Parallel HDF5 found")
else ()
    message(STATUS "Serial HDF5 found")
endif ()


add_subdirectory(backend/QuEST)

add_compile_options(-Wall -Wextra -march=native)
set(CMAKE_C_FLAGS_RELEASE "-O3")

set(PH2_PUBLIC_HEADERS
        ${PH2_INCLUDE_DIRS}/circ.h
        ${PH2_INCLUDE_DIRS}/data.h
        ${PH2_INCLUDE_DIRS}/algos/linen.h
        ${PH2_INCLUDE_DIRS}/algos/rayon.h
)

set(PH2_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
link_libraries(
        QuEST m
        ${HDF5_LIBRARIES}
)
include_directories(
        ${PH2_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
)
if (DISTRIBUTED)
    link_libraries(${MPI_C_LIBRARIES})
    include_directories(${MPI_C_INCLUDE_DIRS})
endif ()


set(DATA2_SRC_DIR ${PROJECT_SOURCE_DIR}/data2)
set(DATA2_SOURCES
        ${DATA2_SRC_DIR}/data2.c
)

set(CIRC_SRC_DIR ${PROJECT_SOURCE_DIR}/circ)
set(CIRC_SOURCES
        ${CIRC_SRC_DIR}/circ.c
        ${CIRC_SRC_DIR}/ops.c
        ${CIRC_SRC_DIR}/hamil.c
)

set(ALGOS_SRC_DIR ${PROJECT_SOURCE_DIR}/algos)
set(ALGOS_SOURCES
        ${ALGOS_SRC_DIR}/linen.c
        ${ALGOS_SRC_DIR}/rayon.c
)

set(PH2_COMMON_SOURCES
        ${DATA2_SOURCES}
        ${CIRC_SOURCES}
        ${ALGOS_SOURCES}
)

set(CLI_SRC_DIR ${PROJECT_SOURCE_DIR}/cli)
set(CLI_SOURCES
        ${CLI_SRC_DIR}/main.c
        ${CLI_SRC_DIR}/log.c
        ${CLI_SRC_DIR}/opt.c
        ${CLI_SRC_DIR}/run.c
)

add_executable(ph2run
        ${PH2_COMMON_SOURCES}
        ${CLI_SOURCES}
)

# -----------------------------------------------------------------------------
# Testing

set(PH2_SIMUL_DATA ${PROJECT_SOURCE_DIR}/simul_data)
add_compile_definitions(PH2_SIMUL_DATA="${PH2_SIMUL_DATA}")

include(CTest)
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif ()