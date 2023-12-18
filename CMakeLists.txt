cmake_minimum_required(VERSION 3.19.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)

set(CMAKE_VERBOSE_MAKEFILE ON)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

add_compile_definitions(
        DATA_DIR="${PROJECT_SOURCE_DIR}/simul_data"
)

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

set(PHASE2_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

link_libraries(
        QuEST m
        ${HDF5_LIBRARIES}
)

include_directories(
        ${PHASE2_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
)

if (DISTRIBUTED)
    link_libraries(
            ${MPI_C_LIBRARIES}
    )
    include_directories(
            ${MPI_C_INCLUDE_DIRS}
    )
endif ()

set(DATA_SRC_DIR ${PROJECT_SOURCE_DIR}/data)
set(DATA_SOURCES
    ${DATA_SRC_DIR}/data.c
)

set(CIRC_SRC_DIR ${PROJECT_SOURCE_DIR}/circ)
set(CIRC_SOURCES
        ${CIRC_SRC_DIR}/circ.c
        ${CIRC_SRC_DIR}/capi.c
        ${CIRC_SRC_DIR}/algos/linen.c
        ${CIRC_SRC_DIR}/algos/rayon.c
)

set(CLI_SRC_DIR ${PROJECT_SOURCE_DIR}/cli)
set(CLI_SOURCES
        ${CLI_SRC_DIR}/main.c
        ${CLI_SRC_DIR}/datio.c
        ${CLI_SRC_DIR}/log.c
        ${CLI_SRC_DIR}/opt.c
        ${CLI_SRC_DIR}/run.c
)

add_executable(ph2run
        ${DATA_SOURCES}
        ${CIRC_SOURCES}
        ${CLI_SOURCES}
)

include(CTest)
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif ()