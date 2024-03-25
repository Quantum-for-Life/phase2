cmake_minimum_required(VERSION 3.19.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)

message(STATUS "Parallel MPI is required")
find_package(MPI REQUIRED)
link_libraries(${MPI_C_LIBRARIES})
include_directories(${MPI_C_INCLUDE_DIRS})

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
if (HDF5_IS_PARALLEL)
    message(STATUS "Parallel HDF5 found")
else ()
    message(STATUS "Parallel HDF5 not found")
    return()
endif ()


link_libraries(m)

add_compile_options(-Wall -Wextra -march=native)
set(CMAKE_C_FLAGS_RELEASE "-O2")

link_libraries(
        ${HDF5_LIBRARIES}
)
include_directories(
        ${HDF5_INCLUDE_DIRS}
)

set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
include_directories(${SRC_DIR})

set(CIRC_SRC_DIR ${SRC_DIR})
set(CIRC_SOURCES
        ${CIRC_SRC_DIR}/circ.c
        ${CIRC_SRC_DIR}/qreg.c
)

set(DATA2_SRC_DIR ${SRC_DIR})
set(DATA2_SOURCES
        ${DATA2_SRC_DIR}/data2.c
)

set(PH2_COMMON_SOURCES
        ${DATA2_SOURCES}
        ${CIRC_SOURCES}
)

set(CLI_SRC_DIR ${SRC_DIR})
set(CLI_SOURCES
        ${CLI_SRC_DIR}/main.c
        ${CLI_SRC_DIR}/log.c
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
