cmake_minimum_required(VERSION 3.19.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)

#set(CMAKE_VERBOSE_MAKEFILE ON)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

#set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

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
    message(STATUS "Serial HDF5 found")
endif ()


link_libraries(m)

#add_subdirectory(backend/QuEST)

add_compile_options(-Wall -Wextra -march=native)
set(CMAKE_C_FLAGS_RELEASE "-O2")

set(PH2_INCLUDE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/backend/paulirot/include/
)

set(PH2_PUBLIC_HEADERS
        ${PH2_INCLUDE_DIRS}/circ.h
        ${PH2_INCLUDE_DIRS}/data.h
        ${PH2_INCLUDE_DIRS}/algos/linen.h
        ${PH2_INCLUDE_DIRS}/algos/rayon.h
        ${PH2_INCLUDE_DIRS}/algos/silk.h
)


link_libraries(
        ${HDF5_LIBRARIES}
)
include_directories(
        ${PH2_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
)

set(PAULIROT_SRC_DIR ${PROJECT_SOURCE_DIR}/backend/paulirot/src)
set(PAULIROT_SOURCES
        ${PAULIROT_SRC_DIR}/ev.c
        ${PAULIROT_SRC_DIR}/paulis.c
        ${PAULIROT_SRC_DIR}/qreg.c
)


set(DATA2_SRC_DIR ${PROJECT_SOURCE_DIR}/data2)
set(DATA2_SOURCES
        ${DATA2_SRC_DIR}/data2.c
)

set(CIRC_SRC_DIR ${PROJECT_SOURCE_DIR}/circ)
set(CIRC_SOURCES
        ${CIRC_SRC_DIR}/circ.c
        ${CIRC_SRC_DIR}/ops.c
        ${CIRC_SRC_DIR}/hamil.c
        ${PAULIROT_SOURCES}
)

set(ALGOS_SRC_DIR ${PROJECT_SOURCE_DIR}/algos)
set(ALGOS_SOURCES
        ${ALGOS_SRC_DIR}/silk.c
)

set(PH2_COMMON_SOURCES
        ${DATA2_SOURCES}
        ${CIRC_SOURCES}
        ${ALGOS_SOURCES}
)

set(CLI_SRC_DIR ${PROJECT_SOURCE_DIR}/ph2run)
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
