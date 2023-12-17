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
endif ()


add_subdirectory(backend/QuEST)


add_compile_options(-Wall -Wextra -march=native)
set(CMAKE_C_FLAGS_RELEASE "-O3")

add_subdirectory(data)
add_subdirectory(circ)
add_subdirectory(ph2run)


include(CTest)
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif ()