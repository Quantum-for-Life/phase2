cmake_minimum_required(VERSION 3.19.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)

add_compile_options(-Wall -Wextra -march=native)
set(CMAKE_C_FLAGS_RELEASE "-O3")


option(DISTRIBUTED "Enable distributed mode (MPI)" OFF)
if (DISTRIBUTED)
    message(STATUS "MPI mode enabled")
    find_package(MPI REQUIRED)
    add_definitions(-DDISTRIBUTED)
endif ()


add_subdirectory(backend/QuEST)
add_subdirectory(data)
add_subdirectory(circ)
add_subdirectory(ph2run)


include(CTest)
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif ()