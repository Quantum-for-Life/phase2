cmake_minimum_required(VERSION 3.19.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)

add_compile_options(-Wall -Wextra -march=native)
set(CMAKE_C_FLAGS_RELEASE "-O3")


option(DISTRIBUTED "enable distributed mode (MPI)" OFF)
if (DISTRIBUTED)
    message(STATUS "MPI mode enabled")
    add_definitions(-DDISTRIBUTED)

    find_package(MPI)
    include_directories(${MPI_C_INCLUDE_DIRS})
    link_libraries(${MPI_C_LIBRARIES})

    set(HDF5_PREFER_PARALLEL TRUE)
endif ()

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
link_libraries(${HDF5_LIBRARIES})
if (HDF5_IS_PARALLEL)
    message(STATUS "Parallel HDF5 found")
else ()
    message(STATUS "Serial HDF5 found")
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