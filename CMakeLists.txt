cmake_minimum_required(VERSION 3.9.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

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

message(STATUS "Configuring QuEST")
add_subdirectory(lib/QuEST)

add_subdirectory(src)

