cmake_minimum_required(VERSION 3.9.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

option(LOG_USE_COLOR "enable color output for logging" ON)
if (LOG_USE_COLOR)
    message("-- Color mode for logging enabled")
    add_definitions(-DLOG_USE_COLOR)
endif ()

option(DISTRIBUTED "enable distributed mode (MPI)" OFF)
if (DISTRIBUTED)
    message("-- MPI mode enabled")
    add_definitions(-DDISTRIBUTED)

    find_package(MPI)
    include_directories(${MPI_C_INCLUDE_DIRS})

    set(HDF5_PREFER_PARALLEL TRUE)
endif ()

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

if (HDF5_IS_PARALLEL)
    message("-- Parallel HDF5 found")
else ()
    message("-- Serial HDF5 found")
endif ()

add_subdirectory(src)

add_subdirectory(lib/QuEST)
