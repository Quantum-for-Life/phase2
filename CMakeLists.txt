cmake_minimum_required(VERSION 3.9.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

option(LOG_USE_COLOR "enable color output for logging" ON)
if (LOG_USE_COLOR)
    add_definitions(-DLOG_USE_COLOR)
endif ()

option(DISTRIBUTED "enable distributed mode (MPI)" OFF)
if (DISTRIBUTED)
    add_definitions(-DDISTRIBUTED)
endif ()

add_subdirectory(lib/QuEST)

add_subdirectory(src)

