message(STATUS "Configuring log")

set(LOG_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

option(LOG_USE_COLOR "enable colored output for logging" ON)
if (LOG_USE_COLOR)
    message(STATUS "Color mode logging enabled")
    add_definitions(-DLOG_USE_COLOR)
endif ()

add_library(log SHARED log.c)
