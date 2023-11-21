message(STATUS "Configuring phase2")

set(CIRC_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

add_library(phase2 SHARED
        circ.c
        circ/linen.c
        circ/rayon.c
        sdat.c
)
target_link_libraries(phase2
        log
        QuEST m
)
target_include_directories(phase2 PRIVATE ${LOG_INCLUDE_DIRS})
