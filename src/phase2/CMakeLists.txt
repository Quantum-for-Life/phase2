message(STATUS "Configuring phase2")

set(PHASE2_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

add_library(phase2 SHARED
        circ.c
        data.c
        algos/linen.c
        algos/rayon.c
)
target_link_libraries(phase2
        QuEST m
)
target_include_directories(phase2 PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${LOG_INCLUDE_DIRS}
)
