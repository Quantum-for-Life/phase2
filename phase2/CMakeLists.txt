set(PHASE2_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PHASE2_INCLUDE_DIRS ${PHASE2_INCLUDE_DIRS} PARENT_SCOPE)

add_library(phase2 SHARED
        circ.c
        capi.c
        data.c
        algos/linen.c
        algos/rayon.c
)

target_link_libraries(phase2
        QuEST m
)

target_include_directories(phase2 PRIVATE
        ${PHASE2_INCLUDE_DIRS}
)
