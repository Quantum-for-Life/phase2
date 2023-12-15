set(CIRC_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CIRC_INCLUDE_DIRS ${CIRC_INCLUDE_DIRS} PARENT_SCOPE)

add_library(circ SHARED
        circ.c
        capi.c
        data.c
        algos/linen.c
        algos/rayon.c
)

target_link_libraries(circ
        QuEST m
)

target_include_directories(circ PRIVATE
        ${CIRC_INCLUDE_DIRS}
)
