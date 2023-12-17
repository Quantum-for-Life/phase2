set(CIRC_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CIRC_INCLUDE_DIRS ${CIRC_INCLUDE_DIRS} PARENT_SCOPE)

add_library(circ SHARED
        circ.c
        capi.c
        algos/linen.c
        algos/rayon.c
)

target_link_libraries(circ
        data
        QuEST m
)

target_include_directories(circ PRIVATE
        ${CIRC_INCLUDE_DIRS}
        ${DATA_INCLUDE_DIRS}
)
