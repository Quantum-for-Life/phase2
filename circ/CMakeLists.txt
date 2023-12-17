add_library(circ STATIC
        circ.c
        capi.c
        algos/linen.c
        algos/rayon.c
)

target_link_libraries(circ
        data
        QuEST m
)
